from __future__ import annotations

import json
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT


MARKERS = {
    "Macrophage": [
        "CD68",
        "CD14",
        "LYZ",
        "C1QA",
        "C1QB",
        "C1QC",
        "APOE",
        "FCGR3A",
        "MRC1",
        "MSR1",
        "SPP1",
    ],
    "DC": ["CD1C", "FCER1A", "CLEC10A", "CD1E", "IRF8"],
    "T_cell": ["CD3D", "CD3E", "CD2", "IL7R", "LTB"],
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "MUC2", "PIGR"],
    "Inflammation": [
        "TNF",
        "IL1B",
        "CXCL8",
        "NFKB1",
        "NFKBIA",
        "CCL20",
        "SAA1",
    ],
}


def main() -> None:
    zip_path = (
        ROOT
        / "data"
        / "external"
        / "experimental"
        / "GSE273194"
        / "Revision_analysis.zip"
    )
    with zipfile.ZipFile(zip_path) as archive:
        with archive.open(
            "Revision_analysis/Data/ukbb_snp_topic_score_w_locusID.csv"
        ) as stream:
            variants = pd.read_csv(stream)
        with archive.open(
            "Revision_analysis/Data/Topic_genescores_all_proteinCoding_Clean.txt"
        ) as stream:
            genes = pd.read_csv(stream, sep=None, engine="python")
    variant = variants[variants["GWAS_snp"].eq("rs2143178")].iloc[0]
    topic_columns = [column for column in variants if column.startswith("topic")]
    percentile_rows = []
    for topic in topic_columns:
        values = variants[topic].dropna()
        value = variant[topic]
        percentile = float((values <= value).mean()) if pd.notna(value) else np.nan
        percentile_rows.append(
            {
                "topic": topic.replace("topic", "Topic"),
                "score": value if pd.notna(value) else np.nan,
                "percentile": percentile,
            }
        )
    percentiles = pd.DataFrame(percentile_rows).sort_values(
        "percentile", ascending=False
    )
    gene_scores = genes.set_index("Gene")
    marker_rows = []
    for topic in topic_columns:
        column = topic.replace("topic", "Topic")
        if column not in gene_scores:
            continue
        for marker_set, marker_genes in MARKERS.items():
            present = [gene for gene in marker_genes if gene in gene_scores.index]
            score = float(gene_scores.loc[present, column].mean())
            marker_rows.append(
                {
                    "topic": column,
                    "marker_set": marker_set,
                    "mean_gene_score": score,
                    "genes_used": len(present),
                }
            )
    markers = pd.DataFrame(marker_rows)
    markers["within_set_rank"] = markers.groupby("marker_set")[
        "mean_gene_score"
    ].rank(ascending=False, method="min")
    table_dir = ROOT / "results" / "tables"
    percentiles.to_csv(
        table_dir / "scatac_pdgfb_variant_topic_percentiles.csv",
        index=False,
    )
    markers.to_csv(
        table_dir / "scatac_topic_marker_enrichment.csv",
        index=False,
    )
    summary = {
        "variant": "rs2143178",
        "cCRE": variant["overlapped_cCRE"],
        "celltype": variant[
            "overlapped_cCRE_differentially_upRegulated_celltype"
        ],
        "top_variant_topics": percentiles.head(10).to_dict("records"),
        "top_macrophage_topics": markers[
            markers["marker_set"].eq("Macrophage")
        ].sort_values("within_set_rank").head(10).to_dict("records"),
        "top_dc_topics": markers[
            markers["marker_set"].eq("DC")
        ].sort_values("within_set_rank").head(10).to_dict("records"),
    }
    (table_dir / "scatac_pdgfb_topic_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
