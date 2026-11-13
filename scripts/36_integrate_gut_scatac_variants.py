from __future__ import annotations

import json
import sys
import zipfile
from pathlib import Path

import pandas as pd
import requests
import urllib3

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT


RECORD = "21508422"
DATA_URL = f"https://zenodo.org/api/records/{RECORD}"
ZIP_NAME = "Revision_analysis.zip"
TABLE_NAME = "Revision_analysis/Data/ukbb_snp_topic_score_w_locusID.csv"
TARGETS = [
    ("PDGFB", "rs2143178", "chr22", 39264824),
    ("PRKCB", "rs7404003", "chr16", 23853326),
    ("CHSY3", "rs971891", "chr5", 130237973),
    ("ERGIC1", "rs564109", "chr5", 172897894),
    ("GCKR", "rs780094", "chr2", 27518370),
    ("PRKCQ", "rs4748198", "chr10", 6607973),
    ("ABCB5", "rs76380568", "chr7", 20542073),
    ("IRF2BP2", "rs10910476", "chr1", 234599210),
    ("IL23R", "rs1008193", "chr1", 67254228),
    ("CAMK2A", "rs1549913", "chr5", 150230409),
]


def ensure_archive(root: Path) -> Path:
    path = root / ZIP_NAME
    if path.exists() and path.stat().st_size > 0:
        return path
    urllib3.disable_warnings()
    metadata = requests.get(DATA_URL, timeout=60, verify=False).json()
    file_record = next(
        item for item in metadata["files"] if item["key"] == ZIP_NAME
    )
    response = requests.get(
        file_record["links"]["self"],
        timeout=180,
        verify=False,
    )
    response.raise_for_status()
    root.mkdir(parents=True, exist_ok=True)
    path.write_bytes(response.content)
    return path


def main() -> None:
    root = ROOT / "data" / "external" / "experimental" / "GSE273194"
    archive = ensure_archive(root)
    with zipfile.ZipFile(archive) as handle:
        with handle.open(TABLE_NAME) as stream:
            data = pd.read_csv(stream)
    rows = []
    for target, rsid, chromosome, position in TARGETS:
        match = data[
            data["GWAS_snp"].eq(rsid)
            | (
                data["Chromosome"].astype(str).eq(chromosome)
                & data["Position"].eq(position)
            )
        ]
        if match.empty:
            rows.append(
                {
                    "target": target,
                    "rsid": rsid,
                    "chromosome": chromosome,
                    "position": position,
                    "present_in_gut_scatac_variant_table": False,
                    "overlapped_cCRE": None,
                    "cCRE_celltype": None,
                    "ABC_max_gene": None,
                    "interpretation": "not_present_in_gut_scATAC_variant_table",
                }
            )
            continue
        row = match.iloc[0]
        celltype = row["overlapped_cCRE_differentially_upRegulated_celltype"]
        ccre = row["overlapped_cCRE"]
        abc_gene = row["overlapped_abcPeak_ABCmaxGene"]
        if isinstance(ccre, str) and ccre:
            interpretation = "exact_variant_in_gut_cCRE"
        else:
            interpretation = "present_without_cCRE_overlap"
        rows.append(
            {
                "target": target,
                "rsid": rsid,
                "chromosome": chromosome,
                "position": position,
                "present_in_gut_scatac_variant_table": True,
                "locus_info": row["locus_info"],
                "finemap_susie_pip": row["finemap_susie_pip"],
                "overlapped_cCRE": ccre,
                "cCRE_celltype": celltype,
                "overlapped_ABCmaxPeak": row["overlapped_ABCmaxPeak"],
                "ABC_max_gene": abc_gene,
                "overlapped_withQTL": row["overlapped_withQTL"],
                "overlapped_withQTL_inCApeak": row[
                    "overlapped_withQTL_inCApeak"
                ],
                "nearest_gene": row["snp_nearest_gene_name"],
                "interpretation": interpretation,
            }
        )
    result = pd.DataFrame(rows)
    table_dir = ROOT / "results" / "tables"
    result.to_csv(
        table_dir / "external_gut_scatac_variant_annotation.csv",
        index=False,
    )
    summary = {
        "resource": "GSE273194 / Zenodo 21508422",
        "variants_queried": len(TARGETS),
        "variants_present": int(
            result["present_in_gut_scatac_variant_table"].sum()
        ),
        "variants_in_cCRE": int(
            result["interpretation"].eq("exact_variant_in_gut_cCRE").sum()
        ),
        "results": result.to_dict("records"),
    }
    (table_dir / "external_gut_scatac_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
