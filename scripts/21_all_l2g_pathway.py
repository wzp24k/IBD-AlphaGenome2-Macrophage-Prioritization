from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT, configure_proxy


LIBRARIES = [
    "KEGG_2021_Human",
    "Reactome_2022",
    "MSigDB_Hallmark_2020",
    "GO_Biological_Process_2023",
]
PATTERN = re.compile(
    r"jak|stat|il.?17|il.?23|interleukin.?23|il6|inflammatory",
    re.IGNORECASE,
)


def main() -> None:
    import gseapy as gp

    configure_proxy()
    table_dir = ROOT / "results" / "tables"
    raw_path = (
        ROOT
        / "data"
        / "raw"
        / "open_targets_v2_GCST90292538.json"
    )
    raw = json.loads(raw_path.read_text(encoding="utf-8"))
    l2g_rows = []
    for page in raw.get("pages", []):
        for outer in page:
            locus_id = outer.get("studyLocusId")
            for prediction in (outer.get("l2GPredictions") or {}).get("rows", []):
                target = prediction.get("target") or {}
                l2g_rows.append(
                    {
                        "studyLocusId": locus_id,
                        "l2g_score": prediction.get("score"),
                        "target_id": target.get("id"),
                        "target_symbol": target.get("approvedSymbol"),
                        "target_name": target.get("approvedName"),
                    }
                )
    l2g = pd.DataFrame(l2g_rows)
    lead = pd.read_csv(table_dir / "ranked_lead_variants_active_ibd.csv")
    lead = lead[
        [
            "studyLocusId",
            "canonical_variant",
            "l2g_target_symbol",
            "priority_score",
            "evidence_tier",
            "active_disease_multicohort_supported",
        ]
    ]
    locus_targets = lead.merge(l2g, on="studyLocusId", how="inner")
    locus_targets["weighted_score"] = (
        locus_targets["priority_score"] * locus_targets["l2g_score"]
    )
    locus_targets["target_rank"] = locus_targets.groupby("studyLocusId")[
        "l2g_score"
    ].rank(ascending=False, method="first")
    locus_targets.to_csv(table_dir / "all_l2g_target_rows.csv", index=False)

    methods = {
        "all_targets": locus_targets[
            (locus_targets["l2g_score"] >= 0.01)
            & locus_targets["target_symbol"].notna()
        ],
        "l2g_ge_0.1": locus_targets[
            (locus_targets["l2g_score"] >= 0.10)
            & locus_targets["target_symbol"].notna()
        ],
        "top3_targets": locus_targets[
            (locus_targets["target_rank"] <= 3)
            & locus_targets["target_symbol"].notna()
        ],
    }
    all_ora = []
    all_rankings = []
    for method, frame in methods.items():
        gene_scores = (
            frame.groupby(["target_id", "target_symbol"], as_index=False)
            .agg(
                weighted_score=("weighted_score", "max"),
                best_l2g=("l2g_score", "max"),
                max_priority=("priority_score", "max"),
                loci=("studyLocusId", "nunique"),
            )
            .sort_values("weighted_score", ascending=False)
        )
        gene_scores["method"] = method
        all_rankings.append(gene_scores)
        background = gene_scores["target_symbol"].astype(str).tolist()
        for k in [10, 20, 28, 50]:
            genes = gene_scores.head(k)["target_symbol"].astype(str).tolist()
            result = gp.enrichr(
                gene_list=genes,
                gene_sets=LIBRARIES,
                background=background,
                outdir=None,
                organism="human",
            ).results
            if result is None or result.empty:
                continue
            result = result.copy()
            result["method"] = method
            result["k"] = k
            all_ora.append(result)
    ora = pd.concat(all_ora, ignore_index=True)
    ora.to_csv(table_dir / "all_l2g_pathway_ora.csv", index=False)
    ora[ora["Term"].astype(str).str.contains(PATTERN)].to_csv(
        table_dir / "all_l2g_pathway_ora_targeted.csv", index=False
    )

    rankings = pd.concat(all_rankings, ignore_index=True)
    rankings.to_csv(table_dir / "all_l2g_gene_rankings.csv", index=False)
    pre_results = []
    for method, frame in rankings.groupby("method"):
        rnk = frame[["target_symbol", "weighted_score"]].copy()
        rnk.columns = ["gene", "score"]
        pre = gp.prerank(
            rnk=rnk,
            gene_sets=LIBRARIES,
            outdir=None,
            permutation_num=1000,
            min_size=3,
            max_size=500,
            seed=20261112,
            threads=4,
            verbose=False,
        ).res2d
        if pre is None or pre.empty:
            continue
        pre["method"] = method
        pre_results.append(pre)
    prerank = pd.concat(pre_results, ignore_index=True)
    prerank.to_csv(table_dir / "all_l2g_pathway_prerank.csv", index=False)
    prerank[prerank["Term"].astype(str).str.contains(PATTERN)].to_csv(
        table_dir / "all_l2g_pathway_prerank_targeted.csv", index=False
    )
    print("All-target ORA targeted, best by adjusted P:")
    cols = ["method", "k", "Term", "P-value", "Adjusted P-value", "Odds Ratio", "Genes"]
    print(
        ora[ora["Term"].astype(str).str.contains(PATTERN)]
        .sort_values(["Adjusted P-value", "P-value"])
        .head(60)[cols]
        .to_string(index=False)
    )
    print("\nAll-target prerank targeted:")
    pre_cols = ["method", "Term", "NES", "NOM p-val", "FDR q-val", "Lead_genes"]
    print(
        prerank[prerank["Term"].astype(str).str.contains(PATTERN)]
        .sort_values(["FDR q-val", "NOM p-val"])
        .head(60)[pre_cols]
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
