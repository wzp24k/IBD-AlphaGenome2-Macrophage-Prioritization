from __future__ import annotations

import json
import sys
from pathlib import Path

import gseapy as gp
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT, configure_proxy


TARGETS = ["PRKCB", "CHSY3"]


def main() -> None:
    configure_proxy()
    path = (
        ROOT
        / "data"
        / "external"
        / "experimental"
        / "GSE175685"
        / "GSE175685_THP1_overexpression_processed_20210506.txt.gz"
    )
    frame = pd.read_csv(path, sep="\t").set_index("Gene.Symbol")
    columns = list(frame.columns)
    empty = [column for column in columns if column.startswith("Empty_")]
    libraries = {
        library: gp.get_library(library, organism="human")
        for library in [
            "MSigDB_Hallmark_2020",
            "KEGG_2021_Human",
            "Reactome_2022",
        ]
    }
    table_dir = ROOT / "results" / "tables"
    all_results = []
    summary = {}
    for target in TARGETS:
        target_columns = [
            column for column in columns if column.startswith(target + "_")
        ]
        target_values = frame[target_columns].to_numpy(dtype=float)
        empty_values = frame[empty].to_numpy(dtype=float)
        statistic, p_value = ttest_ind(
            target_values,
            empty_values,
            axis=1,
            equal_var=False,
            nan_policy="omit",
        )
        ranking = pd.DataFrame(
            {
                "gene": frame.index,
                "stat": statistic,
                "pvalue": p_value,
                "delta": target_values.mean(axis=1) - empty_values.mean(axis=1),
            }
        ).dropna()
        ranking = (
            ranking.sort_values("stat", ascending=False)
            .drop_duplicates("gene")
            .set_index("gene")
        )
        ranking.index = ranking.index.astype(str).astype(object)
        rank_metric = ranking["stat"].astype(float)
        rank_metric.index = rank_metric.index.astype(object)
        for library, gene_sets in libraries.items():
            result = gp.prerank(
                rnk=rank_metric,
                gene_sets=gene_sets,
                min_size=5,
                max_size=500,
                permutation_num=1000,
                threads=4,
                seed=20261113,
                outdir=None,
            ).res2d
            result["target"] = target
            result["library"] = library
            all_results.append(result)
            summary[target] = summary.get(target, {})
            top = result.sort_values("FDR q-val").head(10).to_dict("records")
            summary[target][library] = top
    result = pd.concat(all_results, ignore_index=True)
    result.to_csv(
        table_dir / "gse175685_prkcb_chsy3_gsea.csv",
        index=False,
    )
    (table_dir / "gse175685_prkcb_chsy3_gsea_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    print(
        result.sort_values(["target", "FDR q-val"])
        .groupby("target")
        .head(12)
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
