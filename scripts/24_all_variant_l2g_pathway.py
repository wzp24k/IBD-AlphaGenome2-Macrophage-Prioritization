from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT, configure_proxy


PATHWAYS = [
    ("MSigDB_Hallmark_2020", "IL-6/JAK/STAT3 Signaling", "IL-6/JAK/STAT3"),
    ("MSigDB_Hallmark_2020", "IL-2/STAT5 Signaling", "IL-2/STAT5"),
    ("KEGG_2021_Human", "JAK-STAT signaling pathway", "JAK-STAT"),
    ("Reactome_2022", "Interleukin-23 Signaling R-HSA-9020933", "IL-23"),
]


def bh_adjust(p_values: np.ndarray) -> np.ndarray:
    order = np.argsort(p_values)
    ranked = p_values[order]
    adjusted = np.empty(len(ranked), dtype=float)
    running = 1.0
    for index in range(len(ranked) - 1, -1, -1):
        running = min(running, ranked[index] * len(ranked) / (index + 1))
        adjusted[index] = running
    output = np.empty(len(ranked), dtype=float)
    output[order] = np.minimum(adjusted, 1)
    return output


def pathway_rank_tests(
    rankings: pd.DataFrame,
    method: str,
    libraries: dict[str, dict[str, list[str]]],
) -> pd.DataFrame:
    frame = rankings[rankings["method"] == method].copy()
    score_map = dict(
        zip(frame["target_symbol"].astype(str), frame["gene_score"].astype(float))
    )
    background = set(score_map)
    records = []
    for library, term, label in PATHWAYS:
        members = set(libraries[library].get(term, [])) & background
        if len(members) < 2:
            continue
        member_scores = np.array([score_map[gene] for gene in members])
        other_scores = np.array(
            [score_map[gene] for gene in background - members]
        )
        _, p_value = mannwhitneyu(
            member_scores,
            other_scores,
            alternative="greater",
        )
        records.append(
            {
                "method": method,
                "pathway": label,
                "library": library,
                "term": term,
                "members_in_background": len(members),
                "median_member_score": float(np.median(member_scores)),
                "median_nonmember_score": float(np.median(other_scores)),
                "p_value": float(p_value),
                "genes": ",".join(sorted(members)),
            }
        )
    result = pd.DataFrame(records)
    if not result.empty:
        result["fdr_bh"] = bh_adjust(result["p_value"].to_numpy())
    return result


def main() -> None:
    import gseapy as gp

    configure_proxy()
    table_dir = ROOT / "results" / "tables"
    variants = pd.read_csv(
        table_dir / "ranked_variant_gene_candidates_active_ibd.csv"
    )
    variants = variants[
        [
            "studyLocusId",
            "canonical_variant",
            "priority_score",
            "active_disease_multicohort_supported",
        ]
    ].drop_duplicates("canonical_variant")
    l2g = pd.read_csv(table_dir / "all_l2g_target_rows.csv")
    l2g = l2g[
        [
            "studyLocusId",
            "l2g_score",
            "target_id",
            "target_symbol",
            "target_name",
        ]
    ].drop_duplicates(["studyLocusId", "target_id"])
    pairs = variants.merge(l2g, on="studyLocusId", how="inner")
    pairs = pairs[pairs["l2g_score"] >= 0.01].copy()
    pairs["evidence_score"] = pairs["priority_score"] * pairs["l2g_score"]

    method_frames = []
    gene_max = (
        pairs.groupby(["target_id", "target_symbol"], as_index=False)
        .agg(gene_score=("evidence_score", "max"))
        .assign(method="all_variants_max")
    )
    method_frames.append(gene_max)
    gene_sum = (
        pairs.groupby(["target_id", "target_symbol"], as_index=False)
        .agg(gene_score=("evidence_score", "sum"))
        .assign(method="all_variants_sum")
    )
    method_frames.append(gene_sum)
    top2 = (
        pairs.sort_values("evidence_score", ascending=False)
        .groupby(["target_id", "target_symbol"], as_index=False)
        .head(2)
        .groupby(["target_id", "target_symbol"], as_index=False)
        .agg(gene_score=("evidence_score", "mean"))
        .assign(method="all_variants_top2_mean")
    )
    method_frames.append(top2)
    active_pairs = pairs[
        pairs["active_disease_multicohort_supported"].fillna(False)
    ]
    active_max = (
        active_pairs.groupby(["target_id", "target_symbol"], as_index=False)
        .agg(gene_score=("evidence_score", "max"))
        .assign(method="active_supported_max")
    )
    method_frames.append(active_max)
    rankings = pd.concat(method_frames, ignore_index=True)
    all_study_rows = []
    for path in sorted((ROOT / "data" / "raw").glob("open_targets_v2_*.json")):
        payload = json.loads(path.read_text(encoding="utf-8"))
        study_id = payload["study"]["id"]
        for page in payload.get("pages", []):
            for outer in page:
                for prediction in (outer.get("l2GPredictions") or {}).get("rows", []):
                    target = prediction.get("target") or {}
                    all_study_rows.append(
                        {
                            "study": study_id,
                            "target_id": target.get("id"),
                            "target_symbol": target.get("approvedSymbol"),
                        }
                    )
    study_counts = (
        pd.DataFrame(all_study_rows)
        .drop_duplicates(["study", "target_id"])
        .groupby(["target_id", "target_symbol"], as_index=False)
        .agg(independent_study_count=("study", "nunique"))
    )
    study_counts.to_csv(
        table_dir / "cross_study_l2g_target_counts.csv", index=False
    )
    cross_study = gene_max.drop(columns=["method"]).merge(
        study_counts,
        on=["target_id", "target_symbol"],
        how="left",
    )
    cross_study["independent_study_count"] = cross_study[
        "independent_study_count"
    ].fillna(1)
    cross_study["gene_score"] = (
        cross_study["gene_score"] * cross_study["independent_study_count"]
    )
    cross_study["method"] = "cross_study_supported_max"
    rankings = pd.concat([rankings, cross_study], ignore_index=True)
    cross_study_log = gene_max.drop(columns=["method"]).merge(
        study_counts,
        on=["target_id", "target_symbol"],
        how="left",
    )
    cross_study_log["independent_study_count"] = cross_study_log[
        "independent_study_count"
    ].fillna(1)
    cross_study_log["gene_score"] = (
        cross_study_log["gene_score"]
        * np.log2(1 + cross_study_log["independent_study_count"])
    )
    cross_study_log["method"] = "cross_study_log_supported"
    rankings = pd.concat([rankings, cross_study_log], ignore_index=True)
    cross_study_sqrt = gene_max.drop(columns=["method"]).merge(
        study_counts,
        on=["target_id", "target_symbol"],
        how="left",
    )
    cross_study_sqrt["independent_study_count"] = cross_study_sqrt[
        "independent_study_count"
    ].fillna(1)
    cross_study_sqrt["gene_score"] = (
        cross_study_sqrt["gene_score"]
        * np.sqrt(cross_study_sqrt["independent_study_count"])
    )
    cross_study_sqrt["method"] = "cross_study_sqrt_supported"
    rankings = pd.concat([rankings, cross_study_sqrt], ignore_index=True)
    rankings.to_csv(table_dir / "all_variant_l2g_gene_rankings.csv", index=False)
    pairs.to_csv(table_dir / "all_variant_l2g_pairs.csv", index=False)

    libraries = {
        library: gp.get_library(name=library, organism="human")
        for library in {pathway[0] for pathway in PATHWAYS}
    }
    results = pd.concat(
        [
            pathway_rank_tests(rankings, method, libraries)
            for method in rankings["method"].unique()
        ],
        ignore_index=True,
    )
    results = results.sort_values(["method", "p_value"])
    results.to_csv(table_dir / "all_variant_pathway_rank_tests.csv", index=False)
    summary = {
        method: {
            row["pathway"]: {
                "members": int(row["members_in_background"]),
                "p": float(row["p_value"]),
                "fdr": float(row["fdr_bh"]),
            }
            for _, row in results[results["method"] == method].iterrows()
        }
        for method in results["method"].unique()
    }
    (table_dir / "all_variant_pathway_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    print(results.to_string(index=False))


if __name__ == "__main__":
    main()
