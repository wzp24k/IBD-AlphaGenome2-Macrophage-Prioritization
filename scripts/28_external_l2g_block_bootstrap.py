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
BOOTSTRAPS = 2000


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


def auc_difference(member_scores: np.ndarray, other_scores: np.ndarray) -> float:
    if len(member_scores) == 0 or len(other_scores) == 0:
        return np.nan
    statistic, _ = mannwhitneyu(
        member_scores,
        other_scores,
        alternative="two-sided",
    )
    return float(statistic / (len(member_scores) * len(other_scores)) - 0.5)


def main() -> None:
    import gseapy as gp

    configure_proxy()
    table_dir = ROOT / "results" / "tables"
    pairs = pd.read_csv(table_dir / "all_variant_l2g_pairs.csv")
    external = pd.read_csv(table_dir / "external_ibd_l2g_predictions.csv")
    external_counts = (
        external.dropna(subset=["target_id"])
        .drop_duplicates(["cluster", "target_id"])
        .groupby(["target_id", "target_symbol"], as_index=False)
        .agg(external_cluster_count=("cluster", "nunique"))
    )
    locus_gene = (
        pairs.groupby(
            ["studyLocusId", "target_id", "target_symbol"],
            as_index=False,
        )
        .agg(primary_evidence=("evidence_score", "max"))
    )
    gene_primary = (
        locus_gene.groupby(["target_id", "target_symbol"], as_index=False)
        .agg(primary_score=("primary_evidence", "max"))
    )
    genes = gene_primary.merge(
        external_counts,
        on=["target_id", "target_symbol"],
        how="left",
    )
    genes["external_cluster_count"] = genes["external_cluster_count"].fillna(0)
    genes["external_factor"] = 1 + np.log2(1 + genes["external_cluster_count"])
    genes["external_score"] = genes["primary_score"] * genes["external_factor"]
    genes.to_csv(
        table_dir / "external_gwas_gene_rankings.csv", index=False
    )

    locus_matrix = locus_gene.pivot_table(
        index="studyLocusId",
        columns="target_id",
        values="primary_evidence",
        aggfunc="max",
        fill_value=0,
    )
    factor_by_target = dict(
        zip(genes["target_id"].astype(str), genes["external_factor"])
    )
    symbol_by_target = dict(
        zip(genes["target_id"].astype(str), genes["target_symbol"])
    )
    target_by_symbol = {
        symbol: target
        for target, symbol in symbol_by_target.items()
        if isinstance(symbol, str)
    }
    libraries = {
        library: gp.get_library(library, organism="human")
        for library in {pathway[0] for pathway in PATHWAYS}
    }
    locus_ids = locus_matrix.index.to_numpy()
    rng = np.random.default_rng(20261112)

    bootstrap_effects = {label: [] for _, _, label in PATHWAYS}
    for _ in range(BOOTSTRAPS):
        sampled_loci = rng.choice(locus_ids, size=len(locus_ids), replace=True)
        sampled_matrix = locus_matrix.loc[sampled_loci]
        primary_scores = sampled_matrix.max(axis=0)
        present = primary_scores[primary_scores > 0]
        score_by_symbol = {}
        for target, value in present.items():
            symbol = symbol_by_target.get(str(target))
            if isinstance(symbol, str):
                score_by_symbol[symbol] = float(value) * factor_by_target[str(target)]
        background_symbols = set(score_by_symbol)
        values = np.array(list(score_by_symbol.values()))
        for library, term, label in PATHWAYS:
            members = set(libraries[library].get(term, [])) & background_symbols
            if len(members) < 2:
                bootstrap_effects[label].append(np.nan)
                continue
            member_values = np.array([score_by_symbol[gene] for gene in members])
            other_values = np.array(
                [score_by_symbol[gene] for gene in background_symbols - members]
            )
            bootstrap_effects[label].append(
                auc_difference(member_values, other_values)
            )

    observed_scores = dict(zip(genes["target_symbol"], genes["external_score"]))
    observed_background = {
        symbol: score
        for symbol, score in observed_scores.items()
        if isinstance(symbol, str) and np.isfinite(score) and score > 0
    }
    records = []
    for library, term, label in PATHWAYS:
        members = set(libraries[library].get(term, [])) & set(observed_background)
        member_values = np.array([observed_background[gene] for gene in members])
        other_values = np.array(
            [
                score
                for gene, score in observed_background.items()
                if gene not in members
            ]
        )
        effect = auc_difference(member_values, other_values)
        distribution = np.array(bootstrap_effects[label], dtype=float)
        distribution = distribution[np.isfinite(distribution)]
        p_value = (1 + np.sum(distribution <= 0)) / (1 + len(distribution))
        records.append(
            {
                "pathway": label,
                "library": library,
                "term": term,
                "members": len(members),
                "external_supported_members": int(
                    sum(
                        int(
                            genes.loc[
                                genes["target_symbol"].eq(gene),
                                "external_cluster_count",
                            ].max()
                        )
                        > 0
                        for gene in members
                    )
                ),
                "observed_auc_difference": effect,
                "bootstrap_ci_low": float(np.quantile(distribution, 0.025)),
                "bootstrap_ci_high": float(np.quantile(distribution, 0.975)),
                "bootstrap_p": p_value,
            }
        )
    result = pd.DataFrame(records)
    result["bootstrap_fdr"] = bh_adjust(result["bootstrap_p"].to_numpy())
    result = result.sort_values("bootstrap_p")
    result.to_csv(
        table_dir / "external_gwas_block_bootstrap_pathways.csv", index=False
    )
    summary = {
        "bootstrap_replicates": BOOTSTRAPS,
        "loci": int(len(locus_ids)),
        "genes": int(len(genes)),
        "externally_supported_genes": int(
            (genes["external_cluster_count"] > 0).sum()
        ),
        "external_clusters": ["UKB_WGS_2025", "MVP_2024"],
        "pathways": result.to_dict("records"),
    }
    (table_dir / "external_gwas_block_bootstrap_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    print(result.to_string(index=False))


if __name__ == "__main__":
    main()
