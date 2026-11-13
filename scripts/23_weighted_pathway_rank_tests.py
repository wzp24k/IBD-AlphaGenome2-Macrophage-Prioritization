from __future__ import annotations

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
    n = len(ranked)
    adjusted = np.empty(n, dtype=float)
    running = 1.0
    for index in range(n - 1, -1, -1):
        running = min(running, ranked[index] * n / (index + 1))
        adjusted[index] = running
    output = np.empty(n, dtype=float)
    output[order] = np.minimum(adjusted, 1)
    return output


def main() -> None:
    import gseapy as gp

    configure_proxy()
    table_dir = ROOT / "results" / "tables"
    rankings = pd.read_csv(table_dir / "all_l2g_gene_rankings.csv")
    all_targets = rankings[rankings["method"] == "all_targets"].copy()
    score_map = dict(
        zip(all_targets["target_symbol"].astype(str), all_targets["weighted_score"])
    )
    background = set(score_map)
    libraries = {
        library: gp.get_library(name=library, organism="human")
        for library in {pathway[0] for pathway in PATHWAYS}
    }
    records = []
    for library, term, label in PATHWAYS:
        members = set(libraries[library].get(term, [])) & background
        member_scores = np.array([score_map[gene] for gene in members], dtype=float)
        other_scores = np.array(
            [score_map[gene] for gene in background - members], dtype=float
        )
        statistic, p_value = mannwhitneyu(
            member_scores,
            other_scores,
            alternative="greater",
        )
        records.append(
            {
                "pathway": label,
                "library": library,
                "term": term,
                "members_in_background": len(members),
                "median_member_score": float(np.median(member_scores)),
                "median_nonmember_score": float(np.median(other_scores)),
                "mann_whitney_u": float(statistic),
                "p_value": float(p_value),
                "genes": ",".join(sorted(members)),
            }
        )
    result = pd.DataFrame(records)
    result["fdr_bh"] = bh_adjust(result["p_value"].to_numpy())
    result = result.sort_values("p_value")
    result.to_csv(table_dir / "pathway_weighted_rank_tests.csv", index=False)
    print(result.to_string(index=False))


if __name__ == "__main__":
    main()
