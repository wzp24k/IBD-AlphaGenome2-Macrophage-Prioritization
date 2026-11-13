from __future__ import annotations

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
TERM_PATTERN = re.compile(
    r"jak|stat|il.?17|il.?23|interleukin.?23|il6|inflammatory",
    re.IGNORECASE,
)


def targeted_terms(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty:
        return frame
    return frame[frame["Term"].astype(str).str.contains(TERM_PATTERN)].copy()


def main() -> None:
    import gseapy as gp

    configure_proxy()
    table_dir = ROOT / "results" / "tables"
    lead = pd.read_csv(table_dir / "ranked_lead_variants_active_ibd.csv")
    background = lead["l2g_target_symbol"].dropna().astype(str).unique().tolist()

    top_sets = {
        "tier1_active_supported": lead[
            (lead["evidence_tier"] == "Tier 1")
            & lead["active_disease_multicohort_supported"].fillna(False)
        ]["l2g_target_symbol"].dropna().astype(str).unique().tolist(),
        "top10": lead.head(10)["l2g_target_symbol"].dropna().astype(str).unique().tolist(),
        "top20": lead.head(20)["l2g_target_symbol"].dropna().astype(str).unique().tolist(),
        "top28": lead.head(28)["l2g_target_symbol"].dropna().astype(str).unique().tolist(),
        "active_supported_top28": lead.head(28)[
            lead.head(28)["active_disease_multicohort_supported"].fillna(False)
        ]["l2g_target_symbol"].dropna().astype(str).unique().tolist(),
    }

    all_ora = []
    for set_name, genes in top_sets.items():
        if len(genes) < 2:
            continue
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
        result["top_set"] = set_name
        result["query_size"] = len(genes)
        all_ora.append(result)
    ora = pd.concat(all_ora, ignore_index=True) if all_ora else pd.DataFrame()
    ora.to_csv(table_dir / "pathway_enrichment_full.csv", index=False)
    targeted = targeted_terms(ora)
    targeted.to_csv(table_dir / "pathway_enrichment_targeted.csv", index=False)

    rank = (
        lead.groupby("l2g_target_symbol", as_index=False)
        .agg(score=("priority_score", "max"))
        .dropna()
        .sort_values("score", ascending=False)
    )
    rank = rank[rank["l2g_target_symbol"].astype(str).str.len() > 0]
    pre = gp.prerank(
        rnk=rank,
        gene_sets=LIBRARIES,
        outdir=None,
        permutation_num=1000,
        min_size=3,
        max_size=500,
        seed=20261112,
        threads=4,
        verbose=False,
    ).res2d
    if pre is None:
        pre = pd.DataFrame()
    pre.to_csv(table_dir / "pathway_prerank_full.csv", index=False)
    pre_target = targeted_terms(pre)
    pre_target.to_csv(table_dir / "pathway_prerank_targeted.csv", index=False)

    topk_rows = []
    libraries = {
        library: gp.get_library(name=library, organism="human")
        for library in LIBRARIES
    }
    background_set = set(background)
    for k in [10, 20, 28, 50]:
        top_genes = set(
            lead.head(k)["l2g_target_symbol"].dropna().astype(str).unique()
        )
        for library, gene_sets in libraries.items():
            for term, members in gene_sets.items():
                if not TERM_PATTERN.search(str(term)):
                    continue
                members = set(members) & background_set
                if len(members) < 2:
                    continue
                hits = top_genes & members
                if not hits:
                    continue
                from scipy.stats import fisher_exact

                a = len(hits)
                b = len(top_genes - members)
                c = len((background_set - top_genes) & members)
                d = len((background_set - top_genes) - members)
                odds, p_value = fisher_exact([[a, b], [c, d]], alternative="greater")
                topk_rows.append(
                    {
                        "library": library,
                        "term": term,
                        "k": k,
                        "hit_count": a,
                        "set_size_in_background": len(members),
                        "odds_ratio": odds,
                        "p_value": p_value,
                        "genes": ",".join(sorted(hits)),
                    }
                )
    topk = pd.DataFrame(topk_rows)
    if not topk.empty:
        topk["fdr"] = np.nan
        for k, indices in topk.groupby("k").groups.items():
            p_values = topk.loc[indices, "p_value"].to_numpy()
            order = np.argsort(p_values)
            adjusted = np.empty_like(p_values)
            adjusted[order] = np.minimum.accumulate(
                (p_values[order] * len(p_values) / np.arange(1, len(p_values) + 1))[::-1]
            )[::-1]
            topk.loc[indices, "fdr"] = np.minimum(adjusted, 1)
    topk.to_csv(table_dir / "pathway_topk_targeted.csv", index=False)

    print("Targeted ORA terms:")
    print(targeted.head(50).to_string(index=False) if not targeted.empty else "None")
    print("\nTargeted prerank terms:")
    print(pre_target.head(50).to_string(index=False) if not pre_target.empty else "None")
    print("\nTop targeted top-k terms:")
    print(
        topk.sort_values(["fdr", "p_value"]).head(50).to_string(index=False)
        if not topk.empty
        else "None"
    )


if __name__ == "__main__":
    main()
