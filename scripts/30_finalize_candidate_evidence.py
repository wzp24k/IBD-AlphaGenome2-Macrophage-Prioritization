from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT


def main() -> None:
    table_dir = ROOT / "results" / "tables"
    lead = pd.read_csv(table_dir / "ranked_lead_variants_active_ibd.csv")
    external = pd.read_csv(table_dir / "external_gwas_gene_rankings.csv")
    external_support = external[
        ["target_symbol", "external_cluster_count", "external_score"]
    ].rename(columns={"target_symbol": "l2g_target_symbol"})
    lead = lead.merge(
        external_support,
        on="l2g_target_symbol",
        how="left",
    )
    lead["external_cluster_count"] = (
        lead["external_cluster_count"].fillna(0).astype(int)
    )
    lead["target_confidence"] = np.select(
        [
            lead["l2g_score"] >= 0.70,
            lead["l2g_score"] >= 0.50,
            lead["l2g_score"] >= 0.20,
        ],
        ["high", "moderate", "low"],
        default="very low",
    )
    lead["target_uncertain"] = lead["l2g_score"] < 0.50
    lead["final_tier_display"] = lead["final_evidence_tier"]
    lead.loc[
        lead["final_evidence_tier"].astype(str).str.startswith("Tier 1")
        & lead["target_uncertain"],
        "final_tier_display",
    ] = (
        lead.loc[
            lead["final_evidence_tier"].astype(str).str.startswith("Tier 1")
            & lead["target_uncertain"],
            "final_evidence_tier",
        ]
        + " target-uncertain"
    )
    lead = lead.sort_values(
        ["final_evidence_tier", "priority_score"],
        ascending=[True, False],
    )
    lead = lead.sort_values("priority_score", ascending=False)
    lead.to_csv(
        table_dir / "ranked_lead_variants_final.csv", index=False
    )

    ranked = pd.read_csv(
        table_dir / "ranked_variant_gene_candidates_active_ibd.csv"
    )
    target_meta = lead[
        [
            "canonical_variant",
            "external_cluster_count",
            "external_score",
            "target_confidence",
            "target_uncertain",
            "final_tier_display",
        ]
    ]
    ranked = ranked.merge(
        target_meta,
        on="canonical_variant",
        how="left",
    )
    ranked.to_csv(
        table_dir / "ranked_variant_gene_candidates_final.csv", index=False
    )
    summary = {
        "final_tier1_total": int(
            lead["final_evidence_tier"].astype(str).str.startswith("Tier 1").sum()
        ),
        "tier1_target_uncertain": int(
            (
                lead["final_evidence_tier"].astype(str).str.startswith("Tier 1")
                & lead["target_uncertain"]
            ).sum()
        ),
        "target_confidence_counts": lead["target_confidence"]
        .value_counts()
        .to_dict(),
        "tier1_external_supported": int(
            (
                lead["final_evidence_tier"].astype(str).str.startswith("Tier 1")
                & (lead["external_cluster_count"] > 0)
            ).sum()
        ),
    }
    (table_dir / "candidate_evidence_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))
    print(
        lead[
            [
                "lead_rank",
                "l2g_target_symbol",
                "l2g_score",
                "target_confidence",
                "external_cluster_count",
                "final_tier_display",
            ]
        ]
        .head(20)
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
