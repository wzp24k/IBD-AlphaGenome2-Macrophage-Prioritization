from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT, load_config


OUTPUT_NAMES = {
    "RNA_SEQ": "rna",
    "SC_GENE_EXPRESSION": "sc",
    "SC_GENE_EXPRESSION_DENOISED": "custom",
}


def percentile(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    return values.rank(method="average", pct=True)


def aggregate_condition_effects(scores: pd.DataFrame) -> pd.DataFrame:
    usable = scores[
        scores["imputation_mode"].isin(["observed", "imputed_supported"])
    ].copy()
    grouped = (
        usable.groupby(
            [
                "variant_id",
                "gene_id",
                "gene_name",
                "output_type",
                "condition_role",
            ],
            as_index=False,
        )
        .agg(raw_score=("raw_score", "mean"), track_count=("raw_score", "size"))
    )
    return grouped


def effect_matrix(grouped: pd.DataFrame, output_type: str) -> pd.DataFrame:
    subset = grouped[grouped["output_type"] == output_type]
    if subset.empty:
        return pd.DataFrame()
    pivot = subset.pivot_table(
        index=["variant_id", "gene_id", "gene_name"],
        columns="condition_role",
        values="raw_score",
        aggfunc="mean",
    ).reset_index()
    if "primary_target" not in pivot.columns:
        return pd.DataFrame()
    control_columns = [
        column
        for column in ["immune_control", "epithelial_control", "epithelial_sensitivity"]
        if column in pivot.columns
    ]
    pivot["macrophage_abs"] = pivot.get("primary_target", np.nan).abs()
    pivot["control_abs_max"] = pivot[control_columns].abs().max(axis=1)
    pivot["macrophage_abs_percentile"] = pivot["macrophage_abs"].rank(pct=True)
    pivot["control_abs_percentile"] = pivot["control_abs_max"].rank(pct=True)
    pivot["specificity"] = (
        pivot["macrophage_abs_percentile"] - pivot["control_abs_percentile"]
    )
    return pivot


def merge_target_effect(
    matrix: pd.DataFrame,
    candidates: pd.DataFrame,
    output_type: str,
) -> pd.DataFrame:
    target = matrix.merge(
        candidates[["canonical_variant", "l2g_target_symbol"]],
        left_on=["variant_id", "gene_name"],
        right_on=["canonical_variant", "l2g_target_symbol"],
        how="inner",
    )
    prefix = OUTPUT_NAMES[output_type]
    keep = target[
        [
            "variant_id",
            "gene_name",
            "l2g_target_symbol",
            "primary_target",
            "macrophage_abs",
            "macrophage_abs_percentile",
            "control_abs_max",
            "control_abs_percentile",
            "specificity",
        ]
    ].rename(
        columns={
            "variant_id": "canonical_variant",
            "gene_name": f"{prefix}_target_gene",
            "primary_target": f"{prefix}_macrophage_effect",
            "macrophage_abs": f"{prefix}_macrophage_abs_effect",
            "macrophage_abs_percentile": f"{prefix}_macrophage_abs_percentile",
            "control_abs_max": f"{prefix}_control_abs_max",
            "control_abs_percentile": f"{prefix}_control_abs_percentile",
            "specificity": f"{prefix}_specificity",
        }
    )
    return keep.drop(columns=[f"{prefix}_target_gene"])


def background_calibration(
    score_path: Path,
    candidate_effects: pd.DataFrame,
) -> pd.DataFrame:
    if not score_path.exists():
        return pd.DataFrame(
            {
                "canonical_variant": candidate_effects["canonical_variant"],
                "rna_background_macrophage_p": np.nan,
                "sc_background_macrophage_p": np.nan,
                "rna_background_specificity_p": np.nan,
                "sc_background_specificity_p": np.nan,
            }
        )
    scores = pd.read_parquet(score_path)
    grouped = aggregate_condition_effects(scores)
    result = candidate_effects[["canonical_variant"]].copy()
    for output_type, prefix in (("RNA_SEQ", "rna"), ("SC_GENE_EXPRESSION", "sc")):
        matrix = effect_matrix(grouped, output_type)
        if matrix.empty:
            continue
        mac_null = matrix["macrophage_abs"].dropna().to_numpy()
        spec_null = matrix["specificity"].dropna().to_numpy()
        candidate_column = f"{prefix}_macrophage_abs_effect"
        specificity_column = f"{prefix}_specificity"
        if candidate_column not in candidate_effects:
            continue
        p_mac = []
        p_spec = []
        for value, specificity in zip(
            candidate_effects[candidate_column],
            candidate_effects[specificity_column],
        ):
            if pd.isna(value) or len(mac_null) == 0:
                p_mac.append(np.nan)
            else:
                p_mac.append((1 + np.sum(mac_null >= value)) / (1 + len(mac_null)))
            if pd.isna(specificity) or len(spec_null) == 0:
                p_spec.append(np.nan)
            else:
                p_spec.append(
                    (1 + np.sum(spec_null >= specificity)) / (1 + len(spec_null))
                )
        result[f"{prefix}_background_macrophage_p"] = p_mac
        result[f"{prefix}_background_specificity_p"] = p_spec
    return result


def pathway_enrichment(lead: pd.DataFrame) -> pd.DataFrame:
    gene_sets = {
        "IL23_Th17": {
            "IL23R",
            "IL12B",
            "JAK2",
            "TYK2",
            "STAT3",
            "STAT4",
            "RORC",
            "IL17A",
            "IL17RA",
            "IL17RC",
            "IL21",
            "IL22",
        },
        "JAK_STAT": {
            "JAK1",
            "JAK2",
            "JAK3",
            "TYK2",
            "STAT1",
            "STAT2",
            "STAT3",
            "STAT4",
            "STAT5A",
            "STAT5B",
            "STAT6",
        },
        "Macrophage_innate": {
            "CCL20",
            "LSP1",
            "TNFSF15",
            "TNFRSF6B",
            "PLCG2",
            "NOD2",
            "CARD9",
            "LY96",
            "IRF5",
            "IRF8",
            "FCGR2A",
            "FCGR2B",
        },
    }
    all_genes = set(lead["l2g_target_symbol"].dropna())
    top_count = min(20, len(lead))
    top_genes = set(lead.head(top_count)["l2g_target_symbol"].dropna())
    records = []
    for name, gene_set in gene_sets.items():
        a = len(top_genes & gene_set)
        b = len(top_genes - gene_set)
        c = len((all_genes - top_genes) & gene_set)
        d = len((all_genes - top_genes) - gene_set)
        odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative="greater")
        records.append(
            {
                "gene_set": name,
                "top_genes_in_set": a,
                "top_genes": top_count,
                "background_genes_in_set": c,
                "background_genes": len(all_genes - top_genes),
                "odds_ratio": odds_ratio,
                "p_value": p_value,
                "top_gene_symbols": ",".join(sorted(top_genes & gene_set)),
            }
        )
    return pd.DataFrame(records)


def main() -> None:
    config = load_config()
    table_dir = ROOT / "results" / "tables"
    all_candidates = pd.read_csv(
        ROOT / "data" / "processed" / "ibd_noncoding_credible_variants.csv"
    )
    primary = (
        all_candidates[all_candidates["study_role"] == "primary_ibd"]
        .drop_duplicates("canonical_variant")
        .copy()
    )
    replication = (
        all_candidates[all_candidates["study_role"] != "primary_ibd"]
        .groupby("canonical_variant", as_index=False)
        .agg(
            replication_study_count=("studyId", "nunique"),
            replication_roles=(
                "study_role",
                lambda values: ",".join(sorted(set(values))),
            ),
        )
    )
    primary = primary.merge(replication, on="canonical_variant", how="left")
    primary["replication_study_count"] = (
        primary["replication_study_count"].fillna(0).astype(int)
    )

    atlas = pd.read_csv(table_dir / "atlas_rna_seq_gene_quantiles.csv")
    atlas_target = atlas.merge(
        primary[["canonical_variant", "l2g_target_symbol"]],
        left_on=["variant_id", "gene_name"],
        right_on=["canonical_variant", "l2g_target_symbol"],
        how="inner",
    )
    atlas_target["atlas_abs_q_target"] = atlas_target["q_at_max_abs_raw"].abs()
    atlas_target = atlas_target[
        [
            "variant_id",
            "gene_name",
            "l2g_target_symbol",
            "raw_at_max_abs_raw",
            "q_at_max_abs_raw",
            "atlas_abs_q_target",
            "biosample_at_max_abs_raw",
            "max_abs_q",
            "q_at_max_abs_q",
            "biosample_at_max_abs_q",
            "colon_intestine_max_abs_q",
            "colon_intestine_q_at_max_abs_q",
            "colon_intestine_biosample_at_max_abs_q",
        ]
    ].rename(
        columns={
            "variant_id": "canonical_variant",
            "gene_name": "atlas_target_gene",
            "raw_at_max_abs_raw": "atlas_raw_at_target",
            "q_at_max_abs_raw": "atlas_q_at_target",
            "max_abs_q": "atlas_global_max_abs_q",
            "q_at_max_abs_q": "atlas_global_q_at_max_abs_q",
            "biosample_at_max_abs_q": "atlas_global_track",
            "colon_intestine_max_abs_q": "atlas_gi_max_abs_q",
            "colon_intestine_q_at_max_abs_q": "atlas_gi_q",
            "colon_intestine_biosample_at_max_abs_q": "atlas_gi_track",
        }
    )
    ranked = primary.merge(
        atlas_target, on=["canonical_variant", "l2g_target_symbol"], how="left"
    )

    live = pd.read_parquet(table_dir / "live_condition_gene_scores.parquet")
    live = live[~live["variant_id"].astype(str).str.contains("chr21:39094644", regex=False)]
    grouped = aggregate_condition_effects(live)
    rna = effect_matrix(grouped, "RNA_SEQ")
    sc = effect_matrix(grouped, "SC_GENE_EXPRESSION")
    rna_target = merge_target_effect(rna, primary, "RNA_SEQ")
    sc_target = merge_target_effect(sc, primary, "SC_GENE_EXPRESSION")
    ranked = ranked.merge(rna_target, on="canonical_variant", how="left")
    ranked = ranked.merge(sc_target, on="canonical_variant", how="left")

    custom = pd.read_parquet(
        table_dir / "custom_colorectal_macrophage_gene_scores.parquet"
    )
    custom_grouped = (
        custom.groupby(["variant_id", "gene_name"], as_index=False)
        .agg(raw_score=("raw_score", "mean"))
    )
    custom_grouped["custom_macrophage_abs"] = custom_grouped["raw_score"].abs()
    custom_grouped["custom_macrophage_abs_percentile"] = custom_grouped[
        "custom_macrophage_abs"
    ].rank(pct=True)
    custom_target = custom_grouped.merge(
        primary[["canonical_variant", "l2g_target_symbol"]],
        left_on=["variant_id", "gene_name"],
        right_on=["canonical_variant", "l2g_target_symbol"],
        how="inner",
    ).rename(
        columns={
            "raw_score": "custom_macrophage_effect",
            "custom_macrophage_abs": "custom_macrophage_abs_effect",
            "custom_macrophage_abs_percentile": "custom_macrophage_abs_percentile",
        }
    )
    ranked = ranked.merge(
        custom_target[
            [
                "canonical_variant",
                "custom_macrophage_effect",
                "custom_macrophage_abs_effect",
                "custom_macrophage_abs_percentile",
            ]
        ],
        on="canonical_variant",
        how="left",
    )

    ranked["live_specificity"] = ranked[
        ["rna_specificity", "sc_specificity"]
    ].mean(axis=1, skipna=True)
    ranked["live_macrophage_abs_percentile"] = ranked[
        ["rna_macrophage_abs_percentile", "sc_macrophage_abs_percentile"]
    ].mean(axis=1, skipna=True)
    ranked = ranked.merge(
        background_calibration(
            table_dir / "background_live_condition_gene_scores.parquet",
            ranked,
        ),
        on="canonical_variant",
        how="left",
    )

    percentile_columns = {
        "pip": "pip_percentile",
        "l2g_score": "l2g_percentile",
        "atlas_gi_max_abs_q": "atlas_gi_percentile",
        "live_specificity": "live_specificity_percentile",
        "custom_macrophage_abs_percentile": "custom_macrophage_percentile",
    }
    for source, target in percentile_columns.items():
        ranked[target] = percentile(ranked[source])
    percentile_matrix = ranked[list(percentile_columns.values())]
    ranked["model_evidence_axes"] = ranked[
        [
            "atlas_gi_max_abs_q",
            "rna_specificity",
            "sc_specificity",
            "custom_macrophage_abs_effect",
        ]
    ].notna().sum(axis=1)
    ranked["primary_rank_eligible"] = (
        ranked["l2g_target_symbol"].notna() & (ranked["model_evidence_axes"] >= 2)
    )
    ranked["priority_score"] = 100 * percentile_matrix.sum(
        axis=1, skipna=True
    ) / len(percentile_columns)
    ranked["genetics_weighted_score"] = 100 * (
        0.30 * ranked["pip_percentile"].fillna(0)
        + 0.30 * ranked["l2g_percentile"].fillna(0)
        + 0.20 * ranked["atlas_gi_percentile"].fillna(0)
        + 0.20 * ranked["live_specificity_percentile"].fillna(0)
    )
    ranked["macrophage_weighted_score"] = 100 * (
        0.15 * ranked["pip_percentile"].fillna(0)
        + 0.15 * ranked["l2g_percentile"].fillna(0)
        + 0.20 * ranked["atlas_gi_percentile"].fillna(0)
        + 0.30 * ranked["live_specificity_percentile"].fillna(0)
        + 0.20 * ranked["custom_macrophage_percentile"].fillna(0)
    )
    ranked["priority_rank"] = ranked["priority_score"].rank(
        method="min", ascending=False
    ).astype(int)

    ranked["atlas_supported"] = ranked["atlas_gi_max_abs_q"] >= float(
        config["thresholds"]["macrophage_quantile_high"]
    )
    specificity_threshold = float(config["thresholds"]["specificity_delta_min"])
    ranked["macrophage_selective"] = (
        (ranked["rna_specificity"] >= specificity_threshold)
        | (ranked["sc_specificity"] >= specificity_threshold)
    )
    ranked["patient_macrophage_supported"] = (
        ranked["custom_macrophage_abs_percentile"] >= 0.8
    )
    evidence_count = ranked[
        ["atlas_supported", "macrophage_selective", "patient_macrophage_supported"]
    ].sum(axis=1)
    ranked["evidence_tier"] = np.select(
        [evidence_count == 3, evidence_count == 2, evidence_count == 1],
        ["Tier 1", "Tier 2", "Tier 3"],
        default="Tier 4",
    )
    ranked["evidence_tier_order"] = ranked["evidence_tier"].str.extract(r"(\d)").astype(int)
    ranked = ranked.sort_values("priority_rank")
    lead = (
        ranked.sort_values(
            ["studyLocusId", "evidence_tier_order", "priority_score"],
            ascending=[True, True, False],
        )
        .groupby("studyLocusId", as_index=False)
        .head(1)
        .sort_values(
            ["evidence_tier_order", "priority_score"], ascending=[True, False]
        )
        .copy()
    )
    lead["lead_rank"] = np.arange(1, len(lead) + 1)

    corpus = []
    for column in ["priority_score", "genetics_weighted_score", "macrophage_weighted_score"]:
        rho, p_value = spearmanr(
            lead[column], lead["atlas_gi_max_abs_q"], nan_policy="omit"
        )
        corpus.append(
            {
                "comparison": f"{column}_vs_atlas_gi_q",
                "spearman_rho": rho,
                "p_value": p_value,
                "n": int(lead[[column, "atlas_gi_max_abs_q"]].dropna().shape[0]),
            }
        )
    correlations = pd.DataFrame(corpus)
    enrichment = pathway_enrichment(lead)

    ranked.to_csv(table_dir / "ranked_variant_gene_candidates.csv", index=False)
    lead.to_csv(table_dir / "ranked_lead_variants.csv", index=False)
    legacy_dir = ROOT / "results" / "legacy"
    legacy_dir.mkdir(parents=True, exist_ok=True)
    enrichment.to_csv(legacy_dir / "legacy_top20_pathway_fisher.csv", index=False)
    correlations.to_csv(table_dir / "ranking_sensitivity_correlations.csv", index=False)
    summary = {
        "primary_variants": int(len(primary)),
        "variants_with_atlas_target": int(ranked["atlas_gi_max_abs_q"].notna().sum()),
        "variants_with_live_target": int(ranked["live_specificity"].notna().sum()),
        "variants_with_custom_macrophage": int(
            ranked["custom_macrophage_abs_effect"].notna().sum()
        ),
        "independent_loci": int(lead["studyLocusId"].nunique()),
        "tier_1": int((lead["evidence_tier"] == "Tier 1").sum()),
        "tier_2": int((lead["evidence_tier"] == "Tier 2").sum()),
        "median_priority_score": float(lead["priority_score"].median()),
        "positive_control": "chr21:39094644:G>A (rs2836882, ETS2)",
    }
    (table_dir / "analysis_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))
    print("\nTop lead candidates:")
    columns = [
        "lead_rank",
        "canonical_variant",
        "rsIds",
        "l2g_target_symbol",
        "pip",
        "l2g_score",
        "atlas_gi_max_abs_q",
        "rna_specificity",
        "sc_specificity",
        "custom_macrophage_abs_percentile",
        "evidence_tier",
        "priority_score",
    ]
    print(lead[columns].head(20).to_string(index=False))


if __name__ == "__main__":
    main()
