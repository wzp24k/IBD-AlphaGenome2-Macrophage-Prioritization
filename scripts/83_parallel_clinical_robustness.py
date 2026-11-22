from __future__ import annotations

import importlib.util
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT


TABLE_DIR = ROOT / "results" / "tables"
OUTPUT = TABLE_DIR / "clinical_translation_robustness_summary.json"
N_PERMUTATIONS = 20_000
PERMUTATION_CHUNKS = max(4, min(os.cpu_count() or 4, 16))


def _load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


PQTL = _load_module("pqtl_module", Path(__file__).parent / "73_pqtl_mr_coloc_phewas.py")


def ivw(frame: pd.DataFrame, beta_col: str, se_col: str) -> tuple[float, float]:
    bx = frame[beta_col].to_numpy(float)
    sy = frame[se_col].to_numpy(float)
    by = frame["outcome_beta"].to_numpy(float)
    beta = float(np.sum(bx * by / sy**2) / np.sum((bx / sy) ** 2))
    se = float(np.sqrt(1 / np.sum((bx / sy) ** 2)))
    return beta, se


def weighted_median(
    beta_exposure: np.ndarray,
    beta_outcome: np.ndarray,
    se_outcome: np.ndarray,
) -> float:
    ratio = beta_outcome / beta_exposure
    se_ratio = se_outcome / np.abs(beta_exposure)
    weights = 1 / se_ratio**2
    order = np.argsort(ratio)
    ratio = ratio[order]
    weights = weights[order]
    cumulative = np.cumsum(weights) - 0.5 * weights
    cumulative /= np.sum(weights)
    return float(np.interp(0.5, cumulative, ratio))


def mr_egger(
    beta_exposure: np.ndarray,
    beta_outcome: np.ndarray,
    se_outcome: np.ndarray,
) -> tuple[float, float, float, float]:
    weights = 1 / se_outcome**2
    x = np.column_stack([np.ones(len(beta_exposure)), beta_exposure])
    xtw = x.T * weights
    covariance = np.linalg.pinv(xtw @ x)
    coefficients = covariance @ (xtw @ beta_outcome)
    residual = beta_outcome - x @ coefficients
    df = max(len(beta_exposure) - 2, 1)
    sigma2 = float(np.sum(weights * residual**2) / df)
    covariance *= sigma2
    slope = float(coefficients[1])
    slope_se = float(np.sqrt(max(covariance[1, 1], 0)))
    intercept = float(coefficients[0])
    intercept_p = float(
        2 * stats.norm.sf(abs(intercept / np.sqrt(max(covariance[0, 0], 1e-30))))
    )
    return slope, slope_se, intercept, intercept_p


def cross_ancestry_leave_one_out() -> pd.DataFrame:
    frame = pd.read_csv(TABLE_DIR / "prkcb_cross_ancestry_harmonized.csv")
    rows = []
    for study, group in frame.groupby("study"):
        for omitted in ["none", *group["rsid"].tolist()]:
            subset = group if omitted == "none" else group[group["rsid"].ne(omitted)]
            if subset.empty:
                continue
            beta, se = ivw(subset, "eqtl_beta", "outcome_se")
            rows.append(
                {
                    "study": study,
                    "omitted_variant": omitted,
                    "nsnp": len(subset),
                    "beta": beta,
                    "se": se,
                    "p": float(2 * stats.norm.sf(abs(beta / se))),
                }
            )
    output = pd.DataFrame(rows)
    output.to_csv(
        TABLE_DIR / "prkcb_cross_ancestry_leave_one_out.csv",
        index=False,
    )
    return output


def pqtl_mr_sensitivity() -> pd.DataFrame:
    pqtl = PQTL.load_pqtl_regions()
    outcomes = PQTL.load_gwas_catalog_regions()
    rows = []
    harmonized_tables = []
    for dataset_name, exposure in pqtl.items():
        instruments = PQTL.select_instruments(exposure)
        for outcome_name, outcome in outcomes.items():
            harmonized = PQTL.harmonize(
                instruments,
                outcome,
                merge_key="rsid",
            )
            if harmonized.empty:
                continue
            harmonized.insert(0, "outcome", outcome_name)
            harmonized.insert(0, "pqtl_dataset", dataset_name)
            harmonized_tables.append(harmonized)
            ivw_result = PQTL.mr_estimate(harmonized)
            rows.append(
                {
                    "pqtl_dataset": dataset_name,
                    "outcome": outcome_name,
                    "method": "IVW",
                    "nsnp": ivw_result["nsnp"],
                    "beta": ivw_result["beta"],
                    "se": ivw_result["se"],
                    "p": ivw_result["p"],
                }
            )
            if len(harmonized) >= 3:
                median = weighted_median(
                    harmonized["exposure_beta"].to_numpy(float),
                    harmonized["outcome_beta"].to_numpy(float),
                    harmonized["outcome_se"].to_numpy(float),
                )
                rows.append(
                    {
                        "pqtl_dataset": dataset_name,
                        "outcome": outcome_name,
                        "method": "weighted_median",
                        "nsnp": len(harmonized),
                        "beta": median,
                        "se": np.nan,
                        "p": np.nan,
                    }
                )
                egger, egger_se, intercept, intercept_p = mr_egger(
                    harmonized["exposure_beta"].to_numpy(float),
                    harmonized["outcome_beta"].to_numpy(float),
                    harmonized["outcome_se"].to_numpy(float),
                )
                rows.append(
                    {
                        "pqtl_dataset": dataset_name,
                        "outcome": outcome_name,
                        "method": "MR_Egger",
                        "nsnp": len(harmonized),
                        "beta": egger,
                        "se": egger_se,
                        "p": float(2 * stats.norm.sf(abs(egger / egger_se))),
                        "egger_intercept": intercept,
                        "egger_intercept_p": intercept_p,
                    }
                )
                for omitted in harmonized["rsid"].dropna().unique():
                    subset = harmonized[harmonized["rsid"].ne(omitted)]
                    if len(subset) < 2:
                        continue
                    beta, se = ivw(subset, "exposure_beta", "outcome_se")
                    rows.append(
                        {
                            "pqtl_dataset": dataset_name,
                            "outcome": outcome_name,
                            "method": f"IVW_without_{omitted}",
                            "nsnp": len(subset),
                            "beta": beta,
                            "se": se,
                            "p": float(2 * stats.norm.sf(abs(beta / se))),
                        }
                    )
    output = pd.DataFrame(rows)
    output.to_csv(
        TABLE_DIR / "pqtl_mr_sensitivity.csv",
        index=False,
    )
    if harmonized_tables:
        pd.concat(harmonized_tables, ignore_index=True).to_csv(
            TABLE_DIR / "pqtl_ibd_harmonized.csv",
            index=False,
        )
    return output


def pqtl_coloc_sensitivity() -> pd.DataFrame:
    pqtl = PQTL.load_pqtl_regions()
    outcomes = PQTL.load_gwas_catalog_regions()
    settings = [
        (sd_prior, p12)
        for sd_prior in [0.05, 0.1, 0.2]
        for p12 in [1e-6, 5e-6, 1e-5]
    ]
    rows = []
    for dataset_name, exposure in pqtl.items():
        for outcome_name, outcome in outcomes.items():
            harmonized = PQTL.harmonize(
                exposure,
                outcome,
                merge_key="rsid",
            )
            if len(harmonized) < 50:
                continue
            frame = pd.DataFrame(
                {
                    "eql_beta": harmonized["exposure_beta"],
                    "eql_se": harmonized["exposure_se"],
                    "gwas_beta": harmonized["outcome_beta"],
                    "gwas_se": harmonized["outcome_se"],
                }
            )
            for sd_prior, p12 in settings:
                posterior, _ = PQTL.COLOC.coloc_abf(
                    frame,
                    p1=1e-4,
                    p2=1e-4,
                    p12=p12,
                    sd_prior=sd_prior,
                )
                rows.append(
                    {
                        "pqtl_dataset": dataset_name,
                        "outcome": outcome_name,
                        "sd_prior": sd_prior,
                        "p12": p12,
                        "variants": len(frame),
                        **posterior,
                    }
                )
    output = pd.DataFrame(rows)
    output.to_csv(
        TABLE_DIR / "pqtl_coloc_prior_sensitivity.csv",
        index=False,
    )
    return output


def _smd(values: np.ndarray, labels: np.ndarray) -> float:
    first = values[labels]
    second = values[~labels]
    if len(first) < 2 or len(second) < 2:
        return np.nan
    pooled = (
        (len(first) - 1) * first.var(ddof=1)
        + (len(second) - 1) * second.var(ddof=1)
    ) / (len(first) + len(second) - 2)
    if pooled <= 0:
        return 0.0
    return float((first.mean() - second.mean()) / np.sqrt(pooled))


def _meta_permutation(
    values_a: np.ndarray,
    labels_a: np.ndarray,
    values_b: np.ndarray,
    labels_b: np.ndarray,
    permutations: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    meta = np.empty(permutations, dtype=float)
    for index in range(permutations):
        shuffled_a = rng.permutation(labels_a)
        shuffled_b = rng.permutation(labels_b)
        smd_a = _smd(values_a, shuffled_a)
        smd_b = _smd(values_b, shuffled_b)
        # Equal cohort weights avoid treating the six-sample cohort as highly precise.
        meta[index] = 0.5 * (smd_a + smd_b)
    return meta, rng


def _permutation_chunk(payload):
    values_a, labels_a, values_b, labels_b, permutations, seed = payload
    return _meta_permutation(
        values_a,
        labels_a,
        values_b,
        labels_b,
        permutations,
        seed,
    )[0]


def treatment_response_permutations() -> pd.DataFrame:
    scores = pd.read_csv(TABLE_DIR / "treatment_response_sample_scores.csv")
    meta = pd.read_csv(TABLE_DIR / "treatment_response_meta_analysis.csv")
    features = [
        "PRKCB_expression",
        "PDGFB_expression",
        "PRKCB_network_program",
        "macrophage_marker_program",
    ]
    rows = []
    for feature in features:
        data_a = scores[scores["dataset"].eq("GSE14580")].dropna(
            subset=[feature, "response"]
        )
        data_b = scores[scores["dataset"].eq("GSE111761")].dropna(
            subset=[feature, "response"]
        )
        values_a = data_a[feature].to_numpy(float)
        labels_a = data_a["response"].astype(bool).to_numpy()
        values_b = data_b[feature].to_numpy(float)
        labels_b = data_b["response"].astype(bool).to_numpy()
        observed = 0.5 * (
            _smd(values_a, labels_a) + _smd(values_b, labels_b)
        )
        chunk_size = int(np.ceil(N_PERMUTATIONS / PERMUTATION_CHUNKS))
        payloads = [
            (
                values_a,
                labels_a,
                values_b,
                labels_b,
                min(chunk_size, N_PERMUTATIONS - index * chunk_size),
                104729 + index,
            )
            for index in range(PERMUTATION_CHUNKS)
            if index * chunk_size < N_PERMUTATIONS
        ]
        with ProcessPoolExecutor(max_workers=len(payloads)) as executor:
            chunks = list(executor.map(_permutation_chunk, payloads))
        null = np.concatenate(chunks)
        p = float(
            (1 + np.sum(np.abs(null) >= abs(observed)))
            / (len(null) + 1)
        )
        rows.append(
            {
                "feature": feature,
                "observed_meta_smd": observed,
                "permutations": len(null),
                "permutation_p_two_sided": p,
            }
        )
    output = pd.DataFrame(rows).merge(meta, on="feature", how="left")
    output["meta_ci_excludes_zero"] = (
        output["ci_low"].gt(0) | output["ci_high"].lt(0)
    )
    output["robust_permutation_significant"] = (
        output["permutation_p_two_sided"].lt(0.05)
        & output["meta_ci_excludes_zero"]
    )
    output.to_csv(
        TABLE_DIR / "treatment_response_permutation.csv",
        index=False,
    )
    return output


def readiness_weight_bootstrap() -> tuple[pd.DataFrame, dict[str, float]]:
    frame = pd.read_csv(TABLE_DIR / "clinical_development_readiness_matrix.csv")
    columns = [
        "genetic_transport_score",
        "disease_specific_score",
        "exact_eqtl_score",
        "colocalization_score",
        "public_perturbation_score",
        "pqtl_mr_score",
        "pqtl_colocalization_score",
        "treatment_response_score",
        "clinical_stage_score",
        "direction_compatible_score",
    ]
    components = frame[columns].to_numpy(float)
    n_draws = 100_000
    rng = np.random.default_rng(20261122)
    weights = rng.dirichlet(np.ones(len(columns)), size=n_draws)
    score_draws = weights @ components.T
    ranks = (-score_draws).argsort(axis=1).argsort(axis=1) + 1
    rows = []
    for index, row in frame.iterrows():
        rows.append(
            {
                "target": row["target"],
                "primary_rank": int(
                    frame["readiness_score_0_10"].rank(
                        ascending=False,
                        method="min",
                    ).iloc[index]
                ),
                "probability_rank_1": float(np.mean(ranks[:, index] == 1)),
                "probability_top_2": float(np.mean(ranks[:, index] <= 2)),
                "probability_top_3": float(np.mean(ranks[:, index] <= 3)),
                "median_rank": float(np.median(ranks[:, index])),
            }
        )
    output = pd.DataFrame(rows)
    output.to_csv(
        TABLE_DIR / "clinical_readiness_weight_bootstrap.csv",
        index=False,
    )
    prkcb = output[output["target"].eq("PRKCB")].iloc[0]
    summary = {
        "draws": n_draws,
        "prkcb_probability_rank_1": float(prkcb["probability_rank_1"]),
        "prkcb_probability_top_2": float(prkcb["probability_top_2"]),
        "prkcb_median_rank": float(prkcb["median_rank"]),
    }
    return output, summary


def main() -> None:
    cross_loo = cross_ancestry_leave_one_out()
    pqtl_mr = pqtl_mr_sensitivity()
    pqtl_coloc = pqtl_coloc_sensitivity()
    treatment_perm = treatment_response_permutations()
    readiness_ranks, readiness_summary = readiness_weight_bootstrap()

    summary = {
        "cpu_count": os.cpu_count(),
        "permutation_chunks": PERMUTATION_CHUNKS,
        "treatment_permutations_per_feature": N_PERMUTATIONS,
        "readiness_weight_draws": readiness_summary["draws"],
        "prkcb_cross_ancestry_all_leave_one_out_positive": bool(
            (cross_loo["beta"] > 0).all()
        ),
        "prkcb_cross_ancestry_min_leave_one_out_beta": float(
            cross_loo["beta"].min()
        ),
        "prkcb_cross_ancestry_max_leave_one_out_beta": float(
            cross_loo["beta"].max()
        ),
        "treatment_permutation_significant_features": treatment_perm.loc[
            treatment_perm["robust_permutation_significant"],
            "feature",
        ].tolist(),
        "treatment_permutation_min_p": float(
            treatment_perm["permutation_p_two_sided"].min()
        ),
        "pqtl_sensitivity_tests": int(len(pqtl_mr)),
        "pqtl_positive_sensitivity_fraction": float(
            (pqtl_mr["beta"] > 0).mean()
        ),
        "pqtl_coloc_sensitivity_tests": int(len(pqtl_coloc)),
        "pqtl_coloc_sensitivity_max_h4": float(pqtl_coloc["H4"].max()),
        **readiness_summary,
    }
    OUTPUT.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
