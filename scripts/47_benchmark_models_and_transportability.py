from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    brier_score_loss,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT


BOOTSTRAPS = 500
SEED = 20261114
SCORES = {
    "PIP-only": "pip",
    "L2G-only": "l2g_score",
    "Atlas-only": "atlas_global_max_abs_q",
    "Live-only": "live_specificity",
    "Combined": "framework_pair_score",
}
FEATURES = [
    "pip",
    "l2g_score",
    "atlas_global_max_abs_q",
    "live_specificity",
    "custom_macrophage_percentile",
]


def build_pair_benchmark() -> pd.DataFrame:
    table_dir = ROOT / "results" / "tables"
    candidates = pd.read_csv(
        table_dir / "ranked_variant_gene_candidates_final.csv"
    )
    candidates["query_key"] = (
        candidates["variant_id"].astype(str)
        + "|"
        + candidates["l2g_target_id"].astype(str)
    )
    eqtl = pd.read_csv(table_dir / "external_eqtl_variant_benchmark.csv")
    eqtl = eqtl[eqtl["query_group"].eq("primary")].copy()
    frame = eqtl.merge(
        candidates[
            [
                "query_key",
                "studyLocusId",
                "variant_id",
                "canonical_variant",
                "pip",
                "l2g_score",
                "atlas_global_max_abs_q",
                "atlas_gi_max_abs_q",
                "live_specificity",
                "custom_macrophage_percentile",
                "pip_percentile",
                "atlas_gi_percentile",
                "live_specificity_percentile",
                "priority_score",
                "external_cluster_count",
                "l2g_target_id",
                "l2g_target_symbol",
                "final_evidence_tier",
                "experimental_evidence_level",
                "experimental_support",
            ]
        ],
        on="query_key",
        how="left",
    )
    for column in SCORES.values():
        if column in frame.columns:
            frame[column] = pd.to_numeric(frame[column], errors="coerce")
    frame["l2g_pair_percentile"] = frame["l2g_score"].rank(pct=True)
    frame["framework_pair_score"] = frame[
        [
            "pip_percentile",
            "l2g_pair_percentile",
            "atlas_gi_percentile",
            "live_specificity_percentile",
            "custom_macrophage_percentile",
        ]
    ].fillna(0).mean(axis=1)
    frame["external_any"] = (
        pd.to_numeric(frame["external_cluster_count"], errors="coerce")
        .fillna(0)
        .gt(0)
        .astype(int)
    )
    frame["external_both"] = (
        pd.to_numeric(frame["external_cluster_count"], errors="coerce")
        .fillna(0)
        .ge(2)
        .astype(int)
    )
    frame["eqtl_any"] = frame["any_nominal"].astype(int)
    frame["eqtl_replicated"] = frame["replicated_nominal"].astype(int)

    crispr = pd.read_csv(
        table_dir / "experimental_primary_macrophage_crispr_genes.csv"
    )
    crispr = crispr[crispr["screen"].eq("TNF")][["Gene", "screen_fdr"]]
    crispr = crispr.rename(
        columns={
            "Gene": "l2g_target_symbol",
            "screen_fdr": "crispr_tnf_fdr",
        }
    )
    frame = frame.merge(crispr, on="l2g_target_symbol", how="left")
    frame["crispr_tnf_supported"] = (
        frame["crispr_tnf_fdr"].notna()
        & frame["crispr_tnf_fdr"].lt(0.05)
    ).astype(int)
    frame["crispr_evaluable"] = frame["crispr_tnf_fdr"].notna()
    return frame


def safe_metrics(labels: np.ndarray, scores: np.ndarray) -> tuple[float, float]:
    if len(np.unique(labels)) < 2:
        return np.nan, np.nan
    return float(roc_auc_score(labels, scores)), float(
        average_precision_score(labels, scores)
    )


def bootstrap_metrics(
    frame: pd.DataFrame,
    label_column: str,
    score_column: str,
    group_column: str = "studyLocusId",
) -> dict:
    group_values = list(frame.groupby(group_column).indices.values())
    labels_all = frame[label_column].to_numpy()
    scores_all = frame[score_column].fillna(0).to_numpy()
    rng = np.random.default_rng(SEED)
    aucs = []
    aps = []
    for _ in range(BOOTSTRAPS):
        sampled = rng.integers(0, len(group_values), len(group_values))
        indices = np.concatenate([group_values[index] for index in sampled])
        labels = labels_all[indices]
        scores = scores_all[indices]
        auc, ap = safe_metrics(labels, scores)
        if np.isfinite(auc):
            aucs.append(auc)
            aps.append(ap)
    labels = labels_all
    scores = scores_all
    auc, ap = safe_metrics(labels, scores)
    return {
        "n": len(frame),
        "positives": int(labels.sum()),
        "negatives": int(len(labels) - labels.sum()),
        "auroc": auc,
        "auroc_ci_low": float(np.quantile(aucs, 0.025)) if aucs else np.nan,
        "auroc_ci_high": float(np.quantile(aucs, 0.975)) if aucs else np.nan,
        "auprc": ap,
        "auprc_ci_low": float(np.quantile(aps, 0.025)) if aps else np.nan,
        "auprc_ci_high": float(np.quantile(aps, 0.975)) if aps else np.nan,
        "valid_bootstrap_replicates": len(aucs),
    }


def bootstrap_delta(
    frame: pd.DataFrame,
    label_column: str,
    combined_column: str,
    baseline_column: str,
    group_column: str = "studyLocusId",
) -> dict:
    group_values = list(frame.groupby(group_column).indices.values())
    labels_all = frame[label_column].to_numpy()
    combined_all = frame[combined_column].fillna(0).to_numpy()
    baseline_all = frame[baseline_column].fillna(0).to_numpy()
    rng = np.random.default_rng(SEED + 1)
    auc_deltas = []
    ap_deltas = []
    for _ in range(BOOTSTRAPS):
        sampled = rng.integers(0, len(group_values), len(group_values))
        indices = np.concatenate([group_values[index] for index in sampled])
        labels = labels_all[indices]
        combined_auc, combined_ap = safe_metrics(
            labels,
            combined_all[indices],
        )
        baseline_auc, baseline_ap = safe_metrics(
            labels,
            baseline_all[indices],
        )
        if np.isfinite(combined_auc) and np.isfinite(baseline_auc):
            auc_deltas.append(combined_auc - baseline_auc)
            ap_deltas.append(combined_ap - baseline_ap)
    labels = labels_all
    combined_auc, combined_ap = safe_metrics(
        labels,
        combined_all,
    )
    baseline_auc, baseline_ap = safe_metrics(
        labels,
        baseline_all,
    )
    return {
        "delta_auroc": combined_auc - baseline_auc,
        "delta_auroc_ci_low": float(np.quantile(auc_deltas, 0.025))
        if auc_deltas
        else np.nan,
        "delta_auroc_ci_high": float(np.quantile(auc_deltas, 0.975))
        if auc_deltas
        else np.nan,
        "delta_auprc": combined_ap - baseline_ap,
        "delta_auprc_ci_low": float(np.quantile(ap_deltas, 0.025))
        if ap_deltas
        else np.nan,
        "delta_auprc_ci_high": float(np.quantile(ap_deltas, 0.975))
        if ap_deltas
        else np.nan,
        "valid_bootstrap_replicates": len(auc_deltas),
    }


def expected_calibration_error(
    labels: np.ndarray,
    probabilities: np.ndarray,
    bins: int = 10,
) -> float:
    edges = np.linspace(0, 1, bins + 1)
    total = len(labels)
    error = 0.0
    for index in range(bins):
        lower, upper = edges[index], edges[index + 1]
        if index == bins - 1:
            mask = (probabilities >= lower) & (probabilities <= upper)
        else:
            mask = (probabilities >= lower) & (probabilities < upper)
        if not mask.any():
            continue
        observed = labels[mask].mean()
        predicted = probabilities[mask].mean()
        error += mask.sum() / total * abs(observed - predicted)
    return float(error)


def out_of_fold_calibration(
    frame: pd.DataFrame,
    label_column: str,
    score_column: str,
    group_column: str,
) -> dict | None:
    data = frame.dropna(subset=[label_column, score_column, group_column]).copy()
    labels = data[label_column].to_numpy(dtype=int)
    if labels.sum() < 10 or (1 - labels).sum() < 10:
        return None
    splitter = StratifiedGroupKFold(
        n_splits=5,
        shuffle=True,
        random_state=SEED,
    )
    probabilities = np.full(len(data), np.nan)
    for train, test in splitter.split(
        data[[score_column]],
        labels,
        groups=data[group_column],
    ):
        scaler = StandardScaler()
        train_x = scaler.fit_transform(data.iloc[train][[score_column]])
        test_x = scaler.transform(data.iloc[test][[score_column]])
        model = LogisticRegression(
            C=1.0,
            solver="lbfgs",
            max_iter=2000,
        )
        model.fit(train_x, labels[train])
        probabilities[test] = model.predict_proba(test_x)[:, 1]
    valid = np.isfinite(probabilities)
    clipped = np.clip(probabilities[valid], 1e-6, 1 - 1e-6)
    logits = np.log(clipped / (1 - clipped)).reshape(-1, 1)
    calibration_model = LogisticRegression(
        C=1e6,
        solver="lbfgs",
        max_iter=2000,
    )
    calibration_model.fit(logits, labels[valid])
    return {
        "oof_auroc": float(roc_auc_score(labels[valid], probabilities[valid])),
        "oof_auprc": float(
            average_precision_score(labels[valid], probabilities[valid])
        ),
        "brier_score": float(
            brier_score_loss(labels[valid], probabilities[valid])
        ),
        "expected_calibration_error": expected_calibration_error(
            labels[valid],
            probabilities[valid],
        ),
        "calibration_intercept": float(calibration_model.intercept_[0]),
        "calibration_slope": float(calibration_model.coef_[0, 0]),
        "calibration_n": int(valid.sum()),
    }


def target_level_features(frame: pd.DataFrame) -> pd.DataFrame:
    features = (
        frame.groupby(["l2g_target_id", "l2g_target_symbol"], as_index=False)
        .agg(
            pip=("pip", "max"),
            l2g_score=("l2g_score", "max"),
            atlas_global_max_abs_q=("atlas_global_max_abs_q", "max"),
            live_specificity=("live_specificity", "max"),
            custom_macrophage_percentile=(
                "custom_macrophage_percentile",
                "max",
            ),
            framework_pair_score=("framework_pair_score", "max"),
        )
    )
    external = pd.read_csv(
        ROOT / "results" / "tables" / "external_ibd_l2g_predictions.csv"
    ).dropna(subset=["target_id"])
    labels = (
        external.assign(
            ukb=external["cluster"].eq("UKB_WGS_2025").astype(int),
            mvp=external["cluster"].eq("MVP_2024").astype(int),
        )
        .groupby("target_id", as_index=False)
        .agg(ukb=("ukb", "max"), mvp=("mvp", "max"))
    )
    features["l2g_target_id"] = features["l2g_target_id"].astype(str)
    labels["target_id"] = labels["target_id"].astype(str)
    return features.merge(
        labels,
        left_on="l2g_target_id",
        right_on="target_id",
        how="left",
    ).fillna({"ukb": 0, "mvp": 0})


def logistic_transport(
    frame: pd.DataFrame,
    train_label: str,
    test_label: str,
) -> dict:
    train = frame.dropna(subset=FEATURES).copy()
    train_label_values = train[train_label].astype(int).to_numpy()
    test_label_values = train[test_label].astype(int).to_numpy()
    if train_label_values.sum() < 5 or (1 - train_label_values).sum() < 5:
        return {"error": "insufficient_training_labels"}
    if test_label_values.sum() < 5 or (1 - test_label_values).sum() < 5:
        return {"error": "insufficient_test_labels"}
    scaler = StandardScaler()
    x = scaler.fit_transform(train[FEATURES])
    model = LogisticRegression(
        C=1.0,
        solver="lbfgs",
        max_iter=2000,
        class_weight="balanced",
    )
    model.fit(x, train_label_values)
    probabilities = model.predict_proba(x)[:, 1]
    return {
        "train_label": train_label,
        "test_label": test_label,
        "train_targets": len(train),
        "train_positives": int(train_label_values.sum()),
        "test_targets": len(train),
        "test_positives": int(test_label_values.sum()),
        "auroc": float(roc_auc_score(test_label_values, probabilities)),
        "auprc": float(
            average_precision_score(test_label_values, probabilities)
        ),
        "brier_score": float(
            brier_score_loss(test_label_values, probabilities)
        ),
        "note": (
            "In-sample feature fit with held-out cohort label. UKB and MVP are "
            "cohort-transportability proxies, not a formal ancestry-specific test."
        ),
    }


def heldout_group_validation(
    frame: pd.DataFrame,
    label_column: str,
    group_column: str,
) -> dict:
    data = frame.dropna(
        subset=FEATURES + [label_column, group_column]
    ).copy()
    labels = data[label_column].to_numpy(dtype=int)
    if labels.sum() < 10 or (1 - labels).sum() < 10:
        return {"error": "insufficient_labels"}
    splitter = StratifiedGroupKFold(
        n_splits=5,
        shuffle=True,
        random_state=SEED,
    )
    probabilities = np.full(len(data), np.nan)
    for train, test in splitter.split(
        data[FEATURES],
        labels,
        groups=data[group_column],
    ):
        scaler = StandardScaler()
        train_x = scaler.fit_transform(data.iloc[train][FEATURES])
        test_x = scaler.transform(data.iloc[test][FEATURES])
        model = LogisticRegression(
            C=1.0,
            solver="lbfgs",
            max_iter=2000,
            class_weight="balanced",
        )
        model.fit(train_x, labels[train])
        probabilities[test] = model.predict_proba(test_x)[:, 1]
    valid = np.isfinite(probabilities)
    return {
        "group_column": group_column,
        "n": int(valid.sum()),
        "positives": int(labels[valid].sum()),
        "auroc": float(roc_auc_score(labels[valid], probabilities[valid])),
        "auprc": float(
            average_precision_score(labels[valid], probabilities[valid])
        ),
        "brier_score": float(
            brier_score_loss(labels[valid], probabilities[valid])
        ),
        "note": (
            "Nested grouped logistic benchmark; the original ranking was not "
            "refitted or changed."
        ),
    }


def main() -> None:
    table_dir = ROOT / "results" / "tables"
    frame = build_pair_benchmark()
    frame.to_csv(table_dir / "benchmark_pair_level_dataset.csv", index=False)

    task_specs = [
        ("External genetic support, any cohort", "external_any", None),
        ("External genetic support, both cohorts", "external_both", None),
        ("Exact-variant eQTL, any dataset", "eqtl_any", None),
        ("Exact-variant eQTL, replicated", "eqtl_replicated", None),
        ("Macrophage TNF CRISPR", "crispr_tnf_supported", "crispr_evaluable"),
    ]
    metric_rows = []
    delta_rows = []
    calibration_rows = []
    for task_name, label_column, evaluable_column in task_specs:
        subset = frame.copy()
        if evaluable_column:
            subset = subset[subset[evaluable_column]].copy()
        for model_name, score_column in SCORES.items():
            print(f"Benchmarking: {task_name} / {model_name}", flush=True)
            metrics = bootstrap_metrics(
                subset,
                label_column,
                score_column,
            )
            metric_rows.append(
                {
                    "task": task_name,
                    "model": model_name,
                    "label": label_column,
                    **metrics,
                }
            )
            calibration = out_of_fold_calibration(
                subset,
                label_column,
                score_column,
                group_column="studyLocusId",
            )
            if calibration:
                calibration_rows.append(
                    {
                        "task": task_name,
                        "model": model_name,
                        **calibration,
                    }
                )
            if model_name != "Combined":
                delta_rows.append(
                    {
                        "task": task_name,
                        "baseline": model_name,
                        "comparison": "Combined",
                        **bootstrap_delta(
                            subset,
                            label_column,
                            "framework_pair_score",
                            score_column,
                        ),
                    }
                )

    metrics = pd.DataFrame(metric_rows)
    deltas = pd.DataFrame(delta_rows)
    calibration = pd.DataFrame(calibration_rows)
    metrics.to_csv(
        table_dir / "benchmark_model_metrics.csv",
        index=False,
    )
    deltas.to_csv(
        table_dir / "benchmark_model_deltas.csv",
        index=False,
    )
    calibration.to_csv(
        table_dir / "benchmark_calibration_metrics.csv",
        index=False,
    )

    target_frame = target_level_features(frame)
    transport = pd.DataFrame(
        [
            logistic_transport(target_frame, "ukb", "mvp"),
            logistic_transport(target_frame, "mvp", "ukb"),
        ]
    )
    transport.to_csv(
        table_dir / "benchmark_cohort_transportability.csv",
        index=False,
    )
    heldout = pd.DataFrame(
        [
            {
                "label": "external_both",
                **heldout_group_validation(
                    frame,
                    "external_both",
                    "studyLocusId",
                ),
            },
            {
                "label": "external_both",
                **heldout_group_validation(
                    frame,
                    "external_both",
                    "l2g_target_id",
                ),
            },
        ]
    )
    heldout.to_csv(
        table_dir / "benchmark_heldout_validation.csv",
        index=False,
    )

    eqtl_summary = json.loads(
        (
            table_dir / "external_eqtl_benchmark_summary.json"
        ).read_text(encoding="utf-8")
    )
    perturb = pd.read_csv(table_dir / "perturbseq_gse327644_target_knockdown.csv")
    perturb_targets = perturb.groupby("target", as_index=False).agg(
        minimum_self_log2fc=("mean_log2fc_vs_ntc", "min"),
        minimum_pvalue=("p_value", "min"),
    )
    coverage = {
        "pair_level_rows": len(frame),
        "external_genetic_any_positives": int(frame["external_any"].sum()),
        "external_genetic_both_positives": int(frame["external_both"].sum()),
        "eqtl_rows": int(len(frame)),
        "eqtl_any_positives": int(frame["eqtl_any"].sum()),
        "eqtl_replicated_positives": int(frame["eqtl_replicated"].sum()),
        "crispr_evaluable_pairs": int(frame["crispr_evaluable"].sum()),
        "crispr_tnf_positives": int(frame["crispr_tnf_supported"].sum()),
        "top10_exact_variant_eqtl": eqtl_summary,
        "perturbseq_effective_targets": int(
            (
                perturb_targets["minimum_self_log2fc"].lt(-0.3)
                & perturb_targets["minimum_pvalue"].lt(0.05)
            ).sum()
        ),
        "perturbseq_evaluated_targets": int(len(perturb_targets)),
        "mpra_exact_variant_positives": 1,
        "mpra_exact_variant_note": (
            "One exact-variant macrophage MPRA positive control (ETS2); "
            "insufficient support for an AUROC estimate."
        ),
        "note": (
            "AUROC/AUPRC were calculated only for tasks with at least ten "
            "positive and ten negative labels. Sparse MPRA and Perturb-seq "
            "evidence is reported as coverage, not used to fabricate metrics."
        ),
    }
    (table_dir / "benchmark_coverage_summary.json").write_text(
        json.dumps(coverage, indent=2),
        encoding="utf-8",
    )
    print(metrics.to_string(index=False))
    print(deltas.to_string(index=False))
    print(calibration.to_string(index=False))
    print(transport.to_string(index=False))
    print(heldout.to_string(index=False))
    print(json.dumps(coverage, indent=2))


if __name__ == "__main__":
    main()
