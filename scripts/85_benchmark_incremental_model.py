from __future__ import annotations

import json
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[1]
TABLE_DIR = ROOT / "results" / "tables"
SEED = 20261122
REPEATS = 30
N_SPLITS = 5
WORKERS = max(1, min(os.cpu_count() or 1, 44))


FEATURE_SETS = {
    "PIP-only": ["pip"],
    "PIP+L2G": ["pip", "l2g_score"],
    "PIP+Atlas": ["pip", "atlas_global_max_abs_q"],
    "PIP+Live": ["pip", "live_specificity"],
    "PIP+context": [
        "pip",
        "l2g_score",
        "atlas_global_max_abs_q",
        "live_specificity",
        "custom_macrophage_percentile",
    ],
}

TASKS = {
    "External genetic support, any cohort": "external_any",
    "External genetic support, both cohorts": "external_both",
    "Exact-variant eQTL, replicated": "eqtl_replicated",
}


def make_model(c_value: float) -> Pipeline:
    return Pipeline(
        [
            ("impute", SimpleImputer(strategy="median")),
            ("scale", StandardScaler()),
            (
                "logistic",
                LogisticRegression(
                    C=c_value,
                    class_weight="balanced",
                    solver="lbfgs",
                    max_iter=3000,
                ),
            ),
        ]
    )


def out_of_fold_predictions(
    frame: pd.DataFrame,
    label: str,
    features: list[str],
    c_value: float,
    random_state: int,
) -> np.ndarray:
    labels = frame[label].to_numpy(dtype=int)
    probabilities = np.full(len(frame), np.nan)
    splitter = StratifiedGroupKFold(
        n_splits=N_SPLITS,
        shuffle=True,
        random_state=random_state,
    )
    for train, test in splitter.split(
        frame[features],
        labels,
        groups=frame["studyLocusId"],
    ):
        model = make_model(c_value)
        model.fit(frame.iloc[train][features], labels[train])
        probabilities[test] = model.predict_proba(
            frame.iloc[test][features]
        )[:, 1]
    if not np.isfinite(probabilities).all():
        raise RuntimeError("Out-of-fold predictions contain missing values.")
    return probabilities


def evaluate_configuration(payload) -> dict:
    (
        task_name,
        label,
        feature_name,
        features,
        c_value,
        frame,
        repeats,
    ) = payload
    labels = frame[label].to_numpy(dtype=int)
    repeated_metrics = []
    probability_sum = np.zeros(len(frame), dtype=float)
    for repeat in range(repeats):
        probabilities = out_of_fold_predictions(
            frame,
            label,
            features,
            c_value,
            SEED + repeat,
        )
        probability_sum += probabilities
        repeated_metrics.append(metrics(labels, probabilities))
    repeated_metrics_array = np.asarray(repeated_metrics)
    consensus_probabilities = probability_sum / repeats
    auc, ap = metrics(labels, consensus_probabilities)
    return {
        "task": task_name,
        "feature_set": feature_name,
        "c_value": c_value,
        "n": len(frame),
        "positives": int(labels.sum()),
        "repetitions": repeats,
        "mean_auroc": float(repeated_metrics_array[:, 0].mean()),
        "auroc_sd_across_repeats": float(
            repeated_metrics_array[:, 0].std(ddof=1)
        ),
        "mean_auprc": float(repeated_metrics_array[:, 1].mean()),
        "auprc_sd_across_repeats": float(
            repeated_metrics_array[:, 1].std(ddof=1)
        ),
        "consensus_auroc": auc,
        "consensus_auprc": ap,
    }


def metrics(labels: np.ndarray, scores: np.ndarray) -> tuple[float, float]:
    return (
        float(roc_auc_score(labels, scores)),
        float(average_precision_score(labels, scores)),
    )


def bootstrap_delta(
    frame: pd.DataFrame,
    labels: np.ndarray,
    baseline: np.ndarray,
    combined: np.ndarray,
    repetitions: int = 2000,
) -> dict[str, float]:
    groups = list(frame.groupby("studyLocusId").indices.values())
    rng = np.random.default_rng(SEED + 9)
    auc_deltas = []
    ap_deltas = []
    for _ in range(repetitions):
        sampled = rng.integers(0, len(groups), len(groups))
        indices = np.concatenate([groups[index] for index in sampled])
        if len(np.unique(labels[indices])) < 2:
            continue
        baseline_auc, baseline_ap = metrics(
            labels[indices],
            baseline[indices],
        )
        combined_auc, combined_ap = metrics(
            labels[indices],
            combined[indices],
        )
        auc_deltas.append(combined_auc - baseline_auc)
        ap_deltas.append(combined_ap - baseline_ap)
    return {
        "delta_auroc_ci_low": float(np.quantile(auc_deltas, 0.025)),
        "delta_auroc_ci_high": float(np.quantile(auc_deltas, 0.975)),
        "delta_auprc_ci_low": float(np.quantile(ap_deltas, 0.025)),
        "delta_auprc_ci_high": float(np.quantile(ap_deltas, 0.975)),
        "bootstrap_replicates": len(auc_deltas),
    }


def main() -> None:
    frame_all = pd.read_csv(TABLE_DIR / "benchmark_pair_level_dataset.csv")
    payloads = []
    for task_name, label in TASKS.items():
        frame = frame_all.dropna(
            subset=[label, "studyLocusId", "pip"]
        ).reset_index(drop=True)
        for feature_name, features in FEATURE_SETS.items():
            for c_value in [0.03, 0.1, 0.3, 1.0]:
                payloads.append(
                    (
                        task_name,
                        label,
                        feature_name,
                        features,
                        c_value,
                        frame,
                        REPEATS,
                    )
                )

    with ProcessPoolExecutor(max_workers=WORKERS) as executor:
        rows = list(
            executor.map(
                evaluate_configuration,
                payloads,
                chunksize=1,
            )
        )

    metrics_table = pd.DataFrame(rows)
    metrics_table.to_csv(
        TABLE_DIR / "benchmark_incremental_model_metrics.csv",
        index=False,
    )

    delta_rows = []
    for task_name, label in TASKS.items():
        frame = frame_all.dropna(
            subset=[label, "studyLocusId", "pip"]
        ).reset_index(drop=True)
        labels = frame[label].to_numpy(dtype=int)
        primary_baseline = out_of_fold_predictions(
            frame,
            label,
            FEATURE_SETS["PIP-only"],
            0.1,
            SEED,
        )
        primary_combined = out_of_fold_predictions(
            frame,
            label,
            FEATURE_SETS["PIP+context"],
            0.1,
            SEED,
        )
        primary_baseline_auc, primary_baseline_ap = metrics(
            labels,
            primary_baseline,
        )
        primary_combined_auc, primary_combined_ap = metrics(
            labels,
            primary_combined,
        )
        delta_rows.append(
            {
                "task": task_name,
                "analysis": "prespecified",
                "selected_feature_set": "PIP+context",
                "selected_c_value": 0.1,
                "baseline_auroc": primary_baseline_auc,
                "combined_auroc": primary_combined_auc,
                "delta_auroc": primary_combined_auc - primary_baseline_auc,
                "baseline_auprc": primary_baseline_ap,
                "combined_auprc": primary_combined_ap,
                "delta_auprc": primary_combined_ap - primary_baseline_ap,
                **bootstrap_delta(
                    frame,
                    labels,
                    primary_baseline,
                    primary_combined,
                ),
                "selection_note": (
                    "PIP+context with C=0.1 was specified before repeated "
                    "grouped evaluation."
                ),
            }
        )
        task_metrics = metrics_table[metrics_table["task"].eq(task_name)]
        best = task_metrics.sort_values(
            ["mean_auroc", "mean_auprc"],
            ascending=False,
        ).iloc[0]
        best_features = FEATURE_SETS[best["feature_set"]]
        c_value = float(best["c_value"])
        baseline = out_of_fold_predictions(
            frame,
            label,
            ["pip"],
            c_value,
            SEED,
        )
        combined = out_of_fold_predictions(
            frame,
            label,
            best_features,
            c_value,
            SEED,
        )
        baseline_auc, baseline_ap = metrics(labels, baseline)
        combined_auc, combined_ap = metrics(labels, combined)
        delta_rows.append(
            {
                "task": task_name,
                "analysis": "posthoc_best",
                "selected_feature_set": best["feature_set"],
                "selected_c_value": c_value,
                "baseline_auroc": baseline_auc,
                "combined_auroc": combined_auc,
                "delta_auroc": combined_auc - baseline_auc,
                "baseline_auprc": baseline_ap,
                "combined_auprc": combined_ap,
                "delta_auprc": combined_ap - baseline_ap,
                **bootstrap_delta(
                    frame,
                    labels,
                    baseline,
                    combined,
                ),
                "selection_note": (
                    "Feature set and C were selected post hoc on repeated "
                    "out-of-fold performance; deltas are descriptive and not "
                    "confirmatory."
                ),
            }
        )

    delta_table = pd.DataFrame(delta_rows)
    delta_table.to_csv(
        TABLE_DIR / "benchmark_incremental_model_deltas.csv",
        index=False,
    )
    summary = {
        "repetitions": REPEATS,
        "splits": N_SPLITS,
        "workers": WORKERS,
        "prespecified_by_task": {
            row["task"]: {
                "feature_set": row["selected_feature_set"],
                "c_value": row["selected_c_value"],
                "delta_auroc": row["delta_auroc"],
                "delta_auroc_ci_low": row["delta_auroc_ci_low"],
                "delta_auroc_ci_high": row["delta_auroc_ci_high"],
                "delta_auprc": row["delta_auprc"],
                "delta_auprc_ci_low": row["delta_auprc_ci_low"],
                "delta_auprc_ci_high": row["delta_auprc_ci_high"],
            }
            for row in delta_rows
            if row["analysis"] == "prespecified"
        },
        "posthoc_best_by_task": {
            row["task"]: {
                "feature_set": row["selected_feature_set"],
                "c_value": row["selected_c_value"],
                "delta_auroc": row["delta_auroc"],
                "delta_auroc_ci_low": row["delta_auroc_ci_low"],
                "delta_auroc_ci_high": row["delta_auroc_ci_high"],
                "delta_auprc": row["delta_auprc"],
                "delta_auprc_ci_low": row["delta_auprc_ci_low"],
                "delta_auprc_ci_high": row["delta_auprc_ci_high"],
            }
            for row in delta_rows
            if row["analysis"] == "posthoc_best"
        },
        "interpretation": (
            "The prespecified model tests whether molecular context adds value "
            "beyond PIP. The post hoc best-model results are descriptive only."
        ),
    }
    (TABLE_DIR / "benchmark_incremental_model_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    print(metrics_table.sort_values(
        ["task", "mean_auroc"],
        ascending=[True, False],
    ).to_string(index=False))
    print(delta_table.to_string(index=False))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
