from __future__ import annotations

import json
import os
import sys
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

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT


TABLE_DIR = ROOT / "results" / "tables"
OUTPUT = TABLE_DIR / "expanded_macrophage_benchmark_summary.json"
DATA_OUTPUT = TABLE_DIR / "expanded_macrophage_benchmark_pairs.csv"
SEED = 20261122
REPEATS = 20
SPLITS = 5
WORKERS = max(1, min(os.cpu_count() or 1, 12))

FEATURES = [
    "pip",
    "l2g_score",
    "atlas_global_max_abs_q",
    "live_specificity",
    "custom_macrophage_percentile",
]
PRIMARY_FEATURES = [
    "pip",
    "l2g_score",
    "atlas_global_max_abs_q",
]

MODELS = {
    "PIP-only": ["pip"],
    "PIP+L2G": ["pip", "l2g_score"],
    "PIP+Atlas": ["pip", "atlas_global_max_abs_q"],
    "PIP+L2G+Atlas": PRIMARY_FEATURES,
    "PIP+Live": ["pip", "live_specificity"],
    "PIP+context": FEATURES,
}


def model_pipeline() -> Pipeline:
    return Pipeline(
        [
            ("impute", SimpleImputer(strategy="median")),
            ("scale", StandardScaler()),
            (
                "logistic",
                LogisticRegression(
                    C=0.1,
                    class_weight="balanced",
                    solver="lbfgs",
                    max_iter=3000,
                ),
            ),
        ]
    )


def oof_predictions(
    frame: pd.DataFrame,
    label: str,
    features: list[str],
    seed: int,
) -> np.ndarray:
    y = frame[label].to_numpy(dtype=int)
    probabilities = np.full(len(frame), np.nan)
    splitter = StratifiedGroupKFold(
        n_splits=SPLITS,
        shuffle=True,
        random_state=seed,
    )
    for train, test in splitter.split(
        frame[features],
        y,
        groups=frame["studyLocusId"],
    ):
        model = model_pipeline()
        model.fit(frame.iloc[train][features], y[train])
        probabilities[test] = model.predict_proba(
            frame.iloc[test][features]
        )[:, 1]
    return probabilities


def metrics(labels: np.ndarray, scores: np.ndarray) -> tuple[float, float]:
    return (
        float(roc_auc_score(labels, scores)),
        float(average_precision_score(labels, scores)),
    )


def configuration(payload) -> dict:
    name, features, frame, repeats = payload
    y = frame["active_disease_multicohort_supported"].to_numpy(dtype=int)
    probability_sum = np.zeros(len(frame), dtype=float)
    aucs = []
    aps = []
    for repeat in range(repeats):
        probabilities = oof_predictions(
            frame,
            "active_disease_multicohort_supported",
            features,
            SEED + repeat,
        )
        probability_sum += probabilities
        auc, ap = metrics(y, probabilities)
        aucs.append(auc)
        aps.append(ap)
    consensus = probability_sum / repeats
    consensus_auc, consensus_ap = metrics(y, consensus)
    return {
        "model": name,
        "features": ";".join(features),
        "n": len(frame),
        "positives": int(y.sum()),
        "repetitions": repeats,
        "mean_auroc": float(np.mean(aucs)),
        "mean_auprc": float(np.mean(aps)),
        "consensus_auroc": consensus_auc,
        "consensus_auprc": consensus_ap,
    }


def bootstrap_delta(
    frame: pd.DataFrame,
    baseline: np.ndarray,
    combined: np.ndarray,
    repetitions: int = 2000,
) -> dict[str, float]:
    y = frame["active_disease_multicohort_supported"].to_numpy(dtype=int)
    groups = list(frame.groupby("studyLocusId").indices.values())
    rng = np.random.default_rng(SEED + 2)
    auc_deltas = []
    ap_deltas = []
    for _ in range(repetitions):
        sampled = rng.integers(0, len(groups), len(groups))
        indices = np.concatenate([groups[index] for index in sampled])
        if len(np.unique(y[indices])) < 2:
            continue
        baseline_auc, baseline_ap = metrics(y[indices], baseline[indices])
        combined_auc, combined_ap = metrics(y[indices], combined[indices])
        auc_deltas.append(combined_auc - baseline_auc)
        ap_deltas.append(combined_ap - baseline_ap)
    return {
        "delta_auroc": float(
            np.mean(auc_deltas)
        ),
        "delta_auroc_ci_low": float(np.quantile(auc_deltas, 0.025)),
        "delta_auroc_ci_high": float(np.quantile(auc_deltas, 0.975)),
        "delta_auprc": float(np.mean(ap_deltas)),
        "delta_auprc_ci_low": float(np.quantile(ap_deltas, 0.025)),
        "delta_auprc_ci_high": float(np.quantile(ap_deltas, 0.975)),
        "bootstrap_replicates": len(auc_deltas),
    }


def main() -> None:
    pairs = pd.read_csv(TABLE_DIR / "all_variant_l2g_pairs.csv")
    candidates = pd.read_csv(
        TABLE_DIR / "ranked_variant_gene_candidates_final.csv"
    )
    variant_features = candidates[
        [
            "canonical_variant",
            "pip",
            "atlas_global_max_abs_q",
            "live_specificity",
            "custom_macrophage_percentile",
        ]
    ].drop_duplicates("canonical_variant")
    frame = pairs.merge(
        variant_features,
        on="canonical_variant",
        how="left",
    )
    frame["active_disease_multicohort_supported"] = (
        frame["active_disease_multicohort_supported"]
        .astype(str)
        .str.lower()
        .eq("true")
        .astype(int)
    )
    frame = frame.dropna(
        subset=["studyLocusId", "target_id", "canonical_variant"]
    ).reset_index(drop=True)
    frame.to_csv(DATA_OUTPUT, index=False)

    payloads = [
        (name, features, frame, REPEATS)
        for name, features in MODELS.items()
    ]
    with ProcessPoolExecutor(max_workers=WORKERS) as executor:
        rows = list(executor.map(configuration, payloads))
    metrics_table = pd.DataFrame(rows)
    metrics_table.to_csv(
        TABLE_DIR / "expanded_macrophage_benchmark_metrics.csv",
        index=False,
    )

    pip_predictions = oof_predictions(
        frame,
        "active_disease_multicohort_supported",
        MODELS["PIP-only"],
        SEED,
    )
    combined_predictions = oof_predictions(
        frame,
        "active_disease_multicohort_supported",
        MODELS["PIP+context"],
        SEED,
    )
    delta = bootstrap_delta(frame, pip_predictions, combined_predictions)
    delta["baseline"] = "PIP-only"
    delta["comparison"] = "PIP+context"
    pd.DataFrame([delta]).to_csv(
        TABLE_DIR / "expanded_macrophage_benchmark_delta.csv",
        index=False,
    )
    genetic_predictions = oof_predictions(
        frame,
        "active_disease_multicohort_supported",
        MODELS["PIP+L2G+Atlas"],
        SEED,
    )
    genetic_delta = bootstrap_delta(frame, pip_predictions, genetic_predictions)
    genetic_delta["baseline"] = "PIP-only"
    genetic_delta["comparison"] = "PIP+L2G+Atlas"
    pd.DataFrame([genetic_delta]).to_csv(
        TABLE_DIR / "expanded_macrophage_benchmark_genetic_delta.csv",
        index=False,
    )

    summary = {
        "label": "active_disease_multicohort_supported",
        "n_pairs": len(frame),
        "n_credible_sets": int(frame["studyLocusId"].nunique()),
        "n_targets": int(frame["target_id"].nunique()),
        "positives": int(frame["active_disease_multicohort_supported"].sum()),
        "repetitions": REPEATS,
        "splits": SPLITS,
        "workers": WORKERS,
        "models": metrics_table.to_dict(orient="records"),
        "delta": delta,
        "genetic_delta": genetic_delta,
        "primary_model": "PIP+context",
        "genetic_model": "PIP+L2G+Atlas",
        "model_note": (
            "PIP+context is the prespecified primary model. The genetic/model-only "
            "PIP+L2G+Atlas model is reported as a complementary comparison."
        ),
        "note": (
            "The expanded benchmark uses the active-disease multicohort "
            "macrophage-state label available for all 901 variant-target pairs. "
            "Both the full macrophage-context model and the genetic/model-only "
            "comparator are reported."
        ),
    }
    OUTPUT.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
