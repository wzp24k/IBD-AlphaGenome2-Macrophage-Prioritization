from __future__ import annotations

import gzip
import io
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT


DATA_DIR = ROOT / "data" / "external" / "treatment_response"
TABLE_DIR = ROOT / "results" / "tables"

DATASETS = [
    {
        "id": "GSE14580",
        "platform": "GPL570",
        "disease": "UC",
        "treatment": "infliximab",
    },
    {
        "id": "GSE111761",
        "platform": "GPL13497",
        "disease": "CD",
        "treatment": "anti-TNF",
    },
]

MACROPHAGE_GENES = [
    "CD68",
    "CD14",
    "FCGR3A",
    "AIF1",
    "C1QA",
    "C1QB",
    "C1QC",
    "LYZ",
    "TYROBP",
    "ITGAM",
    "CSF1R",
    "MRC1",
    "MSR1",
    "TLR4",
    "CCL2",
    "CXCL10",
]


def parse_series_matrix(path: Path) -> tuple[list[str], list[str], pd.DataFrame]:
    sample_titles = []
    characteristics = []
    table_lines = []
    in_table = False
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!Sample_title"):
                sample_titles = [
                    value.strip('"')
                    for value in line.rstrip("\n").split("\t")[1:]
                ]
            elif line.startswith("!Sample_characteristics_ch1"):
                characteristics.append(
                    [
                        value.strip('"')
                        for value in line.rstrip("\n").split("\t")[1:]
                    ]
                )
            elif line.startswith("!series_matrix_table_begin"):
                in_table = True
            elif line.startswith("!series_matrix_table_end"):
                break
            elif in_table:
                table_lines.append(line)
    matrix = pd.read_csv(
        io.StringIO("".join(table_lines)),
        sep="\t",
        quotechar='"',
    )
    matrix = matrix.rename(columns={"ID_REF": "probe_id"})
    return sample_titles, characteristics, matrix


def parse_platform(path: Path) -> pd.DataFrame:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    marker = next(
        i for i, line in enumerate(lines)
        if line.lower() == "!platform_table_begin"
    )
    end = next(
        i for i, line in enumerate(lines)
        if line.lower() == "!platform_table_end"
    )
    platform = pd.read_csv(
        io.StringIO("\n".join(lines[marker + 1 : end])),
        sep="\t",
        dtype=str,
    )
    symbol_column = next(
        column for column in platform.columns
        if column.lower().replace("_", " ") in {"gene symbol", "gene_symbol"}
    )
    output = platform[["ID", symbol_column]].rename(
        columns={"ID": "probe_id", symbol_column: "gene_symbol"}
    )
    output["gene_symbol"] = (
        output["gene_symbol"]
        .fillna("")
        .str.split("///")
        .str[0]
        .str.strip()
    )
    return output[output["gene_symbol"].ne("")]


def collapse_to_genes(
    matrix: pd.DataFrame,
    platform: pd.DataFrame,
) -> pd.DataFrame:
    merged = matrix.merge(platform, on="probe_id", how="inner")
    sample_columns = [
        column for column in matrix.columns if column != "probe_id"
    ]
    merged["probe_mean"] = merged[sample_columns].mean(axis=1, skipna=True)
    merged = merged.sort_values("probe_mean", ascending=False).drop_duplicates(
        "gene_symbol"
    )
    return merged.set_index("gene_symbol")[sample_columns].apply(
        pd.to_numeric,
        errors="coerce",
    )


def response_labels(dataset_id: str, titles: list[str], characteristics: list[list[str]]) -> list[bool | None]:
    labels = []
    for index, title in enumerate(titles):
        if dataset_id == "GSE111761":
            lower = title.lower()
            if "non-responder" in lower or "non- responder" in lower:
                labels.append(False)
            elif "responder" in lower:
                labels.append(True)
            else:
                labels.append(None)
            continue
        labels.append(None)
        for row in characteristics:
            if index >= len(row):
                continue
            value = row[index].lower()
            if "response to infliximab:" in value:
                if value.endswith("yes") or "responder" in value:
                    labels[-1] = True
                elif value.endswith("no") or "nonresponder" in value:
                    labels[-1] = False
    return labels


def zscore_rows(frame: pd.DataFrame) -> pd.DataFrame:
    means = frame.mean(axis=1)
    sds = frame.std(axis=1, ddof=1).replace(0, np.nan)
    return frame.sub(means, axis=0).div(sds, axis=0)


def program_score(
    expression: pd.DataFrame,
    genes: list[str],
    minimum_genes: int = 4,
) -> tuple[pd.Series, int]:
    present = [gene for gene in genes if gene in expression.index]
    if len(present) < minimum_genes:
        return pd.Series(np.nan, index=expression.columns), len(present)
    return expression.loc[present].mean(axis=0), len(present)


def smd_analysis(
    values: pd.Series,
    labels: list[bool | None],
) -> dict[str, float]:
    valid = [
        (float(value), bool(label))
        for value, label in zip(values, labels)
        if label is not None and pd.notna(value)
    ]
    responder = np.array([value for value, label in valid if label], dtype=float)
    nonresponder = np.array([value for value, label in valid if not label], dtype=float)
    if len(responder) < 2 or len(nonresponder) < 2:
        return {
            "n_responder": len(responder),
            "n_nonresponder": len(nonresponder),
            "mean_responder": np.nan,
            "mean_nonresponder": np.nan,
            "smd": np.nan,
            "smd_se": np.nan,
            "p": np.nan,
            "auc": np.nan,
        }
    pooled_var = (
        (len(responder) - 1) * responder.var(ddof=1)
        + (len(nonresponder) - 1) * nonresponder.var(ddof=1)
    ) / (len(responder) + len(nonresponder) - 2)
    smd = (responder.mean() - nonresponder.mean()) / np.sqrt(pooled_var)
    smd_se = np.sqrt(
        (len(responder) + len(nonresponder))
        / (len(responder) * len(nonresponder))
        + smd**2 / (2 * (len(responder) + len(nonresponder)))
    )
    p = stats.ttest_ind(responder, nonresponder, equal_var=False).pvalue
    ranks = stats.rankdata(np.concatenate([responder, nonresponder]))
    rank_sum = ranks[: len(responder)].sum()
    auc = (
        rank_sum - len(responder) * (len(responder) + 1) / 2
    ) / (len(responder) * len(nonresponder))
    return {
        "n_responder": len(responder),
        "n_nonresponder": len(nonresponder),
        "mean_responder": float(responder.mean()),
        "mean_nonresponder": float(nonresponder.mean()),
        "smd": float(smd),
        "smd_se": float(smd_se),
        "p": float(p),
        "auc": float(auc),
    }


def random_effect_meta(frame: pd.DataFrame) -> dict[str, float]:
    valid = frame.dropna(subset=["smd", "smd_se"]).copy()
    if valid.empty:
        return {
            "datasets": 0,
            "smd": np.nan,
            "se": np.nan,
            "ci_low": np.nan,
            "ci_high": np.nan,
            "p": np.nan,
            "i2": np.nan,
        }
    weights = 1 / valid["smd_se"] ** 2
    fixed = float((valid["smd"] * weights).sum() / weights.sum())
    q = float((weights * (valid["smd"] - fixed) ** 2).sum())
    df = max(len(valid) - 1, 1)
    c = float(weights.sum() - (weights**2).sum() / weights.sum())
    tau2 = max(0.0, (q - df) / c) if c > 0 else 0.0
    random_weights = 1 / (valid["smd_se"] ** 2 + tau2)
    beta = float((valid["smd"] * random_weights).sum() / random_weights.sum())
    se = float(np.sqrt(1 / random_weights.sum()))
    return {
        "datasets": len(valid),
        "smd": beta,
        "se": se,
        "ci_low": beta - 1.96 * se,
        "ci_high": beta + 1.96 * se,
        "p": float(2 * stats.norm.sf(abs(beta / se))),
        "i2": float(max(0, (q - df) / q) * 100) if q > 0 else 0.0,
    }


def main() -> None:
    virtual = pd.read_csv(TABLE_DIR / "prkcb_virtual_knockdown_genes.csv")
    network_genes = sorted(set(virtual["gene"].astype(str)))

    rows = []
    score_tables = []
    for dataset in DATASETS:
        title, characteristics, matrix = parse_series_matrix(
            DATA_DIR / f'{dataset["id"]}_series_matrix.txt.gz'
        )
        platform = parse_platform(ROOT / "work" / f'{dataset["platform"]}_data.txt')
        expression = collapse_to_genes(matrix, platform)
        labels = response_labels(dataset["id"], title, characteristics)
        selected_columns = [column for column, label in zip(expression.columns, labels) if label is not None]
        selected_labels = [label for label in labels if label is not None]
        expression = expression[selected_columns]
        normalized = zscore_rows(expression)
        prkcb_program, prkcb_n = program_score(normalized, network_genes, 5)
        macrophage_program, macrophage_n = program_score(normalized, MACROPHAGE_GENES, 5)
        features = {
            "PRKCB_expression": normalized.loc["PRKCB"] if "PRKCB" in normalized.index else pd.Series(np.nan, index=normalized.columns),
            "PDGFB_expression": normalized.loc["PDGFB"] if "PDGFB" in normalized.index else pd.Series(np.nan, index=normalized.columns),
            "PRKCB_network_program": prkcb_program,
            "macrophage_marker_program": macrophage_program,
        }
        for feature, values in features.items():
            result = smd_analysis(values, selected_labels)
            rows.append(
                {
                    "dataset": dataset["id"],
                    "disease": dataset["disease"],
                    "treatment": dataset["treatment"],
                    "feature": feature,
                    "genes_in_score": prkcb_n
                    if feature == "PRKCB_network_program"
                    else macrophage_n
                    if feature == "macrophage_marker_program"
                    else 1,
                    **result,
                }
            )
        score_tables.append(
            pd.DataFrame(
                {
                    "sample": normalized.columns,
                    "response": selected_labels,
                    "PRKCB_expression": features["PRKCB_expression"].values,
                    "PDGFB_expression": features["PDGFB_expression"].values,
                    "PRKCB_network_program": features["PRKCB_network_program"].values,
                    "macrophage_marker_program": features["macrophage_marker_program"].values,
                }
            ).assign(dataset=dataset["id"], disease=dataset["disease"], treatment=dataset["treatment"])
        )

    long = pd.DataFrame(rows)
    long.to_csv(
        TABLE_DIR / "treatment_response_program_associations.csv",
        index=False,
    )
    pd.concat(score_tables, ignore_index=True).to_csv(
        TABLE_DIR / "treatment_response_sample_scores.csv",
        index=False,
    )
    meta_rows = []
    for feature, group in long.groupby("feature"):
        meta = random_effect_meta(group)
        meta["feature"] = feature
        meta["datasets_evaluable"] = int(group["smd"].notna().sum())
        meta_rows.append(meta)
    meta = pd.DataFrame(meta_rows)
    meta.to_csv(
        TABLE_DIR / "treatment_response_meta_analysis.csv",
        index=False,
    )
    summary = {
        "datasets": [dataset["id"] for dataset in DATASETS],
        "samples": {row["dataset"]: int(row["n_responder"] + row["n_nonresponder"]) for _, row in long.iterrows() if row["feature"] == "PRKCB_expression"},
        "macrophage_genes_present": int((long.loc[long["feature"].eq("macrophage_marker_program"), "genes_in_score"]).max()),
        "prkcb_network_genes_present": int((long.loc[long["feature"].eq("PRKCB_network_program"), "genes_in_score"]).max()),
        "features_tested": int(long["feature"].nunique()),
        "meta_ci_excluding_zero": int(
            ((meta["ci_low"] > 0) | (meta["ci_high"] < 0)).sum()
        ),
    }
    (TABLE_DIR / "treatment_response_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))
    print(long.to_string(index=False))
    print(meta.to_string(index=False))


if __name__ == "__main__":
    main()
