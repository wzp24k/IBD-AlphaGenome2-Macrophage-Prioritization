from __future__ import annotations

import json
import math
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT
from remote_tabix import query_region


TABLE_DIR = ROOT / "results" / "tables"
CACHE_DIR = ROOT / "work" / "tabix_cache" / "cross_ancestry_v2"

STUDIES = {
    "IBD_EAS": (
        "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
        "GCST90432001-GCST90433000/GCST90432182/harmonised/"
        "GCST90432182.h.tsv.gz"
    ),
    "IBD_combined": (
        "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
        "GCST90292001-GCST90293000/GCST90292538/harmonised/"
        "GCST90292538.h.tsv.gz"
    ),
    "Crohn": (
        "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
        "GCST90446001-GCST90447000/GCST90446792/harmonised/"
        "GCST90446792.h.tsv.gz"
    ),
    "UC": (
        "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
        "GCST90446001-GCST90447000/GCST90446794/harmonised/"
        "GCST90446794.h.tsv.gz"
    ),
}


def parse_ssf(line: str) -> dict[str, object] | None:
    fields = line.split("\t")
    if len(fields) < 10:
        return None
    try:
        return {
            "chromosome": fields[0].removeprefix("chr"),
            "position": int(fields[1]),
            "effect_allele": fields[2].upper(),
            "other_allele": fields[3].upper(),
            "beta": float(fields[4]),
            "se": float(fields[5]),
            "eaf": float(fields[6]) if fields[6] not in {"", "NA", "."} else np.nan,
            "p": float(fields[7]),
            "rsid": fields[9] if fields[9] not in {"", ".", "NA"} else "",
        }
    except (ValueError, IndexError):
        return None


def query_study(name: str, url: str) -> tuple[str, list[dict[str, object]]]:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_path = CACHE_DIR / f"{name}.csv"
    if cache_path.exists():
        return name, pd.read_csv(cache_path).to_dict("records")
    raw = query_region(
        url,
        url + ".tbi",
        "16",
        23_820_000,
        23_870_000,
        CACHE_DIR,
    )
    rows = [row for line in raw if (row := parse_ssf(line)) is not None]
    frame = pd.DataFrame(rows)
    frame.to_csv(cache_path, index=False)
    return name, rows


def harmonize_outcome(
    instruments: pd.DataFrame,
    outcome_rows: list[dict[str, object]],
) -> pd.DataFrame:
    outcomes = pd.DataFrame(outcome_rows).rename(
        columns={
            "effect_allele": "effect_allele_outcome",
            "other_allele": "other_allele_outcome",
            "rsid": "rsid_outcome",
        }
    )
    merged = instruments.merge(
        outcomes,
        left_on="rsid",
        right_on="rsid_outcome",
        suffixes=("_eqtl", "_outcome"),
    )
    kept = []
    for _, row in merged.iterrows():
        effect = str(row["effect_allele_eqtl"]).upper()
        other = str(row["other_allele_eqtl"]).upper()
        outcome_effect = str(row["effect_allele_outcome"]).upper()
        outcome_other = str(row["other_allele_outcome"]).upper()
        same_effect = (
            outcome_effect == effect
            or outcome_effect.startswith(effect)
            or effect.startswith(outcome_effect)
        )
        same_other = (
            outcome_other == other
            or outcome_other.startswith(other)
            or other.startswith(outcome_other)
        )
        swapped_effect = (
            outcome_effect == other
            or outcome_effect.startswith(other)
            or other.startswith(outcome_effect)
        )
        swapped_other = (
            outcome_other == effect
            or outcome_other.startswith(effect)
            or effect.startswith(outcome_other)
        )
        if same_effect and same_other:
            beta = float(row["beta"])
        elif swapped_effect and swapped_other:
            beta = -float(row["beta"])
        else:
            continue
        kept.append(
            {
                "rsid": row["rsid"],
                "position": int(row["position_eqtl"]),
                "effect_allele": effect,
                "other_allele": other,
                "eqtl_beta": float(row["eql_beta"]),
                "eqtl_se": float(row["eql_se"]),
                "outcome_beta": beta,
                "outcome_se": float(row["se"]),
                "outcome_p": float(row["p"]),
                "outcome_rsid": row.get("rsid_outcome", ""),
            }
        )
    return pd.DataFrame(kept).sort_values("position")


def ivw_analysis(frame: pd.DataFrame) -> dict[str, object]:
    if frame.empty:
        return {
            "nsnp": 0,
            "beta": np.nan,
            "se": np.nan,
            "p": np.nan,
            "q": np.nan,
            "q_p": np.nan,
        }
    bx = frame["eqtl_beta"].to_numpy(float)
    by = frame["outcome_beta"].to_numpy(float)
    sy = frame["outcome_se"].to_numpy(float)
    weights = (bx / sy) ** 2
    beta = float(np.sum(bx * by / sy**2) / np.sum((bx / sy) ** 2))
    se = float(np.sqrt(1 / np.sum((bx / sy) ** 2)))
    z = beta / se
    if len(frame) > 1:
        ratios = by / bx
        ratio_se = sy / np.abs(bx)
        q = float(np.sum(((ratios - beta) / ratio_se) ** 2))
        q_p = float(stats.chi2.sf(q, len(frame) - 1))
    else:
        q = np.nan
        q_p = np.nan
    return {
        "nsnp": len(frame),
        "beta": beta,
        "se": se,
        "p": float(2 * stats.norm.sf(abs(z))),
        "q": q,
        "q_p": q_p,
    }


def target_disease_support(
    targets: list[str],
    study_path: Path,
    disease: str,
) -> pd.DataFrame:
    payload = json.loads(study_path.read_text(encoding="utf-8"))
    rows = []
    for target in targets:
        scores = []
        pips = []
        pvalues = []
        loci = []
        for page in payload["pages"]:
            for locus in page:
                target_rows = [
                    row
                    for row in (locus.get("l2GPredictions") or {}).get("rows", [])
                    if row.get("target", {}).get("approvedSymbol") == target
                ]
                if not target_rows:
                    continue
                scores.append(max(float(row["score"]) for row in target_rows))
                loci.append(locus.get("studyLocusId"))
                locus_rows = (locus.get("locus") or {}).get("rows", [])
                pips.extend(
                    float(row["posteriorProbability"])
                    for row in locus_rows
                    if row.get("posteriorProbability") is not None
                )
                try:
                    pvalues.append(
                        float(locus["pValueMantissa"])
                        * 10 ** int(locus["pValueExponent"])
                    )
                except (KeyError, TypeError, ValueError):
                    pass
        rows.append(
            {
                "target": target,
                "disease": disease,
                "credible_sets": len(set(loci)),
                "max_l2g": max(scores, default=np.nan),
                "max_pip": max(pips, default=np.nan),
                "min_p": min(pvalues, default=np.nan),
                "l2g_ge_0_5": int(sum(score >= 0.5 for score in scores)),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    top = pd.read_csv(TABLE_DIR / "top10_clinical_translation_annotation.csv")
    targets = top.sort_values("priority_rank")["target"].tolist()

    disease_tables = []
    for disease, path in [
        ("Crohn", ROOT / "data/raw/open_targets_v2_GCST90446792.json"),
        ("UC", ROOT / "data/raw/open_targets_v2_GCST90446794.json"),
    ]:
        disease_tables.append(target_disease_support(targets, path, disease))
    disease = pd.concat(disease_tables, ignore_index=True)
    disease_wide = disease.pivot(index="target", columns="disease")
    disease_wide.columns = [
        f"{column}_{suffix.lower()}" for column, suffix in disease_wide.columns
    ]
    disease_wide = disease_wide.reset_index()

    external = pd.read_csv(TABLE_DIR / "external_ibd_l2g_predictions.csv")
    external_counts = (
        external.groupby("target_symbol")
        .agg(
            external_clusters=("cluster", "nunique"),
            external_l2g_max=("l2g_score", "max"),
        )
        .reset_index()
        .rename(columns={"target_symbol": "target"})
    )
    support = (
        top[["priority_rank", "target", "variant"]]
        .merge(disease_wide, on="target", how="left")
        .merge(external_counts, on="target", how="left")
        .fillna({"external_clusters": 0})
    )
    support["cd_uc_any"] = (
        support["credible_sets_crohn"].fillna(0).gt(0)
        | support["credible_sets_uc"].fillna(0).gt(0)
    )
    support["cd_uc_both"] = (
        support["credible_sets_crohn"].fillna(0).gt(0)
        & support["credible_sets_uc"].fillna(0).gt(0)
    )
    support.to_csv(
        TABLE_DIR / "top10_cd_uc_target_support.csv",
        index=False,
    )

    instruments = pd.read_csv(TABLE_DIR / "prkcb_mr_instruments.csv")
    instruments = instruments[
        instruments["clump_r2_threshold"].eq(0.01)
    ].copy()
    instruments["rsid"] = instruments["rsid"].astype(str)
    instruments["position"] = instruments["position"].astype(int)
    instruments["effect_allele_eqtl"] = instruments["variant_key"].str.split(":").str[1]
    instruments["other_allele_eqtl"] = (
        instruments["variant_key"].str.split(":").str[0]
    )
    # variant_key is encoded as position:other:effect; recover the actual alleles.
    allele_parts = instruments["variant_key"].str.split(":", expand=True)
    instruments["other_allele_eqtl"] = allele_parts[1]
    instruments["effect_allele_eqtl"] = allele_parts[2]

    outcomes = {}
    with ThreadPoolExecutor(max_workers=4) as executor:
        futures = {
            name: executor.submit(query_study, name, url)
            for name, url in STUDIES.items()
        }
        for name, future in futures.items():
            out_name, rows = future.result()
            outcomes[out_name] = rows

    harmonized_tables = []
    mr_rows = []
    for study_name, rows in outcomes.items():
        harmonized = harmonize_outcome(instruments, rows)
        harmonized.insert(0, "study", study_name)
        harmonized_tables.append(harmonized)
        result = ivw_analysis(harmonized)
        result.update(
            {
                "study": study_name,
                "ancestry": "East Asian" if study_name == "IBD_EAS" else "Multi/European",
            }
        )
        mr_rows.append(result)
    harmonized_all = pd.concat(harmonized_tables, ignore_index=True)
    harmonized_all.to_csv(
        TABLE_DIR / "prkcb_cross_ancestry_harmonized.csv",
        index=False,
    )
    mr = pd.DataFrame(mr_rows)
    mr.to_csv(
        TABLE_DIR / "prkcb_cross_ancestry_mr.csv",
        index=False,
    )

    summary = {
        "targets": len(targets),
        "cd_supported_targets_l2g_ge_0_5": int(
            disease.loc[
                disease["disease"].eq("Crohn") & disease["max_l2g"].ge(0.5)
            ].shape[0]
        ),
        "uc_supported_targets_l2g_ge_0_5": int(
            disease.loc[
                disease["disease"].eq("UC") & disease["max_l2g"].ge(0.5)
            ].shape[0]
        ),
        "both_disease_targets": int(support["cd_uc_both"].sum()),
        "prkcb_tests": int(len(mr)),
        "prkcb_positive_direction": int((mr["beta"] > 0).sum()),
        "prkcb_nominal_p_lt_0_05": int((mr["p"] < 0.05).sum()),
        "prkcb_eas_beta": float(
            mr.loc[mr["study"].eq("IBD_EAS"), "beta"].iloc[0]
        ),
        "prkcb_eas_p": float(
            mr.loc[mr["study"].eq("IBD_EAS"), "p"].iloc[0]
        ),
    }
    (TABLE_DIR / "cd_uc_cross_ancestry_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))
    print(support.to_string(index=False))
    print(mr.to_string(index=False))


if __name__ == "__main__":
    main()
