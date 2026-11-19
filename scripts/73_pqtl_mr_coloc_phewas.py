from __future__ import annotations

import importlib.util
import json
import math
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT
from remote_tabix import query_region


TABLE_DIR = ROOT / "results" / "tables"
OPEN_GWAS_INFO = ROOT / "data" / "external" / "opengwas" / "gwasinfo.json"
CACHE_DIR = ROOT / "work" / "tabix_cache" / "pqtl_translation"
TRIAL_URL = "https://api.opengwas.io/api/gwasinfo/files/trial"

PQTL_DATASETS = {
    "SomaScan_PKC_beta": "prot-a-2370",
    "SomaScan_PKCB2": "prot-c-5475_10_3",
}

CURATED_PQTL = pd.DataFrame(
    [
        {
            "chromosome": "16",
            "position": 23_850_005,
            "rsid": "rs12919956",
            "other_allele": "C",
            "effect_allele": "T",
            "beta": 0.14460348,
            "se": (0.205 - 0.085) / 3.92,
            "p": 2e-6,
            "eaf": 0.80578125,
        }
    ]
)

GWAS_CATALOG = {
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
    "IBD_EAS": (
        "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
        "GCST90432001-GCST90433000/GCST90432182/harmonised/"
        "GCST90432182.h.tsv.gz"
    ),
}

TRAITS = {
    "Coronary artery disease": ["coronary artery disease"],
    "Body mass index": ["body mass index"],
    "Rheumatoid arthritis": ["rheumatoid arthritis"],
}

GRCH37_POSITION = {
    23_850_005: 23_861_326,
    23_850_892: 23_862_213,
    23_853_326: 23_864_647,
    23_836_768: 23_848_089,
}


def _load_coloc_module():
    path = Path(__file__).resolve().parent / "48_top10_coloc_and_negative_controls.py"
    spec = importlib.util.spec_from_file_location("coloc48", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


COLOC = _load_coloc_module()


def parse_gwas_catalog(line: str) -> dict[str, object] | None:
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
            "eaf": float(fields[6]) if fields[6] not in {"", ".", "NA"} else np.nan,
            "p": float(fields[7]),
            "rsid": fields[9] if fields[9] not in {"", ".", "NA"} else "",
        }
    except (ValueError, IndexError):
        return None


def parse_vcf(line: str) -> dict[str, object] | None:
    if line.startswith("#"):
        return None
    fields = line.split("\t")
    if len(fields) < 10:
        return None
    info = {}
    for entry in fields[7].split(";"):
        if "=" in entry:
            key, value = entry.split("=", 1)
            info[key] = value
    fmt = fields[8].split(":")
    values = fields[9].split(":")
    record = dict(zip(fmt, values))
    try:
        return {
            "chromosome": fields[0].removeprefix("chr"),
            "position": int(fields[1]),
            "rsid": fields[2],
            "other_allele": fields[3].upper(),
            "effect_allele": fields[4].upper(),
            "beta": float(record.get("ES", "nan")),
            "se": float(record.get("SE", "nan")),
            "p": 10 ** (-float(record.get("LP", "nan"))),
            "eaf": float(record.get("AF", info.get("AF", "nan"))),
        }
    except (ValueError, TypeError):
        return None


def trial_urls(dataset_id: str) -> tuple[str, str]:
    response = requests.post(
        TRIAL_URL,
        json={"id": dataset_id, "type": "vcf"},
        timeout=60,
    )
    response.raise_for_status()
    urls = response.json()
    if len(urls) < 2:
        raise RuntimeError(f"VCF and index unavailable for {dataset_id}")
    return urls[0], urls[1]


def load_gwas_catalog_regions() -> dict[str, pd.DataFrame]:
    def load_one(name: str, url: str) -> tuple[str, pd.DataFrame]:
        cache_path = CACHE_DIR / f"gwas_catalog_{name}.csv"
        if cache_path.exists():
            return name, pd.read_csv(cache_path)
        lines = query_region(
            url,
            url + ".tbi",
            "16",
            23_300_000,
            23_950_000,
            CACHE_DIR,
        )
        rows = [row for line in lines if (row := parse_gwas_catalog(line)) is not None]
        frame = pd.DataFrame(rows)
        frame.to_csv(cache_path, index=False)
        return name, frame

    with ThreadPoolExecutor(max_workers=4) as executor:
        futures = [
            executor.submit(load_one, name, url)
            for name, url in GWAS_CATALOG.items()
        ]
        return dict(future.result() for future in futures)


def load_pqtl_regions() -> dict[str, pd.DataFrame]:
    def load_one(name: str, dataset_id: str) -> tuple[str, pd.DataFrame]:
        cache_path = CACHE_DIR / f"pqtl_{name}.csv"
        if cache_path.exists():
            return name, pd.read_csv(cache_path)
        data_url, index_url = trial_urls(dataset_id)
        lines = query_region(
            data_url,
            index_url,
            "16",
            23_300_000,
            23_950_000,
            CACHE_DIR,
        )
        rows = [row for line in lines if (row := parse_vcf(line)) is not None]
        frame = pd.DataFrame(rows)
        frame.to_csv(cache_path, index=False)
        return name, frame

    with ThreadPoolExecutor(max_workers=2) as executor:
        futures = [
            executor.submit(load_one, name, dataset_id)
            for name, dataset_id in PQTL_DATASETS.items()
        ]
        return dict(future.result() for future in futures)


def select_instruments(frame: pd.DataFrame) -> pd.DataFrame:
    candidates = frame[
        frame["p"].lt(0.05)
        & frame["position"].between(23_350_000, 23_900_000)
        & frame["beta"].replace(0, np.nan).notna()
        & frame["se"].gt(0)
    ].sort_values("p")
    selected = []
    for _, row in candidates.iterrows():
        if any(abs(int(row["position"]) - int(old["position"])) < 100_000 for old in selected):
            continue
        selected.append(row.to_dict())
    return pd.DataFrame(selected[:5])


def harmonize(
    exposure: pd.DataFrame,
    outcome: pd.DataFrame,
    *,
    exposure_prefix: str = "exposure",
    outcome_prefix: str = "outcome",
    merge_key: str = "position",
) -> pd.DataFrame:
    if merge_key == "rsid":
        left = exposure.rename(columns={"position": "position_exposure"}).drop_duplicates("rsid")
        right = outcome.rename(columns={"position": "position_outcome"}).drop_duplicates("rsid")
        merged = left.merge(
            right,
            on="rsid",
            suffixes=(f"_{exposure_prefix}", f"_{outcome_prefix}"),
        )
        positions = merged["position_exposure"]
    else:
        merged = exposure.merge(
            outcome,
            on="position",
            suffixes=(f"_{exposure_prefix}", f"_{outcome_prefix}"),
        )
        positions = merged["position"]
    rows = []
    for position, (_, row) in zip(positions, merged.iterrows()):
        exposure_effect = row[f"effect_allele_{exposure_prefix}"]
        exposure_other = row[f"other_allele_{exposure_prefix}"]
        outcome_effect = row[f"effect_allele_{outcome_prefix}"]
        outcome_other = row[f"other_allele_{outcome_prefix}"]
        if exposure_effect == outcome_effect and exposure_other == outcome_other:
            outcome_beta = float(row[f"beta_{outcome_prefix}"])
        elif exposure_effect == outcome_other and exposure_other == outcome_effect:
            outcome_beta = -float(row[f"beta_{outcome_prefix}"])
        else:
            continue
        rows.append(
            {
                "rsid": row.get(f"rsid_{exposure_prefix}", ""),
                "position": int(position),
                "effect_allele": exposure_effect,
                "other_allele": exposure_other,
                "exposure_beta": float(row[f"beta_{exposure_prefix}"]),
                "exposure_se": float(row[f"se_{exposure_prefix}"]),
                "exposure_p": float(row[f"p_{exposure_prefix}"]),
                "outcome_beta": outcome_beta,
                "outcome_se": float(row[f"se_{outcome_prefix}"]),
                "outcome_p": float(row[f"p_{outcome_prefix}"]),
            }
        )
    return pd.DataFrame(rows).drop_duplicates("position")


def mr_estimate(frame: pd.DataFrame) -> dict[str, float | int]:
    if frame.empty:
        return {"nsnp": 0, "beta": np.nan, "se": np.nan, "p": np.nan, "q_p": np.nan}
    bx = frame["exposure_beta"].to_numpy(float)
    by = frame["outcome_beta"].to_numpy(float)
    sy = frame["outcome_se"].to_numpy(float)
    beta = float(np.sum(bx * by / sy**2) / np.sum((bx / sy) ** 2))
    se = float(np.sqrt(1 / np.sum((bx / sy) ** 2)))
    p = float(2 * stats.norm.sf(abs(beta / se)))
    if len(frame) > 1:
        ratios = by / bx
        ratio_se = sy / np.abs(bx)
        q = float(np.sum(((ratios - beta) / ratio_se) ** 2))
        q_p = float(stats.chi2.sf(q, len(frame) - 1))
    else:
        q_p = np.nan
    return {"nsnp": len(frame), "beta": beta, "se": se, "p": p, "q_p": q_p}


def run_coloc(
    exposure: pd.DataFrame,
    outcome: pd.DataFrame,
    exposure_name: str,
    outcome_name: str,
) -> dict[str, object]:
    harmonized = harmonize(
        exposure,
        outcome,
        exposure_prefix="exposure",
        outcome_prefix="outcome",
        merge_key="rsid",
    )
    # Convert the common harmonized columns into the naming expected by coloc_abf.
    coloc_frame = pd.DataFrame(
        {
            "eql_beta": harmonized["exposure_beta"],
            "eql_se": harmonized["exposure_se"],
            "gwas_beta": harmonized["outcome_beta"],
            "gwas_se": harmonized["outcome_se"],
        }
    )
    if len(coloc_frame) < 50:
        return {
            "exposure_dataset": exposure_name,
            "outcome": outcome_name,
            "variants": len(coloc_frame),
            "H4": np.nan,
            "status": "insufficient_common_variants",
        }
    posteriors, _ = COLOC.coloc_abf(
        coloc_frame,
        p1=1e-4,
        p2=1e-4,
        p12=5e-6,
        sd_prior=0.1,
    )
    return {
        "exposure_dataset": exposure_name,
        "outcome": outcome_name,
        "variants": len(coloc_frame),
        **posteriors,
        "status": "evaluated",
    }


def choose_phewas_datasets(info: dict) -> dict[str, dict[str, object]]:
    datasets = info["datasets_compressed"]
    chosen = {}
    for label, patterns in TRAITS.items():
        candidates = []
        for dataset_id, values in datasets.items():
            record = values[0]
            trait = str(record[0]).lower()
            if not any(pattern in trait for pattern in patterns):
                continue
            sample_size = record[10] if len(record) > 10 else None
            if not isinstance(sample_size, (int, float)) or not np.isfinite(sample_size):
                continue
            if trait.startswith("response") or "treatment" in trait:
                continue
            candidates.append((float(sample_size), dataset_id, record))
        if candidates:
            sample_size, dataset_id, record = max(candidates)
            chosen[label] = {
                "dataset_id": dataset_id,
                "trait": record[0],
                "sample_size": sample_size,
            }
    return chosen


def query_phewas_outcome(
    label: str,
    dataset_id: str,
    positions: tuple[int, int],
) -> tuple[str, pd.DataFrame]:
    cache_path = CACHE_DIR / f"phewas_grch37_{label.replace(' ', '_')}.csv"
    if cache_path.exists():
        return label, pd.read_csv(cache_path)
    try:
        data_url, index_url = trial_urls(dataset_id)
        lines = query_region(
            data_url,
            index_url,
            "16",
            positions[0] - 5000,
            positions[1] + 5000,
            CACHE_DIR,
        )
        rows = [row for line in lines if (row := parse_vcf(line)) is not None]
        frame = pd.DataFrame(rows)
    except Exception:
        frame = pd.DataFrame()
    frame.to_csv(cache_path, index=False)
    return label, frame


def main() -> None:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    pqtl_region = load_pqtl_regions()
    pqtl = dict(pqtl_region)
    pqtl["deCODE_curated"] = CURATED_PQTL.copy()
    ibd_outcomes = load_gwas_catalog_regions()

    instrument_tables = []
    mr_rows = []
    coloc_rows = []
    for dataset_name, frame in pqtl.items():
        instruments = select_instruments(frame)
        instruments.insert(0, "pqtl_dataset", dataset_name)
        instrument_tables.append(instruments)
        for outcome_name, outcome in ibd_outcomes.items():
            harmonized = harmonize(instruments, outcome, merge_key="rsid")
            result = mr_estimate(harmonized)
            result.update(
                {
                    "pqtl_dataset": dataset_name,
                    "outcome": outcome_name,
                }
            )
            mr_rows.append(result)
            if dataset_name in pqtl_region:
                coloc_rows.append(
                    run_coloc(frame, outcome, dataset_name, outcome_name)
                )

    instruments_all = pd.concat(instrument_tables, ignore_index=True)
    instruments_all.to_csv(
        TABLE_DIR / "pqtl_cis_instruments.csv",
        index=False,
    )
    mr = pd.DataFrame(mr_rows).sort_values(["pqtl_dataset", "outcome"])
    mr.to_csv(
        TABLE_DIR / "pqtl_cis_mr_ibd.csv",
        index=False,
    )
    coloc = pd.DataFrame(coloc_rows).sort_values(["exposure_dataset", "outcome"])
    coloc["coloc_supported"] = coloc["H4"].ge(0.8)
    coloc.to_csv(
        TABLE_DIR / "pqtl_ibd_colocalization.csv",
        index=False,
    )

    # PheW-MR uses the strongest pQTL dataset with the largest instrument set.
    primary_name = "deCODE_curated"
    primary_instruments = CURATED_PQTL.copy()
    primary_instruments["position"] = primary_instruments["position"].map(
        GRCH37_POSITION
    )
    primary_instruments = primary_instruments.dropna(subset=["position"])
    primary_instruments["position"] = primary_instruments["position"].astype(int)
    info = json.loads(OPEN_GWAS_INFO.read_text(encoding="utf-8"))
    phewas_datasets = choose_phewas_datasets(info)
    min_position = int(primary_instruments["position"].min())
    max_position = int(primary_instruments["position"].max())

    phewas_frames = {}
    with ThreadPoolExecutor(max_workers=6) as executor:
        futures = {
            label: executor.submit(
                query_phewas_outcome,
                label,
                metadata["dataset_id"],
                (min_position, max_position),
            )
            for label, metadata in phewas_datasets.items()
        }
        for label, future in futures.items():
            out_label, frame = future.result()
            phewas_frames[out_label] = frame

    phewas_rows = []
    for label, outcome in phewas_frames.items():
        harmonized = harmonize(primary_instruments, outcome, merge_key="rsid")
        result = mr_estimate(harmonized)
        result.update(
            {
                "phenotype": label,
                "gwas_id": phewas_datasets[label]["dataset_id"],
                "sample_size": phewas_datasets[label]["sample_size"],
                "pqtl_dataset": primary_name,
            }
        )
        phewas_rows.append(result)
    phewas = pd.DataFrame(phewas_rows)
    if not phewas.empty:
        from statsmodels.stats.multitest import multipletests

        valid = phewas["p"].notna()
        phewas.loc[valid, "fdr"] = multipletests(
            phewas.loc[valid, "p"],
            method="fdr_bh",
        )[1]
    phewas.to_csv(
        TABLE_DIR / "prkcb_pqtl_phewas_mr.csv",
        index=False,
    )

    summary = {
        "pqtl_datasets": list(pqtl),
        "pqtl_instruments": {
            name: int(select_instruments(frame).shape[0])
            for name, frame in pqtl.items()
        },
        "ibd_mr_tests": int(len(mr)),
        "ibd_positive_direction": int((mr["beta"] > 0).sum()),
        "ibd_nominal_p_lt_0_05": int((mr["p"] < 0.05).sum()),
        "coloc_tests": int(len(coloc)),
        "coloc_h4_ge_0_8": int(coloc["coloc_supported"].sum()),
        "phewas_phenotypes_queried": int(len(phewas)),
        "phewas_phenotypes_evaluable": int((phewas["nsnp"] > 0).sum()),
        "phewas_nominal_p_lt_0_05": int((phewas["p"] < 0.05).sum()),
        "phewas_fdr_lt_0_05": int((phewas["fdr"] < 0.05).sum()),
        "primary_phewas_pqtl_dataset": primary_name,
    }
    (TABLE_DIR / "pqtl_translation_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))
    print(mr.to_string(index=False))
    print(coloc.to_string(index=False))
    print(phewas.to_string(index=False))


if __name__ == "__main__":
    main()
