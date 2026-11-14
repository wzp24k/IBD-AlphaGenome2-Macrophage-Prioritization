from __future__ import annotations

import json
import io
import math
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from scipy import stats
from statsmodels.regression.linear_model import WLS

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT
from remote_tabix import query_region


SEED = 20261114
PRKCB_RSID = "rs7404003"
PRKCB_VARIANT_ID = "16_23853326_C_A"
PRKCB_GTEX_VARIANT = "chr16_23853326_C_A_b38"
PRKCB_ENSEMBL = "ENSG00000166501"
REGION_CHR = "16"
REGION_START = 23_700_000
REGION_END = 23_950_000
GWAS_N = 398_668
EQTL_N = 191
GWAS_CASE_FRACTION = 45_106 / GWAS_N
VCF_URL = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
    "1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/"
    "1kGP_high_coverage_Illumina.chr16.filtered.SNV_INDEL_SV_"
    "phased_panel.vcf.gz"
)
VCF_INDEX_URL = VCF_URL + ".tbi"
OPEN_TARGETS_URL = "https://api.platform.opentargets.org/api/v4/graphql"
EQL_COLUMNS = [
    "molecular_trait_id",
    "chromosome",
    "position",
    "ref",
    "alt",
    "variant",
    "ma_samples",
    "maf",
    "pvalue",
    "beta",
    "se",
    "type",
    "ac",
    "an",
    "r2",
    "molecular_trait_object_id",
    "gene_id",
    "median_tpm",
    "rsid",
]
EQL_BASE = (
    "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/"
    "QTS000002/QTD000021/QTD000021"
)


def bh_adjust(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    result = np.full(values.shape, np.nan)
    valid = np.isfinite(values)
    if not valid.any():
        return result
    selected = values[valid]
    order = np.argsort(selected)
    ranked = selected[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    restored = np.empty_like(adjusted)
    restored[order] = adjusted
    result[valid] = restored
    return result


def load_panel() -> pd.DataFrame:
    panel = pd.read_csv(
        ROOT
        / "data"
        / "external"
        / "1000G_3202_samples_ped_population.txt",
        sep=r"\s+",
    )
    panel["sample"] = panel["SampleID"].astype(str).str.strip()
    return panel


def fetch_1000g_region(cache_dir: Path) -> Path:
    output = (
        ROOT
        / "data"
        / "external"
        / "1000G"
        / "1kGP_high_coverage_chr16_PRKCB_region.vcf"
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists() and output.stat().st_size > 0:
        return output
    lines = query_region(
        VCF_URL,
        VCF_INDEX_URL,
        REGION_CHR,
        REGION_START,
        REGION_END,
        cache_dir,
    )
    if not lines:
        raise RuntimeError("No 1000 Genomes variants were returned for PRKCB.")
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return output


def fetch_eqtl_alleles(cache_dir: Path) -> pd.DataFrame:
    output = (
        ROOT
        / "data"
        / "external"
        / "PRKCB_BLUEPRINT_monocyte_eqtl_region.csv"
    )
    if output.exists() and output.stat().st_size > 0:
        return pd.read_csv(output)
    lines = query_region(
        EQL_BASE + ".all.tsv.gz",
        EQL_BASE + ".all.tsv.gz.tbi",
        REGION_CHR,
        REGION_START,
        REGION_END,
        cache_dir,
    )
    if not lines:
        raise RuntimeError("No BLUEPRINT monocyte eQTL rows were returned.")
    frame = pd.read_csv(
        io.StringIO("\n".join(lines)),
        sep="\t",
        names=EQL_COLUMNS,
    )
    frame = frame[
        frame["molecular_trait_id"].eq(PRKCB_ENSEMBL)
        & frame["rsid"].notna()
    ].copy()
    for column in ["position", "pvalue", "beta", "se"]:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    frame = (
        frame.dropna(
            subset=["rsid", "position", "ref", "alt", "beta", "se"]
        )
        .drop_duplicates("rsid")
        .sort_values("position")
    )
    frame.to_csv(output, index=False)
    return frame


def parse_1000g(
    vcf_path: Path,
    eur_indices: list[int],
) -> tuple[pd.DataFrame, np.ndarray]:
    records = []
    genotypes = []
    with vcf_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            sample_columns = fields[9:]
            if not records and max(eur_indices) >= len(sample_columns):
                raise RuntimeError(
                    "The 1000G VCF sample order does not match the panel."
                )
            dosage = np.empty(len(eur_indices), dtype=np.float32)
            for output_index, sample_index in enumerate(eur_indices):
                value = fields[9 + int(sample_index)]
                alleles = value.replace("|", "/").split("/")[:2]
                called = [allele for allele in alleles if allele in {"0", "1"}]
                if not called:
                    dosage[output_index] = np.nan
                else:
                    dosage[output_index] = sum(int(allele) for allele in called)
            records.append(
                {
                    "chromosome": fields[0],
                    "position": int(fields[1]),
                    "ref_vcf": fields[3].upper(),
                    "alt_vcf": fields[4].upper(),
                }
            )
            genotypes.append(dosage)
    metadata = pd.DataFrame(records)
    if metadata.empty:
        raise RuntimeError("No 1000G variants were parsed.")
    metadata["variant_key"] = (
        metadata["position"].astype(str)
        + ":"
        + metadata["ref_vcf"]
        + ":"
        + metadata["alt_vcf"]
    )
    keep = ~metadata["variant_key"].duplicated(keep="first")
    metadata = metadata.loc[keep].reset_index(drop=True)
    matrix = np.vstack(genotypes).T[:, keep.to_numpy()].astype(np.float32)
    return metadata, matrix


def clean_genotype_matrix(
    matrix: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    result = matrix.astype(np.float32, copy=True)
    call_rate = np.mean(np.isfinite(result), axis=0)
    keep = call_rate >= 0.95
    result = result[:, keep]
    for column in range(result.shape[1]):
        values = result[:, column]
        valid = np.isfinite(values)
        if valid.any():
            values[~valid] = np.nanmean(values)
        else:
            values[:] = 0
    return result, keep


def allele_frequencies(matrix: np.ndarray) -> np.ndarray:
    return matrix.mean(axis=0) / 2


def ld_matrix(matrix: np.ndarray) -> np.ndarray:
    standard = matrix - matrix.mean(axis=0, keepdims=True)
    scale = standard.std(axis=0, ddof=0, keepdims=True)
    scale[scale == 0] = 1
    standard = standard / scale
    correlation = (standard.T @ standard) / matrix.shape[0]
    return np.clip(correlation, -1, 1)


def ld_clump(
    variants: pd.DataFrame,
    r2: np.ndarray,
    threshold: float,
    p_threshold: float = 5e-8,
) -> list[str]:
    eligible = variants[variants["eql_p"] < p_threshold].sort_values(
        ["eql_p", "gwas_p"]
    )
    selected: list[str] = []
    available = set(variants["rsid"])
    for rsid in eligible["rsid"]:
        if rsid not in available:
            continue
        if not selected:
            selected.append(rsid)
            continue
        max_r2 = np.max(
            r2[
                variants.index[variants["rsid"].eq(rsid)][0],
                variants.index[variants["rsid"].isin(selected)],
            ]
            ** 2
        )
        if max_r2 < threshold:
            selected.append(rsid)
    if PRKCB_RSID in available and PRKCB_RSID not in selected:
        selected.insert(0, PRKCB_RSID)
    return selected


def weighted_median(values: np.ndarray, weights: np.ndarray) -> float:
    order = np.argsort(values)
    values = values[order]
    weights = weights[order]
    cumulative = np.cumsum(weights)
    cutoff = cumulative[-1] / 2
    return float(values[np.searchsorted(cumulative, cutoff)])


def mr_ivw(bx: np.ndarray, by: np.ndarray, sy: np.ndarray) -> dict:
    weights = 1 / sy**2
    denominator = np.sum(weights * bx**2)
    estimate = np.sum(weights * bx * by) / denominator
    fixed_se = math.sqrt(1 / denominator)
    residual = by - estimate * bx
    q = float(np.sum(weights * residual**2))
    df = max(1, len(bx) - 1)
    scale = max(1.0, math.sqrt(q / df))
    random_se = fixed_se * scale
    return {
        "method": "IVW_random_effects",
        "nsnp": len(bx),
        "estimate": float(estimate),
        "se": float(random_se),
        "p": float(2 * stats.norm.sf(abs(estimate / random_se))),
        "q": q,
        "q_df": df,
        "q_p": float(stats.chi2.sf(q, df)),
    }


def mr_weighted_median(
    bx: np.ndarray,
    by: np.ndarray,
    sy: np.ndarray,
    rng: np.random.Generator,
) -> dict:
    ratio = by / bx
    ratio_se = sy / np.abs(bx)
    weights = 1 / ratio_se**2
    estimate = weighted_median(ratio, weights)
    boot = np.empty(5000)
    for index in range(len(boot)):
        boot[index] = weighted_median(
            rng.normal(ratio, ratio_se),
            weights,
        )
    se = float(np.std(boot, ddof=1))
    return {
        "method": "weighted_median",
        "nsnp": len(bx),
        "estimate": float(estimate),
        "se": se,
        "p": float(2 * stats.norm.sf(abs(estimate / se))),
        "q": np.nan,
        "q_df": np.nan,
        "q_p": np.nan,
    }


def mr_egger(bx: np.ndarray, by: np.ndarray, sy: np.ndarray) -> dict:
    if len(bx) < 3:
        return {
            "method": "MR_Egger",
            "nsnp": len(bx),
            "estimate": np.nan,
            "se": np.nan,
            "p": np.nan,
            "q": np.nan,
            "q_df": np.nan,
            "q_p": np.nan,
            "intercept": np.nan,
            "intercept_p": np.nan,
        }
    orientation = np.where(bx < 0, -1, 1)
    bx = bx * orientation
    by = by * orientation
    model = WLS(by, np.column_stack([np.ones(len(bx)), bx]), weights=1 / sy**2)
    fit = model.fit()
    estimate = float(fit.params[1])
    se = float(fit.bse[1])
    return {
        "method": "MR_Egger",
        "nsnp": len(bx),
        "estimate": estimate,
        "se": se,
        "p": float(2 * stats.t.sf(abs(estimate / se), fit.df_resid)),
        "q": np.nan,
        "q_df": np.nan,
        "q_p": np.nan,
        "intercept": float(fit.params[0]),
        "intercept_p": float(fit.pvalues[0]),
    }


def mr_instrument_table(
    variants: pd.DataFrame,
    r2: np.ndarray,
    threshold: float,
) -> pd.DataFrame:
    selected = ld_clump(variants, r2, threshold=threshold)
    frame = variants[variants["rsid"].isin(selected)].copy()
    frame["r2_with_lead"] = [
        r2[
            variants.index[variants["rsid"].eq(rsid)][0],
            variants.index[variants["rsid"].eq(PRKCB_RSID)][0],
        ]
        ** 2
        for rsid in frame["rsid"]
    ]
    frame["f_statistic"] = (frame["eql_beta"] / frame["eql_se"]) ** 2
    frame["ratio"] = frame["gwas_beta"] / frame["eql_beta"]
    frame["ratio_se"] = frame["gwas_se"] / frame["eql_beta"].abs()
    return frame.sort_values("eql_p")


def smr_and_heidi(
    variants: pd.DataFrame,
    r: np.ndarray,
    top_rsid: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    top = variants[variants["rsid"].eq(top_rsid)].iloc[0]
    beta_smr = float(top["gwas_beta"] / top["eql_beta"])
    var_smr = (
        top["gwas_se"] ** 2 / top["eql_beta"] ** 2
        + top["gwas_beta"] ** 2 * top["eql_se"] ** 2 / top["eql_beta"] ** 4
    )
    se_smr = math.sqrt(var_smr)
    chi_smr = beta_smr**2 / var_smr
    smr = pd.DataFrame(
        [
            {
                "method": "SMR_style_single_instrument",
                "variant": top_rsid,
                "eqtl_beta": float(top["eql_beta"]),
                "eqtl_se": float(top["eql_se"]),
                "gwas_beta": float(top["gwas_beta"]),
                "gwas_se": float(top["gwas_se"]),
                "smr_beta": beta_smr,
                "smr_se": se_smr,
                "smr_p": float(stats.chi2.sf(chi_smr, 1)),
            }
        ]
    )
    rows = []
    for label, frame in [
        ("all_ld_matched", variants),
        (
            "ld_pruned_r2_0.1",
            variants[
                variants["rsid"].isin(
                    ld_clump(variants, r, threshold=0.1)
                )
            ],
        ),
        (
            "ld_pruned_r2_0.01",
            variants[
                variants["rsid"].isin(
                    ld_clump(variants, r, threshold=0.01)
                )
            ],
        ),
    ]:
        frame = frame.sort_values("position").reset_index(drop=True)
        if len(frame) < 3:
            continue
        indices = [
            variants.index[variants["rsid"].eq(rsid)][0]
            for rsid in frame["rsid"]
        ]
        sub_r = r[np.ix_(indices, indices)]
        bx = frame["eql_beta"].to_numpy(dtype=float)
        by = frame["gwas_beta"].to_numpy(dtype=float)
        sy = frame["gwas_se"].to_numpy(dtype=float)
        residual = by - beta_smr * bx
        covariance = np.outer(sy, sy) * sub_r
        ridge = max(1e-6, float(np.trace(covariance) / len(frame) * 1e-6))
        covariance = covariance + np.eye(len(frame)) * ridge
        inv_cov = np.linalg.pinv(covariance)
        q = float(residual @ inv_cov @ residual)
        h = np.outer(bx, bx @ inv_cov) / float(bx @ inv_cov @ bx)
        df = max(1.0, float(np.trace(np.eye(len(frame)) - h)))
        rows.append(
            {
                "method": "HEIDI_style_LD_aware_heterogeneity",
                "variant_set": label,
                "nsnp": len(frame),
                "smr_beta_reference": beta_smr,
                "q": q,
                "effective_df": df,
                "p": float(stats.chi2.sf(q, df)),
                "median_abs_residual": float(np.median(np.abs(residual))),
                "heterogeneity_status": (
                    "consistent"
                    if stats.chi2.sf(q, df) >= 0.05
                    else "heterogeneous"
                ),
            }
        )
    return smr, pd.DataFrame(rows)


def steiger_table(instruments: pd.DataFrame) -> pd.DataFrame:
    frame = instruments.copy()
    z_eqtl = frame["eql_beta"] / frame["eql_se"]
    z_gwas = frame["gwas_beta"] / frame["gwas_se"]
    r2_eqtl = z_eqtl**2 / (z_eqtl**2 + EQTL_N)
    r2_gwas = z_gwas**2 / (
        z_gwas**2
        + GWAS_N * GWAS_CASE_FRACTION * (1 - GWAS_CASE_FRACTION)
    )
    frame["r2_exposure"] = r2_eqtl
    frame["r2_outcome"] = r2_gwas
    frame["correct_direction"] = r2_eqtl > r2_gwas
    frame["direction_p"] = [
        stats.norm.sf(
            abs(a - b) / math.sqrt(a * (1 - a) / EQTL_N + b * (1 - b) / GWAS_N)
        )
        for a, b in zip(r2_eqtl, r2_gwas)
    ]
    return frame[
        [
            "rsid",
            "r2_exposure",
            "r2_outcome",
            "correct_direction",
            "direction_p",
        ]
    ]


def query_open_targets(variant_id: str) -> dict:
    query = """
    query VariantEvidence($id: String!) {
      variant(variantId: $id) {
        rsIds
        evidences(size: 500) {
          count
          rows {
            datasourceId
            datatypeId
            target { approvedSymbol }
            studyId
            pValueMantissa
            pValueExponent
            beta
            directionOnTarget
            cellType
          }
        }
      }
    }
    """
    response = requests.post(
        OPEN_TARGETS_URL,
        json={"query": query, "variables": {"id": variant_id}},
        timeout=90,
    )
    response.raise_for_status()
    payload = response.json()
    if payload.get("errors"):
        raise RuntimeError(str(payload["errors"]))
    return payload["data"]["variant"]


def query_gtex_sqtl() -> pd.DataFrame:
    endpoint = (
        "http://gtexportal.org/api/v2/association/singleTissueSqtl"
    )
    tissues = [None, "Colon_Sigmoid", "Colon_Transverse"]
    rows = []
    for tissue in tissues:
        params = {
            "variantId": PRKCB_GTEX_VARIANT,
            "itemsPerPage": 250,
        }
        if tissue:
            params["tissueSiteDetailId"] = tissue
        try:
            response = requests.get(endpoint, params=params, timeout=60)
            response.raise_for_status()
            items = response.json().get("data", [])
            status = "available"
        except Exception as exc:
            items = []
            status = f"query_error:{type(exc).__name__}"
        for item in items:
            rows.append(
                {
                    "source": "GTEx_v10_sQTL",
                    "variant": PRKCB_RSID,
                    "target": item.get("geneSymbol"),
                    "tissue": item.get(
                        "tissueSiteDetailId",
                        tissue or "all_tissues",
                    ),
                    "effect": item.get("nes"),
                    "p_value": item.get("pValue"),
                    "status": status,
                }
            )
        if not items:
            rows.append(
                {
                    "source": "GTEx_v10_sQTL",
                    "variant": PRKCB_RSID,
                    "target": "PRKCB",
                    "tissue": tissue or "all_tissues",
                    "effect": np.nan,
                    "p_value": np.nan,
                    "status": status,
                }
            )
    return pd.DataFrame(rows)


def qtl_coverage(
    gtex_eqtl: pd.DataFrame,
    eqtl_triangulation: pd.DataFrame,
) -> pd.DataFrame:
    open_targets = query_open_targets(PRKCB_VARIANT_ID)
    open_rows = open_targets["evidences"]["rows"]
    datasets = sorted(
        {
            f"{row.get('datasourceId')}::{row.get('datatypeId')}"
            for row in open_rows
        }
    )
    sqtl = query_gtex_sqtl()
    sqtl_hit = sqtl["p_value"].notna().any() if len(sqtl) else False
    gtex_hit = (
        gtex_eqtl["pValue"].lt(0.05).fillna(False).any()
        if len(gtex_eqtl)
        else False
    )
    eqtl_row = eqtl_triangulation[
        eqtl_triangulation["query_group"].eq("top10_matched")
        & eqtl_triangulation["target"].eq("PRKCB")
    ]
    myeloid_hit = (
        eqtl_row["any_nominal"].fillna(False).astype(bool).any()
        if len(eqtl_row)
        else False
    )
    rows = [
        {
            "molecular_layer": "mRNA expression",
            "source": "GTEx_v10_colon_eQTL",
            "exact_variant_coverage": bool(len(gtex_eqtl) > 0),
            "association_supported": bool(gtex_hit),
            "observed_direction": (
                float(gtex_eqtl.iloc[0]["nes"])
                if len(gtex_eqtl)
                else np.nan
            ),
            "interpretation": "colon direction available",
        },
        {
            "molecular_layer": "mRNA expression",
            "source": "eight_myeloid_eQTL_datasets",
            "exact_variant_coverage": bool(
                eqtl_row["datasets_queried"].fillna(0).max() > 0
            ),
            "association_supported": bool(myeloid_hit),
            "observed_direction": "positive_in_all_eight",
            "interpretation": "myeloid eQTL direction",
        },
        {
            "molecular_layer": "splicing",
            "source": "GTEx_v10_sQTL",
            "exact_variant_coverage": bool(
                sqtl["status"].eq("available").any()
            ),
            "association_supported": bool(sqtl_hit),
            "observed_direction": np.nan,
            "interpretation": "no exact-variant sQTL returned",
        },
        {
            "molecular_layer": "protein",
            "source": "Open_Targets_integrated_evidence",
            "exact_variant_coverage": bool(len(open_rows) > 0),
            "association_supported": False,
            "observed_direction": np.nan,
            "interpretation": "no exact-variant pQTL record",
        },
        {
            "molecular_layer": "chromatin accessibility",
            "source": "Open_Targets_integrated_evidence",
            "exact_variant_coverage": any(
                "chromatin" in dataset.lower()
                or "caqtl" in dataset.lower()
                for dataset in datasets
            ),
            "association_supported": False,
            "observed_direction": np.nan,
            "interpretation": "no exact-variant caQTL record",
        },
    ]
    sqtl["observed_direction"] = sqtl["effect"]
    return pd.DataFrame(rows), sqtl


def main() -> None:
    rng = np.random.default_rng(SEED)
    table_dir = ROOT / "results" / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = ROOT / "work" / "tabix_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)

    variant_source = pd.read_csv(
        table_dir / "top10_eqtl_coloc_variants.csv"
    )
    variants = variant_source[variant_source["target"].eq("PRKCB")].copy()
    variants = variants.drop_duplicates("rsid").sort_values("position")
    if PRKCB_RSID not in set(variants["rsid"]):
        raise RuntimeError("PRKCB lead variant is absent from coloc input.")

    panel = load_panel()
    eur_indices = panel.index[
        panel["Superpopulation"].eq("EUR")
    ].tolist()
    if len(panel) != 3202:
        raise RuntimeError("The 1000G panel does not contain 3,202 samples.")
    print(f"Fetching 1000G PRKCB region for {len(eur_indices)} EUR samples")
    vcf_path = fetch_1000g_region(cache_dir)
    metadata, genotype = parse_1000g(
        vcf_path,
        eur_indices,
    )
    eqtl_alleles = fetch_eqtl_alleles(cache_dir)
    eqtl_alleles["variant_key"] = (
        eqtl_alleles["position"].astype(str)
        + ":"
        + eqtl_alleles["ref"].astype(str).str.upper()
        + ":"
        + eqtl_alleles["alt"].astype(str).str.upper()
    )
    genotype, call_rate_keep = clean_genotype_matrix(genotype)
    metadata = metadata.iloc[call_rate_keep].reset_index(drop=True)
    frequency = genotype.mean(axis=0) / 2
    maf_keep = (frequency >= 0.01) & (frequency <= 0.99)
    genotype = genotype[:, maf_keep]
    metadata = metadata.iloc[maf_keep].reset_index(drop=True)
    r_all = ld_matrix(genotype)
    ld_lookup = metadata[["variant_key"]].copy()
    ld_lookup["maf_eur"] = np.minimum(
        allele_frequencies(genotype),
        1 - allele_frequencies(genotype),
    )

    analysis = variants.merge(
        eqtl_alleles[["rsid", "variant_key"]],
        on="rsid",
        how="inner",
    ).merge(
        ld_lookup,
        on="variant_key",
        how="inner",
    )
    analysis = (
        analysis.drop_duplicates("rsid")
        .sort_values("position")
        .reset_index(drop=True)
    )
    analysis["z_gwas"] = analysis["gwas_beta"] / analysis["gwas_se"]
    analysis["z_eqtl"] = analysis["eql_beta"] / analysis["eql_se"]
    indices = [
        metadata.index[metadata["variant_key"].eq(key)][0]
        for key in analysis["variant_key"]
    ]
    r = r_all[np.ix_(indices, indices)]
    analysis.to_csv(table_dir / "prkcb_ld_analysis_variants.csv", index=False)
    pd.DataFrame(r, index=analysis["rsid"], columns=analysis["rsid"]).to_csv(
        table_dir / "prkcb_ld_matrix.csv"
    )
    ld_long = []
    lead_index = analysis.index[analysis["rsid"].eq(PRKCB_RSID)][0]
    for index, row in analysis.iterrows():
        ld_long.append(
            {
                "rsid": row["rsid"],
                "position": row["position"],
                "r2_with_rs7404003": float(
                    r[index, lead_index] ** 2
                ),
                "d_prime_with_rs7404003": float(
                    abs(r[index, lead_index])
                ),
            }
        )
    pd.DataFrame(ld_long).to_csv(
        table_dir / "prkcb_ld_with_lead.csv",
        index=False,
    )

    instrument_frames = []
    mr_rows = []
    for threshold in [0.001, 0.01, 0.1]:
        frame = mr_instrument_table(analysis, r, threshold=threshold)
        frame["clump_r2_threshold"] = threshold
        instrument_frames.append(frame)
        bx = frame["eql_beta"].to_numpy(dtype=float)
        by = frame["gwas_beta"].to_numpy(dtype=float)
        sy = frame["gwas_se"].to_numpy(dtype=float)
        for result in [
            mr_ivw(bx, by, sy),
            mr_weighted_median(bx, by, sy, rng),
            mr_egger(bx, by, sy),
        ]:
            result["clump_r2_threshold"] = threshold
            result["mean_f_statistic"] = float(
                frame["f_statistic"].mean()
            )
            mr_rows.append(result)
    instruments = pd.concat(instrument_frames, ignore_index=True)
    instruments.to_csv(
        table_dir / "prkcb_mr_instruments.csv",
        index=False,
    )
    mr_estimates = pd.DataFrame(mr_rows)
    mr_estimates.to_csv(
        table_dir / "prkcb_mr_estimates.csv",
        index=False,
    )

    primary = pd.concat(
        [
            frame
            for frame in instrument_frames
            if float(frame["clump_r2_threshold"].iloc[0]) == 0.01
        ],
        ignore_index=True,
    )
    steiger = steiger_table(primary)
    steiger.to_csv(
        table_dir / "prkcb_mr_steiger.csv",
        index=False,
    )
    smr, heidi = smr_and_heidi(analysis, r, PRKCB_RSID)
    smr.to_csv(table_dir / "prkcb_smr_estimates.csv", index=False)
    heidi.to_csv(table_dir / "prkcb_heidi_sensitivity.csv", index=False)

    gtex_eqtl = pd.read_csv(
        table_dir / "experimental_gtex_prkcb_eqtl.csv"
    )
    triangulation = pd.read_csv(
        table_dir / "external_eqtl_top10_triangulation.csv"
    )
    coverage, sqtl = qtl_coverage(gtex_eqtl, triangulation)
    coverage.to_csv(
        table_dir / "prkcb_molecular_qtl_coverage.csv",
        index=False,
    )
    sqtl.to_csv(table_dir / "prkcb_gtex_sqtl_lookup.csv", index=False)

    primary_mr = mr_estimates[
        mr_estimates["clump_r2_threshold"].eq(0.01)
    ].set_index("method")
    summary = {
        "prkcb_variant": "chr16:23853326:C>A",
        "prkcb_rsid": PRKCB_RSID,
        "gwas_sample_size": GWAS_N,
        "eqtl_sample_size": EQTL_N,
        "eur_ld_samples": len(eur_indices),
        "ld_matched_variants": len(analysis),
        "lead_variant_ld_available": True,
        "primary_clump_r2": 0.01,
        "primary_instruments": int(
            primary_mr["nsnp"].max()
        ),
        "mean_f_statistic": float(
            primary_mr["mean_f_statistic"].max()
        ),
        "ivw_beta": float(primary_mr.loc["IVW_random_effects", "estimate"]),
        "ivw_se": float(primary_mr.loc["IVW_random_effects", "se"]),
        "ivw_p": float(primary_mr.loc["IVW_random_effects", "p"]),
        "weighted_median_beta": float(
            primary_mr.loc["weighted_median", "estimate"]
        ),
        "weighted_median_p": float(
            primary_mr.loc["weighted_median", "p"]
        ),
        "mr_egger_beta": float(
            primary_mr.loc["MR_Egger", "estimate"]
        ),
        "mr_egger_p": float(primary_mr.loc["MR_Egger", "p"]),
        "mr_egger_intercept_p": float(
            primary_mr.loc["MR_Egger", "intercept_p"]
        ),
        "steiger_correct_fraction": float(
            steiger["correct_direction"].mean()
        ),
        "smr_beta": float(smr.iloc[0]["smr_beta"]),
        "smr_p": float(smr.iloc[0]["smr_p"]),
        "heidi_style_primary_p": float(
            heidi.loc[
                heidi["variant_set"].eq("all_ld_matched"),
                "p",
            ].iloc[0]
        ),
        "heidi_style_pruned_p": float(
            heidi.loc[
                heidi["variant_set"].eq("ld_pruned_r2_0.1"),
                "p",
            ].iloc[0]
            if heidi["variant_set"].eq("ld_pruned_r2_0.1").any()
            else heidi.loc[
                heidi["variant_set"].eq("ld_pruned_r2_0.01"),
                "p",
            ].iloc[0]
        ),
        "qtl_layers_with_exact_variant_support": int(
            coverage["association_supported"].sum()
        ),
        "pqtl_exact_variant_support": False,
        "sqtl_exact_variant_support": False,
        "caqtl_exact_variant_support": False,
        "method_boundary": (
            "SMR/HEIDI and MR were implemented here as LD-aware "
            "summary-statistic analyses, not official SMR software. "
            "They support a directional mediation hypothesis but do not "
            "replace allele-specific perturbation or experimental causality."
        ),
    }
    (table_dir / "prkcb_causal_mechanism_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))
    print(mr_estimates.to_string(index=False))
    print(coverage.to_string(index=False))


if __name__ == "__main__":
    main()
