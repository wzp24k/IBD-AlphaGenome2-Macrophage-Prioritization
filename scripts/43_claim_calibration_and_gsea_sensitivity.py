from __future__ import annotations

import json
import sys
from pathlib import Path

import gseapy as gp
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT, configure_proxy


TARGETS = ["PRKCB", "CHSY3"]
GROWTH_TERMS = [
    "E2F Targets",
    "G2-M Checkpoint",
    "Myc Targets V1",
    "Myc Targets V2",
]
PATHWAY_LIBRARIES = [
    "MSigDB_Hallmark_2020",
    "KEGG_2021_Human",
    "Reactome_2022",
]


def deduplicate_frame(frame: pd.DataFrame) -> pd.DataFrame:
    values = frame.groupby(level=0).mean()
    values.index = values.index.astype(str)
    return values


def growth_score(frame: pd.DataFrame, growth_genes: set[str]) -> tuple[np.ndarray, dict]:
    shared = [gene for gene in growth_genes if gene in frame.index]
    if len(shared) < 50:
        raise RuntimeError("Too few proliferation genes for sensitivity analysis.")
    matrix = frame.loc[shared].to_numpy(dtype=float)
    means = matrix.mean(axis=1, keepdims=True)
    standard_deviations = matrix.std(axis=1, ddof=1, keepdims=True)
    valid = standard_deviations[:, 0] > 0
    matrix = matrix[valid]
    means = means[valid]
    standard_deviations = standard_deviations[valid]
    z = (matrix - means) / standard_deviations
    u, singular_values, _ = np.linalg.svd(z.T, full_matrices=False)
    pc1 = u[:, 0] * singular_values[0]
    if np.corrcoef(pc1, z.mean(axis=0))[0, 1] < 0:
        pc1 = -pc1
    variance_explained = float(
        singular_values[0] ** 2 / np.sum(singular_values**2)
    )
    return pc1, {
        "genes_used": len(shared),
        "pc1_variance_explained": variance_explained,
    }


def residualize(frame: pd.DataFrame, growth_pc1: np.ndarray) -> pd.DataFrame:
    design = np.column_stack([np.ones(len(growth_pc1)), growth_pc1])
    values = frame.to_numpy(dtype=float)
    coefficients, _, _, _ = np.linalg.lstsq(design, values.T, rcond=None)
    fitted = (design @ coefficients).T
    residual = values - fitted
    return pd.DataFrame(
        residual,
        index=frame.index,
        columns=frame.columns,
    )


def target_ranking(
    frame: pd.DataFrame,
    target: str,
    empty_columns: list[str],
) -> pd.DataFrame:
    target_columns = [
        column for column in frame.columns if column.startswith(target + "_")
    ]
    target_values = frame[target_columns].to_numpy(dtype=float)
    empty_values = frame[empty_columns].to_numpy(dtype=float)
    statistic, p_value = ttest_ind(
        target_values,
        empty_values,
        axis=1,
        equal_var=False,
        nan_policy="omit",
    )
    ranking = pd.DataFrame(
        {
            "gene": frame.index,
            "stat": statistic,
            "pvalue": p_value,
            "delta": target_values.mean(axis=1)
            - empty_values.mean(axis=1),
        }
    ).dropna()
    ranking = (
        ranking.sort_values("stat", ascending=False)
        .drop_duplicates("gene")
        .set_index("gene")
    )
    ranking.index = ranking.index.astype(str).astype(object)
    return ranking


def claim_table() -> pd.DataFrame:
    rows = [
        {
            "issue": "External experimental coverage",
            "observed_evidence": (
                "1/10 Tier 1 variants has exact-variant experimental support; "
                "5/10 targets have gene-level perturbation data; 5/10 remain "
                "experimentally uncovered."
            ),
            "permitted_claim": (
                "Report coverage explicitly as one exact-variant, five "
                "gene-level and five uncovered candidates."
            ),
            "prohibited_claim": (
                "Do not state that all Top 10 variants were functionally "
                "validated."
            ),
        },
        {
            "issue": "PDGFB external target support",
            "observed_evidence": (
                "PDGFB has high L2G, macrophage specificity, exact cCRE/topic "
                "annotation and effective Perturb-seq knockdown, but no external "
                "UKB or MVP L2G support."
            ),
            "permitted_claim": (
                "Describe multimodal internal and public experimental support "
                "without external target-level replication."
            ),
            "prohibited_claim": (
                "Do not describe the PDGFB target as externally replicated."
            ),
        },
        {
            "issue": "PRKCQ target engagement",
            "observed_evidence": (
                "PRKCQ transcript knockdown is inconsistent; LPS IL-23 and "
                "IL-6/JAK/STAT3 pathway suppression is reproducible across guides."
            ),
            "permitted_claim": (
                "Report attenuated transcript evidence and LPS pathway-level "
                "functional support."
            ),
            "prohibited_claim": (
                "Do not call the transcript-level result strong target engagement."
            ),
        },
        {
            "issue": "IL23R perturbation",
            "observed_evidence": (
                "IL23R self log2FC is approximately +0.107 with FDR 0.85; the "
                "target transcript was not reduced."
            ),
            "permitted_claim": (
                "Classify the Perturb-seq result as inconclusive because editing "
                "was not verified."
            ),
            "prohibited_claim": (
                "Do not call IL23R experimentally negative or experimentally "
                "supported."
            ),
        },
        {
            "issue": "CAMK2A colocalization",
            "observed_evidence": (
                "Primary CAMK2A coloc H4 is 0.09; the maximum sensitivity H4 is "
                "approximately 0.33 and never reaches 0.80."
            ),
            "permitted_claim": (
                "Report no colocalization support and retain CAMK2A as a "
                "model-prioritized, non-positive candidate."
            ),
            "prohibited_claim": (
                "Do not describe CAMK2A as colocalized, replicated or causal."
            ),
        },
        {
            "issue": "GSEA confounding",
            "observed_evidence": (
                "PRKCB and CHSY3 overexpression enrich inflammatory programs, "
                "but E2F, MYC, ribosome, translation and metabolic programs are "
                "also prominent."
            ),
            "permitted_claim": (
                "Report unadjusted enrichment with a growth and metabolic "
                "confounding caveat and growth-residualized sensitivity."
            ),
            "prohibited_claim": (
                "Do not interpret GSEA as variant-level causality or as evidence "
                "that growth effects are absent."
            ),
        },
        {
            "issue": "Bootstrap interval interpretation",
            "observed_evidence": (
                "JAK-STAT FDR is 0.037 and IL-6/JAK/STAT3 FDR is 0.037, but "
                "their 95% bootstrap intervals cross zero."
            ),
            "permitted_claim": (
                "Report FDR-supported one-sided bootstrap tails and disclose "
                "that the interval estimates are not strictly positive."
            ),
            "prohibited_claim": (
                "Do not describe either pathway as a strictly positive "
                "interval-confirmed effect."
            ),
        },
    ]
    return pd.DataFrame(rows)


def main() -> None:
    configure_proxy()
    path = (
        ROOT
        / "data"
        / "external"
        / "experimental"
        / "GSE175685"
        / "GSE175685_THP1_overexpression_processed_20210506.txt.gz"
    )
    table_dir = ROOT / "results" / "tables"
    frame = pd.read_csv(path, sep="\t").set_index("Gene.Symbol")
    frame = deduplicate_frame(frame)
    empty_columns = [
        column for column in frame.columns if column.startswith("Empty_")
    ]
    libraries = {
        library: gp.get_library(library, organism="human")
        for library in PATHWAY_LIBRARIES
    }
    growth_genes = set()
    hallmark = libraries["MSigDB_Hallmark_2020"]
    for term in GROWTH_TERMS:
        growth_genes.update(hallmark[term])
    growth_pc1, growth_metadata = growth_score(frame, growth_genes)
    residuals = residualize(frame, growth_pc1)

    all_results = []
    summary = {
        "growth_pc1_metadata": growth_metadata,
        "targets": {},
    }
    for target in TARGETS:
        target_columns = [
            column
            for column in frame.columns
            if column.startswith(target + "_")
        ]
        target_pc = float(np.mean(growth_pc1[frame.columns.isin(target_columns)]))
        empty_pc = float(np.mean(growth_pc1[frame.columns.isin(empty_columns)]))
        summary["targets"][target] = {
            "growth_pc1_target_mean": target_pc,
            "growth_pc1_empty_mean": empty_pc,
            "growth_pc1_difference": target_pc - empty_pc,
            "pathways": {},
        }
        ranking = target_ranking(residuals, target, empty_columns)
        rank_metric = ranking["stat"].astype(float)
        rank_metric.index = rank_metric.index.astype(object)
        if rank_metric.duplicated(keep=False).any():
            rank_metric = rank_metric + np.random.default_rng(20261113).normal(
                0,
                1e-10,
                len(rank_metric),
            )
        for library, gene_sets in libraries.items():
            result = gp.prerank(
                rnk=rank_metric,
                gene_sets=gene_sets,
                min_size=5,
                max_size=500,
                permutation_num=1000,
                threads=4,
                seed=20261113,
                outdir=None,
            ).res2d
            result["target"] = target
            result["library"] = library
            all_results.append(result)
            summary["targets"][target]["pathways"][library] = (
                result.sort_values("FDR q-val").head(10).to_dict("records")
            )
    adjusted = pd.concat(all_results, ignore_index=True)
    adjusted.to_csv(
        table_dir / "gse175685_growth_residualized_gsea.csv",
        index=False,
    )

    original = pd.read_csv(table_dir / "gse175685_prkcb_chsy3_gsea.csv")
    selected_terms = {
        "JAK-STAT signaling pathway",
        "IL-6/JAK/STAT3 Signaling",
        "TNF-alpha Signaling via NF-kB",
        "IL-2/STAT5 Signaling",
    }
    original_selected = original[
        original["Term"].isin(selected_terms)
    ][["target", "library", "Term", "NES", "FDR q-val"]].rename(
        columns={"NES": "original_NES", "FDR q-val": "original_FDR"}
    )
    adjusted_selected = adjusted[
        adjusted["Term"].isin(selected_terms)
    ][["target", "library", "Term", "NES", "FDR q-val"]].rename(
        columns={"NES": "residualized_NES", "FDR q-val": "residualized_FDR"}
    )
    sensitivity = original_selected.merge(
        adjusted_selected,
        on=["target", "library", "Term"],
        how="inner",
    )
    sensitivity["direction_retained"] = (
        np.sign(sensitivity["original_NES"])
        == np.sign(sensitivity["residualized_NES"])
    )
    sensitivity["fdr_supported_after_residualization"] = (
        sensitivity["residualized_FDR"] < 0.05
    )
    sensitivity.to_csv(
        table_dir / "gse175685_pathway_confounding_sensitivity.csv",
        index=False,
    )
    claims = claim_table()
    claims.to_csv(table_dir / "claim_calibration.csv", index=False)
    summary["selected_pathway_sensitivity"] = sensitivity.to_dict("records")
    summary["claim_calibration_rows"] = len(claims)
    (table_dir / "gse175685_growth_residualized_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(summary["growth_pc1_metadata"], indent=2))
    print(sensitivity.to_string(index=False))


if __name__ == "__main__":
    main()
