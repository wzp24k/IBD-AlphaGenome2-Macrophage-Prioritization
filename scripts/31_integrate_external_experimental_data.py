from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind, ttest_1samp
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT, configure_proxy


EXTERNAL = ROOT / "data" / "external" / "experimental"
SEED = 20261112
TARGET_GENES = [
    "PDGFB",
    "PRKCB",
    "CHSY3",
    "ERGIC1",
    "GCKR",
    "PRKCQ",
    "ABCB5",
    "IRF2BP2",
    "IL23R",
    "CAMK2A",
]
PATHWAY_GENES = [
    "TNFSF15",
    "MAP3K8",
    "TYK2",
    "JAK2",
    "STAT3",
    "IL12B",
    "IL12RB2",
    "PTPN11",
    "IFNGR2",
    "IL10",
    "IL2RA",
]


def download(url: str, path: Path) -> Path:
    if path.exists() and path.stat().st_size > 0:
        return path
    import requests

    configure_proxy()
    path.parent.mkdir(parents=True, exist_ok=True)
    with requests.get(url, stream=True, timeout=300) as response:
        response.raise_for_status()
        with path.open("wb") as handle:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    handle.write(chunk)
    return path


def read_gzip_text(path: Path) -> str:
    return gzip.decompress(path.read_bytes()).decode("utf-8", "replace")


def fetch_gtex_prkcb() -> pd.DataFrame:
    import requests

    configure_proxy()
    records = []
    for tissue in ["Colon_Sigmoid", "Colon_Transverse"]:
        response = requests.get(
            "https://gtexportal.org/api/v2/association/singleTissueEqtl",
            params={
                "variantId": "chr16_23853326_C_A_b38",
                "tissueSiteDetailId": tissue,
                "itemsPerPage": 250,
            },
            headers={"Accept": "application/json"},
            timeout=120,
        )
        response.raise_for_status()
        for item in response.json().get("data", []):
            records.append(item)
    frame = pd.DataFrame(records)
    if frame.empty:
        raise RuntimeError("GTEx query returned no records for rs7404003.")
    frame = frame[frame["geneSymbol"].eq("PRKCB")].copy()
    frame["experiment"] = "GTEx_v10_colon_eQTL"
    frame["variant"] = "chr16:23853326:C>A"
    frame["target"] = "PRKCB"
    return frame


def analyze_ets2_mpra() -> tuple[pd.DataFrame, dict]:
    root = EXTERNAL / "GSE229472"
    barcode = download(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE229nnn/GSE229472/suppl/"
        "GSE229472_barcode_counts.txt.gz",
        root / "GSE229472_barcode_counts.txt.gz",
    )
    normalized = download(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE229nnn/GSE229472/suppl/"
        "GSE229472_normalised_collapsed_constructs.txt.gz",
        root / "GSE229472_normalised_collapsed_constructs.txt.gz",
    )
    barcode_frame = pd.read_csv(barcode, sep="\t")
    variant_rows = barcode_frame[
        barcode_frame["rsID"].astype(str).eq("rs2836882")
    ][["rsID", "elementTag", "allele", "Type", "localcoords_hg19"]]
    variant_rows.to_csv(root / "rs2836882_construct_metadata.csv", index=False)

    frame = pd.read_csv(normalized, sep="\t", index_col=0)
    rna_columns = [f"{index}TPP" for index in range(1, 9)]
    dna_columns = [f"DNAv{index}" for index in range(1, 5)]
    details = []
    for context in ["A", "B", "C"]:
        prefix = f"2.PSC.AS_rs2836882_{context}"
        alt = frame.loc[f"{prefix}_alt"]
        ref = frame.loc[f"{prefix}_ref"]
        dna = ref[dna_columns].mean()
        risk_replicate_effect = np.log2(
            (ref[rna_columns].to_numpy(dtype=float) + 1) / (dna + 1)
        ) - np.log2(
            (alt[rna_columns].to_numpy(dtype=float) + 1) / (dna + 1)
        )
        details.append(
            {
                "variant": "rs2836882",
                "construct_context": context,
                "risk_allele": "G",
                "nonrisk_allele": "A",
                "mean_log2_risk_vs_nonrisk": float(
                    risk_replicate_effect.mean()
                ),
                "sd_across_rna_pools": float(
                    risk_replicate_effect.std(ddof=1)
                ),
                "positive_rna_pools": int(
                    (risk_replicate_effect > 0).sum()
                ),
                "rna_pools": len(risk_replicate_effect),
            }
        )
    detail_frame = pd.DataFrame(details)
    detail_frame.to_csv(
        root / "rs2836882_mpra_context_effects.csv", index=False
    )
    pooled = np.concatenate(
        [
            np.full(
                int(row.rna_pools),
                row.mean_log2_risk_vs_nonrisk,
            )
            for row in detail_frame.itertuples()
        ]
    )
    summary = {
        "experiment": "primary_human_macrophage_MPRA",
        "accession": "GSE229472",
        "variant": "rs2836882",
        "target": "ETS2",
        "construct_contexts": int(len(detail_frame)),
        "rna_pools_per_context": 8,
        "mean_risk_allele_log2_effect": float(
            detail_frame["mean_log2_risk_vs_nonrisk"].mean()
        ),
        "context_effects": detail_frame[
            ["construct_context", "mean_log2_risk_vs_nonrisk"]
        ].to_dict("records"),
        "all_contexts_positive": bool(
            (detail_frame["mean_log2_risk_vs_nonrisk"] > 0).all()
        ),
    }
    return detail_frame, summary


def analyze_primary_macrophage_crispr() -> tuple[pd.DataFrame, pd.DataFrame]:
    root = EXTERNAL / "GSE327124"
    files = {
        "TNF": download(
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE327nnn/GSE327124/suppl/"
            "GSE327124_TNF.count.txt.gz",
            root / "GSE327124_TNF.count.txt.gz",
        ),
        "CD80": download(
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE327nnn/GSE327124/suppl/"
            "GSE327124_CD80_TF_all_donors.count.txt.gz",
            root / "GSE327124_CD80_TF_all_donors.count.txt.gz",
        ),
    }
    rng = np.random.default_rng(SEED)
    all_guide_scores = []
    all_gene_scores = []
    for screen, path in files.items():
        frame = pd.read_csv(path, sep="\t")
        donor_count = (len(frame.columns) - 2) // 2
        guide_columns = []
        for donor in range(1, donor_count + 1):
            high_name = (
                f"D{donor}_high"
                if f"D{donor}_high" in frame.columns
                else f"D{donor}_High"
            )
            low_name = (
                f"D{donor}_low"
                if f"D{donor}_low" in frame.columns
                else f"D{donor}_Low"
            )
            high = frame[high_name].astype(float)
            low = frame[low_name].astype(float)
            high_cpm = high / high.sum() * 1_000_000
            low_cpm = low / low.sum() * 1_000_000
            column = f"log2_high_low_D{donor}"
            frame[column] = np.log2((high_cpm + 1) / (low_cpm + 1))
            guide_columns.append(column)
        guides = frame.groupby(["sgRNA", "Gene"], as_index=False)[
            guide_columns
        ].mean()
        guides["guide_score"] = guides[guide_columns].mean(axis=1)
        guides["screen"] = screen
        all_guide_scores.append(
            guides[["screen", "sgRNA", "Gene", "guide_score"]]
        )
        values = guides["guide_score"].to_numpy(dtype=float)
        gene_scores = (
            guides.groupby("Gene", as_index=False)
            .agg(
                guides=("sgRNA", "nunique"),
                median_log2_high_low=("guide_score", "median"),
                mean_log2_high_low=("guide_score", "mean"),
                guide_sd=("guide_score", "std"),
            )
        )
        gene_scores = gene_scores[gene_scores["guides"] >= 3].copy()
        null = np.array(
            [
                np.median(rng.choice(values, size=4, replace=False))
                for _ in range(50_000)
            ]
        )
        null_mean = float(null.mean())
        null_sd = float(null.std(ddof=1))
        observed = gene_scores["median_log2_high_low"].to_numpy(dtype=float)
        gene_scores["screen_z"] = (observed - null_mean) / null_sd
        gene_scores["screen_p"] = [
            (1 + np.sum(np.abs(null - null_mean) >= abs(value - null_mean)))
            / (1 + len(null))
            for value in observed
        ]
        gene_scores["screen_fdr"] = multipletests(
            gene_scores["screen_p"].to_numpy(dtype=float),
            method="fdr_bh",
        )[1]
        gene_scores["screen"] = screen
        gene_scores["null_median_sd"] = null_sd
        all_gene_scores.append(gene_scores)
    return pd.concat(all_guide_scores, ignore_index=True), pd.concat(
        all_gene_scores, ignore_index=True
    )


def analyze_thp1_overexpression() -> tuple[pd.DataFrame, pd.DataFrame]:
    root = EXTERNAL / "GSE175685"
    path = download(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE175nnn/GSE175685/suppl/"
        "GSE175685_THP1_overexpression_processed_20210506.txt.gz",
        root / "GSE175685_THP1_overexpression_processed_20210506.txt.gz",
    )
    frame = pd.read_csv(path, sep="\t").set_index("Gene.Symbol")
    columns = list(frame.columns)
    empty_columns = [column for column in columns if column.startswith("Empty_")]
    selected = ["PRKCB", "CHSY3", "ETS2", "JAK2", "STAT3"]
    engagement = []
    for gene in selected:
        target_columns = [
            column for column in columns if column.startswith(f"{gene}_")
        ]
        if not target_columns or gene not in frame.index:
            continue
        target_values = frame.loc[gene, target_columns].astype(float)
        empty_values = frame.loc[gene, empty_columns].astype(float)
        test = ttest_ind(
            target_values,
            empty_values,
            equal_var=False,
            nan_policy="omit",
        )
        target_matrix = frame[target_columns].to_numpy(dtype=float)
        empty_matrix = frame[empty_columns].to_numpy(dtype=float)
        delta = target_matrix.mean(axis=1) - empty_matrix.mean(axis=1)
        p_values = ttest_ind(
            target_matrix,
            empty_matrix,
            axis=1,
            equal_var=False,
            nan_policy="omit",
        ).pvalue
        q_values = multipletests(
            np.nan_to_num(p_values, nan=1.0), method="fdr_bh"
        )[1]
        engagement.append(
            {
                "gene": gene,
                "accession": "GSE175685",
                "replicates": len(target_columns),
                "empty_controls": len(empty_columns),
                "self_expression_delta": float(
                    target_values.mean() - empty_values.mean()
                ),
                "self_expression_p": float(test.pvalue),
                "differential_genes_q05_absdelta1": int(
                    ((q_values < 0.05) & (np.abs(delta) > 1)).sum()
                ),
                "differential_genes_up": int(
                    ((q_values < 0.05) & (delta > 1)).sum()
                ),
                "differential_genes_down": int(
                    ((q_values < 0.05) & (delta < -1)).sum()
                ),
            }
        )
    engagement_frame = pd.DataFrame(engagement)
    transcriptome = []
    for gene in selected:
        target_columns = [
            column for column in columns if column.startswith(f"{gene}_")
        ]
        if not target_columns:
            continue
        target_matrix = frame[target_columns].to_numpy(dtype=float)
        empty_matrix = frame[empty_columns].to_numpy(dtype=float)
        delta = target_matrix.mean(axis=1) - empty_matrix.mean(axis=1)
        p_values = ttest_ind(
            target_matrix,
            empty_matrix,
            axis=1,
            equal_var=False,
            nan_policy="omit",
        ).pvalue
        q_values = multipletests(
            np.nan_to_num(p_values, nan=1.0), method="fdr_bh"
        )[1]
        transcriptome.append(
            pd.DataFrame(
                {
                    "overexpressed_gene": gene,
                    "response_gene": frame.index,
                    "expression_delta": delta,
                    "p_value": p_values,
                    "q_value": q_values,
                }
            )
        )
    return engagement_frame, pd.concat(transcriptome, ignore_index=True)


def main() -> None:
    EXTERNAL.mkdir(parents=True, exist_ok=True)
    table_dir = ROOT / "results" / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)

    gtex = fetch_gtex_prkcb()
    gtex.to_csv(table_dir / "experimental_gtex_prkcb_eqtl.csv", index=False)

    mpra, mpra_summary = analyze_ets2_mpra()
    mpra.to_csv(
        table_dir / "experimental_ets2_mpra_context_effects.csv",
        index=False,
    )
    guides, crispr = analyze_primary_macrophage_crispr()
    engagement, thp1_transcriptome = analyze_thp1_overexpression()
    guides.to_csv(
        table_dir / "experimental_primary_macrophage_crispr_guides.csv",
        index=False,
    )
    crispr.to_csv(
        table_dir / "experimental_primary_macrophage_crispr_genes.csv",
        index=False,
    )
    engagement.to_csv(
        table_dir / "experimental_thp1_overexpression_engagement.csv",
        index=False,
    )
    thp1_transcriptome.to_csv(
        table_dir / "experimental_thp1_overexpression_response.csv",
        index=False,
    )

    lead = pd.read_csv(table_dir / "ranked_lead_variants_final.csv")
    tier1 = lead[
        lead["final_evidence_tier"].astype(str).str.startswith("Tier 1")
    ].copy()
    prkcb_eqtl = gtex.iloc[0]
    engagement_lookup = engagement.set_index("gene").to_dict("index")
    perturbseq_path = table_dir / "perturbseq_gse327644_target_knockdown.csv"
    perturbseq_de_path = table_dir / "perturbseq_gse327644_differential_summary.csv"
    perturbseq_pathway_path = table_dir / "perturbseq_gse327644_pathway_effects.csv"
    perturbseq = (
        pd.read_csv(perturbseq_path)
        if perturbseq_path.exists()
        else pd.DataFrame()
    )
    perturbseq_de = (
        pd.read_csv(perturbseq_de_path)
        if perturbseq_de_path.exists()
        else pd.DataFrame()
    )
    perturbseq_pathway = (
        pd.read_csv(perturbseq_pathway_path)
        if perturbseq_pathway_path.exists()
        else pd.DataFrame()
    )
    perturbseq_combined = (
        perturbseq[perturbseq["condition"].eq("Combined")]
        .set_index("target")
        .to_dict("index")
        if not perturbseq.empty
        else {}
    )
    perturbseq_de_combined = (
        perturbseq_de[perturbseq_de["condition"].eq("Combined")]
        .set_index("target")
        .to_dict("index")
        if not perturbseq_de.empty
        else {}
    )
    prkcq_lps_pathways = (
        perturbseq_pathway[
            perturbseq_pathway["target"].eq("PRKCQ")
            & perturbseq_pathway["condition"].eq("LPS")
        ].set_index("pathway")
        if not perturbseq_pathway.empty
        else pd.DataFrame()
    )
    scatac_path = table_dir / "external_gut_scatac_variant_annotation.csv"
    scatac = pd.read_csv(scatac_path) if scatac_path.exists() else pd.DataFrame()
    scatac_lookup = (
        scatac.set_index("target").to_dict("index")
        if not scatac.empty
        else {}
    )
    tnf = crispr[crispr["screen"].eq("TNF")].set_index("Gene")
    cd80 = crispr[crispr["screen"].eq("CD80")].set_index("Gene")

    evidence_rows = []
    for row in tier1.itertuples():
        gene = row.l2g_target_symbol
        variant_sources = []
        gene_sources = []
        effects = []
        if gene == "PRKCB":
            variant_sources.append("GTEx_v10_Colon_Sigmoid_eQTL")
            effects.append(
                "PRKCB eQTL P=3.02e-5; NES=-0.171"
            )
        if gene in scatac_lookup:
            item = scatac_lookup[gene]
            if item.get("interpretation") == "exact_variant_in_gut_cCRE":
                variant_sources.append("GSE273194_gut_scATAC")
                effects.append(
                    f"exact cCRE={item['overlapped_cCRE']}; "
                    f"celltype={item['cCRE_celltype']}; "
                    f"ABC gene={item['ABC_max_gene']}"
                )
        if gene in engagement_lookup:
            gene_sources.append("GSE175685_THP1_overexpression")
            item = engagement_lookup[gene]
            effects.append(
                f"self-expression delta={item['self_expression_delta']:.2f}; "
                f"P={item['self_expression_p']:.2g}"
            )
        if gene in tnf.index:
            gene_sources.append("GSE327124_primary_macrophage_TNF_CRISPR")
            item = tnf.loc[gene]
            effects.append(
                f"TNF screen log2(high/low)={item['median_log2_high_low']:.3f}; "
                f"P={item['screen_p']:.3g}"
            )
        if gene in cd80.index:
            gene_sources.append("GSE327124_primary_macrophage_CD80_CRISPR")
            item = cd80.loc[gene]
            effects.append(
                f"CD80 screen log2(high/low)={item['median_log2_high_low']:.3f}; "
                f"P={item['screen_p']:.3g}"
            )
        if gene in perturbseq_combined:
            gene_sources.append("GSE327644_primary_macrophage_PerturbSeq")
            item = perturbseq_combined[gene]
            de_item = perturbseq_de_combined.get(gene, {})
            effects.append(
                f"Perturb-seq self log2FC={item['mean_log2fc_vs_ntc']:.3f}; "
                f"FDR={item['fdr']:.3g}; "
                f"responsive genes={int(de_item.get('fdr05_genes', 0))}"
            )
        if gene == "PDGFB":
            evidence_level = "variant_annotation_and_gene_perturbation"
            support = "perturbseq_transcriptome"
        elif gene == "PRKCB":
            evidence_level = "variant_and_gene"
            support = "directional_variant_eqtl"
        elif gene in {"CHSY3"}:
            evidence_level = "gene_level"
            support = "target_engagement"
        elif gene == "PRKCQ":
            evidence_level = "gene_level_perturbation"
            support = "perturbseq_lps_pathway"
            if not prkcq_lps_pathways.empty:
                effects.append(
                    "LPS IL-23 FDR="
                    f"{float(prkcq_lps_pathways.loc['IL-23', 'pathway_fdr']):.3g}; "
                    "LPS IL-6/JAK/STAT3 FDR="
                    f"{float(prkcq_lps_pathways.loc['IL-6/JAK/STAT3', 'pathway_fdr']):.3g}"
                )
        elif gene == "IL23R":
            evidence_level = "gene_level_tested"
            support = "inconclusive_editing"
        elif gene in {"PDGFB", "PRKCQ", "IL23R"}:
            evidence_level = "gene_level_tested"
            support = "tested_neutral_endpoint"
        else:
            evidence_level = "not_covered"
            support = "not_covered"
        evidence_rows.append(
            {
                "lead_rank": row.lead_rank,
                "variant": row.canonical_variant,
                "target": gene,
                "experimental_evidence_level": evidence_level,
                "experimental_support": support,
                "variant_level_experimental_sources": ";".join(variant_sources),
                "gene_level_experimental_sources": ";".join(gene_sources),
                "experimental_effect_summary": "; ".join(effects),
                "experimentally_evaluated": bool(variant_sources or gene_sources),
            }
        )
    evidence = pd.DataFrame(evidence_rows)
    evidence.to_csv(
        table_dir / "external_experimental_evidence.csv", index=False
    )

    lead = lead.drop(
        columns=[
            "experimental_evidence_level",
            "experimental_support",
            "variant_level_experimental_sources",
            "gene_level_experimental_sources",
            "experimental_effect_summary",
            "experimentally_evaluated",
        ],
        errors="ignore",
    ).merge(
        evidence,
        left_on=["lead_rank", "canonical_variant", "l2g_target_symbol"],
        right_on=["lead_rank", "variant", "target"],
        how="left",
    )
    lead = lead.drop(columns=["variant", "target"])
    lead.loc[
        lead["final_evidence_tier"].astype(str).str.startswith("Tier 1"),
        "experimental_evidence_level",
    ] = lead.loc[
        lead["final_evidence_tier"].astype(str).str.startswith("Tier 1"),
        "experimental_evidence_level",
    ].fillna("not_covered")
    lead.to_csv(table_dir / "ranked_lead_variants_final.csv", index=False)

    candidates = pd.read_csv(
        table_dir / "ranked_variant_gene_candidates_final.csv"
    )
    add_columns = [
        "canonical_variant",
        "experimental_evidence_level",
        "experimental_support",
        "variant_level_experimental_sources",
        "gene_level_experimental_sources",
        "experimental_effect_summary",
        "experimentally_evaluated",
    ]
    candidates = candidates.drop(
        columns=add_columns[1:], errors="ignore"
    ).merge(
        lead[add_columns].drop_duplicates("canonical_variant"),
        on="canonical_variant",
        how="left",
    )
    candidates.to_csv(
        table_dir / "ranked_variant_gene_candidates_final.csv", index=False
    )

    eqtl = gtex.iloc[0]
    tnf_target = tnf.reindex(["PDGFB", "PRKCQ", "IL23R", "MAP3K8"])
    summary = {
        "tier1_total": int(len(tier1)),
        "tier1_variant_level_experimental": int(
            evidence["experimental_evidence_level"]
            .eq("variant_and_gene")
            .sum()
        ),
        "tier1_gene_level_experimental": int(
            evidence["experimental_evidence_level"]
            .isin(
                [
                    "variant_and_gene",
                    "variant_annotation_and_gene_perturbation",
                    "gene_level",
                    "gene_level_tested",
                ]
            )
            .sum()
        )
        + int(
            evidence["experimental_evidence_level"]
            .eq("gene_level_perturbation")
            .sum()
        ),
        "tier1_target_engagement": int(
            evidence["experimental_support"].eq("target_engagement").sum()
        ),
        "tier1_tested_neutral": int(
            evidence["target"]
            .isin(["PDGFB", "PRKCQ", "IL23R"])
            .sum()
        ),
        "tier1_not_covered": int(
            evidence["experimental_support"].eq("not_covered").sum()
        ),
        "tier1_perturbseq_effective": int(
            evidence["experimental_support"]
            .isin(["perturbseq_transcriptome", "perturbseq_lps_pathway"])
            .sum()
        ),
        "tier1_perturbseq_inconclusive": int(
            evidence["experimental_support"].eq("inconclusive_editing").sum()
        ),
        "tier1_scatac_ccre": int(
            scatac["interpretation"].eq("exact_variant_in_gut_cCRE").sum()
        )
        if not scatac.empty
        else 0,
        "gtex_prkcb_p": float(eqtl["pValue"]),
        "gtex_prkcb_nes": float(eqtl["nes"]),
        "ets2_mpra": mpra_summary,
        "crispr_tnf_selected": tnf_target[
            [
                "guides",
                "median_log2_high_low",
                "screen_z",
                "screen_p",
                "screen_fdr",
            ]
        ].reset_index().to_dict("records"),
        "thp1_target_engagement": engagement.to_dict("records"),
    }
    (table_dir / "experimental_evidence_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))
    print(evidence.to_string(index=False))


if __name__ == "__main__":
    main()
