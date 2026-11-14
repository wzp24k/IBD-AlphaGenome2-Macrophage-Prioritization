# IBD Non-coding Variant Prioritization with AlphaGenome 2

This project implements a reproducible pipeline for prioritizing inflammatory
bowel disease (IBD) non-coding credible-set variants by gut
macrophage-associated regulatory effects.

## Analysis frame

- Primary discovery study: `GCST90292538` (combined IBD).
- Locus-description companions: disease-specific Crohn disease and ulcerative
  colitis studies plus the overlapping 2017 phenotype analyses. These are not
  counted as independent replication.
- External genetic support: UK Biobank WGS 2025 and VA Million Veteran Program
  2024, each counted once as a non-overlapping cohort cluster.
- Candidate set: 95% credible-set variants with posterior inclusion probability
  (PIP) >= 0.1, excluding coding and canonical splice consequences.
- Primary AlphaGenome 2 condition: colon macrophage
  (`UBERON:0001155_CL:0000235`).
- Controls: colon T cell (`UBERON:0001155_CL:0000084`) and colon epithelial
  cell (`UBERON:0001155_CL:0000066`).
- Replication conditions: matching ileum macrophage, T-cell, and epithelial
  conditions.
- Positive control: rs2836882 (`chr21:39094644:G>A`) with ETS2 as the target
  gene.

## Score interpretation

The AlphaGenome Atlas `quantile_score` is treated as the rank of a variant's
signed RNA-seq effect against the genome-wide distribution of common variants
(MAF > 0.01 in gnomAD v3), using the June 2026 genome-wide calibration. A
signed RNA-seq quantile spans [-1, 1]. Live-model, cell-conditioned expression
scores are reported separately as signed log-fold changes because the current
live API response does not expose the Atlas quantile layer.

## Reproduce

```powershell
$py = "work\python312\python.exe"
& $py scripts\01_fetch_credible_sets.py
& $py scripts\02_prepare_variants.py
powershell -ExecutionPolicy Bypass -File .\run_core_pipeline.ps1
```

The API key can be supplied through `ALPHAGENOME_API_KEY` or `API.txt`. The key
is never included in generated outputs.

## Main deliverables

- `manuscript\Manuscript_AlphaGenome2_IBD.docx`
- `manuscript\manuscript.md`
- `Supplementary_Tables.docx`
- `manuscript\Supplementary_Tables.docx`
- `results\IBD_AlphaGenome2_supplementary_tables.xlsx`
- `results\tables\supplementary_table_index.csv`
- `results\supplementary_tables_audit.json`
- `Cover_Letter_Generic.docx`
- `results\tables\ranked_lead_variants_final.csv`
- `results\tables\ranked_variant_gene_candidates_final.csv`
- `results\figures\figure1_analysis_flow.png` through `figure12_prkcb_mechanism.png`
- `results\figures\figure1_analysis_flow.drawio`
- `Figure_1.drawio`
- `results\tables\ranked_lead_variants_active_ibd.csv`
- `results\tables\active_ibd_validation_summary.json`
- `results\tables\external_gwas_block_bootstrap_pathways.csv`
- `results\tables\benchmark_model_metrics.csv`
- `results\tables\benchmark_model_deltas.csv`
- `results\tables\benchmark_calibration_metrics.csv`
- `results\tables\benchmark_cohort_transportability.csv`
- `results\tables\benchmark_heldout_validation.csv`
- `results\tables\external_eqtl_variant_benchmark.csv`
- `results\tables\external_eqtl_top10_triangulation.csv`
- `results\tables\top10_eqtl_coloc_results.csv`
- `results\tables\top10_eqtl_coloc_sensitivity_summary.csv`
- `results\tables\top10_clinical_translation_annotation.csv`
- `results\tables\prkcq_repeated_pool_summary.csv`
- `results\tables\prkcb_mr_estimates.csv`
- `results\tables\prkcb_smr_estimates.csv`
- `results\tables\prkcb_heidi_sensitivity.csv`
- `results\tables\prkcb_susie_coloc_sensitivity.csv`
- `results\tables\prkcb_molecular_qtl_coverage.csv`
- `results\tables\prkcb_signature_reversal.csv`
- `results\tables\prkcb_virtual_pathway_bootstrap_summary.csv`
- `results\tables\prkcb_mr_leave_one_out.csv`
- `results\tables\prkcb_heidi_window_sensitivity.csv`
- `results\tables\prkcb_virtual_shrinkage_sensitivity.csv`
- `results\tables\prkcb_virtual_bootstrap_1000_summary.csv`
- `results\tables\candidate_evidence_summary.json`
- `results\tables\external_experimental_evidence.csv`
- `results\tables\experimental_gtex_prkcb_eqtl.csv`
- `results\tables\experimental_primary_macrophage_crispr_genes.csv`
- `results\tables\experimental_eqtl_catalogue_coverage.csv`
- `results\tables\perturbseq_gse327644_target_knockdown.csv`
- `results\tables\perturbseq_gse327644_pathway_effects.csv`
- `results\tables\perturbseq_umi_threshold_sensitivity.csv`
- `results\tables\perturbseq_leave_one_guide_out.csv`
- `results\tables\perturbseq_leave_one_donor_out.csv`
- `results\tables\perturbseq_mixed_model_effects.csv`
- `results\tables\camk2a_coloc_sensitivity.csv`
- `results\tables\scatac_pdgfb_variant_topic_percentiles.csv`
- `results\tables\gse175685_prkcb_chsy3_gsea.csv`
- `results\tables\external_gut_scatac_variant_annotation.csv`
- `results\EXPERIMENTAL_FOLLOWUP_PLAN.md`
- `results\PERTURBSEQ_AUDIT.md`
- `results\ROBUSTNESS_AND_DIRECTION_AUDIT.md`
- `results\FIGURE_STYLE_AUDIT.md`
- `results\CITATION_CROSSREF_AUDIT.md`
- `results\BENCHMARK_AND_TRANSLATION_AUDIT.md`
- `results\PRKCB_MAP3K8_CLINICAL_FOLLOWUP_AUDIT.md`
- `results\PRKCB_CAUSAL_MECHANISM_AUDIT.md`
- `results\MAIN_NARRATIVE_REBUILD_AUDIT.md`
- `results\tables\prkcb_ibd_trial_results.csv`
- `results\tables\prkcb_map3k8_drug_selectivity.csv`
- `results\tables\all_l2g_pathway_ora_targeted.csv`
- `results\PATHWAY_AUDIT.md`
- `results\OUTLINE_AUDIT.md`

## Method boundary

Atlas `quantile_score` is interpreted as the signed RNA-seq effect rank against
the genome-wide gnomAD v3 common-variant background (MAF > 0.01), using the
June 2026 genome-wide calibration. The live cell-conditioned API does not
return that quantile layer; live results are therefore retained as signed raw
effects and analyzed with an explicit macrophage-versus-control percentile
delta. UK Biobank and MVP provide external human genetic support, not
wet-laboratory evidence that a variant changes gene expression. The PRKCQ
repeated-pool analysis quantifies sampling stability within the available
donors and does not create independent biological replicates. Public
experimental data are integrated separately: GTEx v10 supplies an exact
PRKCB colon eQTL, GSE229472 supplies the ETS2 macrophage MPRA positive control,
GSE327124 supplies primary macrophage CRISPR screens, and GSE175685 supplies
THP-1 overexpression target engagement. GSE327644 supplies primary macrophage
Perturb-seq across resting and LPS states. GSE273194 supplies exact-variant
Crohn intestinal scATAC/cCRE and ABC target annotation. Neutral and inconclusive
experimental results are retained as neutral or inconclusive.

The PRKCB mechanism extension uses 1000 Genomes high-coverage LD, two-sample
MR, an explicitly labeled SMR/HEIDI-style implementation, SuSiE-colocalization,
public signature reversal and a linear-Gaussian virtual knockdown. These
analyses support directional mediation and experimental prioritization. They
do not replace allele-specific MPRA, primary macrophage CRISPRi or a formal
causal experiment.
