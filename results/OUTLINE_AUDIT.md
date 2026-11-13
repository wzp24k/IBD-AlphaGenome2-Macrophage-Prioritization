# Outline and evidence audit

Audit date: 2026-11-12

## Completed outline requirements

- Primary non-coding credible-set analysis: 395 variants from 281 independent
  SuSiE credible sets in GCST90292538.
- Colon macrophage AlphaGenome 2 scoring with colon T-cell and epithelial
  controls, ileum replication conditions, and donor-derived macrophage
  conditioning.
- Healthy and active macrophage conditioning from Oliver 2024, Kong 2023 and
  Garrido-Trigo 2023.
- Literature-anchored ETS2 positive control and IL2RA non-macrophage
  comparator.
- Frequency-matched 71-variant external background.
- All-target L2G pathway reconstruction using 2,328 predictions across 293
  credible-set loci.
- Locus-block bootstrap with 2,000 resamples and Benjamini-Hochberg adjustment
  across four pre-specified pathways.
- Six resampled donor-balanced PRKCQ pools per condition.
- Explicit target-confidence and target-uncertain flags for Tier 1 candidates.
- Exact-variant experimental integration from a GTEx v10 sigmoid-colon eQTL
  for rs7404003 and PRKCB.
- Primary human macrophage MPRA positive-control analysis from GSE229472.
- Primary human myeloid CRISPR screen integration from GSE327124.
- THP-1 myeloid overexpression target-engagement analysis from GSE175685.
- Exact-variant eQTL Catalogue screen of eight monocyte or macrophage datasets
  for the five remaining Tier 1 candidates.
- GSE327644 primary human macrophage Perturb-seq analysis across four donors,
  resting and LPS states.
- Perturb-seq UMI-threshold, leave-one-guide, leave-one-donor and donor-aware
  mixed-model robustness analyses using full gene-library depth.
- Exact-variant Crohn intestinal scATAC/ABC annotation from GSE273194.
- Formal single-causal-variant colocalization sensitivity for CAMK2A.
- PRKCB and CHSY3 transcriptome-wide GSEA after myeloid overexpression.
- Growth-PC1-residualized GSEA sensitivity and a seven-boundary claim
  calibration table.
- scATAC topic specificity analysis for the PDGFB cCRE.
- Reproducible scripts, tables, ten figures, eleven Word tables and final QC.

## Corrected evidence interpretation

- The six initial Open Targets studies are not independent. The three Liu 2023
  phenotype analyses and three de Lange 2017 analyses contain overlapping
  samples.
- External genetic weighting now uses only UKB WGS 2025 and MVP 2024. Multiple
  CD, UC or combined-IBD results from one cohort are counted once.
- PDGFB is present in overlapping 2017 phenotype-companion analyses. Its target
  had no L2G support in either external cohort and is not described as
  externally replicated.
- PRKCQ is classified as attenuated. The six repeated pools quantify sampling
  stability within available donors and are not independent biological donors.
- ABCB5, IRF2BP2 and IL23R are marked target-uncertain because L2G < 0.50.
  IL23R is the only Tier 1 target supported by both external cohorts.
- UKB and MVP provide external human genetic evidence. They are not wet-lab
  experimental validation.
- One of 10 Tier 1 variants has exact variant-level experimental support,
  five targets have gene-level perturbation data, three genes are explicitly
  tested-neutral for cytokine selection, and five remain experimentally
  uncovered.
- Perturb-seq provides effective editing and transcriptome support for PDGFB,
  LPS pathway support for PRKCQ, and an inconclusive editing result for IL23R.
- CAMK2A lacks colocalization support: primary H4 = 0.09 and maximum sensitivity
  H4 is approximately 0.33, below the prespecified 0.80 threshold.

## Final statistical results

| Pathway | AUC difference | 95% bootstrap CI | Bootstrap P | FDR |
|---|---:|---|---:|---:|
| Reactome IL-23 | 0.322 | 0.236 to 0.480 | 0.00050 | 0.0020 |
| IL-2/STAT5 | 0.145 | 0.031 to 0.226 | 0.0085 | 0.017 |
| KEGG JAK-STAT | 0.097 | -0.003 to 0.205 | 0.031 | 0.037 |
| IL-6/JAK/STAT3 | 0.176 | -0.012 to 0.241 | 0.037 | 0.037 |

JAK-STAT and IL-6/JAK/STAT3 have FDR-significant bootstrap tails but narrowly
crossing confidence intervals. The manuscript reports this distinction.

## Remaining author actions

- Confirm target-journal requirements and article-processing-charge choice.
- Confirm funding, ethics and final CRediT contribution statements.
- Perform variant-specific MPRA, CRISPRi/CRISPRa or primary macrophage
  perturbation for the remaining Tier 1 candidates before making causal claims.
