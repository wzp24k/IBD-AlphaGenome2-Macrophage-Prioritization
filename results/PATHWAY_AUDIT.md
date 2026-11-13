# Pathway analysis audit

This file supersedes the earlier pathway audit that used overlapping Open
Targets studies.

## Design

- Primary evidence: all L2G targets from 395 variants across 293 credible-set
  identifiers. Of these, 276 loci had usable target evidence and were used for
  locus-block resampling.
- External clusters: UKB WGS 2025 and MVP 2024, each counted once.
- Gene score: maximum primary priority times L2G, multiplied by
  `1 + log2(1 + external cluster count)`.
- Resampling unit: credible-set locus.
- Replicates: 2,000.
- P value: one-sided bootstrap tail against an AUC difference of zero.
- Multiplicity: Benjamini-Hochberg across four pre-specified pathways.

## Results

| Pathway | Members | External members | AUC difference | 95% CI | P | FDR |
|---|---:|---:|---:|---|---:|---:|
| Reactome IL-23 | 4 | 3 | 0.322 | 0.236 to 0.480 | 0.00050 | 0.0020 |
| IL-2/STAT5 | 19 | 1 | 0.145 | 0.031 to 0.226 | 0.0085 | 0.017 |
| KEGG JAK-STAT | 26 | 5 | 0.097 | -0.003 to 0.205 | 0.031 | 0.037 |
| IL-6/JAK/STAT3 | 15 | 1 | 0.176 | -0.012 to 0.241 | 0.037 | 0.037 |

## Interpretation

All four pathways have positive AUC differences and FDR < 0.05. The IL-23 and
IL-2/STAT5 confidence intervals exclude zero. The JAK-STAT and
IL-6/JAK/STAT3 intervals narrowly include zero, so the manuscript reports
FDR-significant bootstrap tails without claiming strictly positive interval
estimates.

The earlier six-study result is not used because the studies share samples and
therefore cannot provide independent external weighting.

## External experimental support

The GSE327124 primary human macrophage CRISPR TNF screen independently supports
MAP3K8 as a positive regulator of the TNF-high state. MAP3K8 knockout produced a
median log2 high/low effect of -0.385 (permutation P < 0.0001; FDR = 0.0044),
whereas TNF knockout (-1.447) and TNFAIP3 knockout (0.711) confirmed the expected
screen direction.

The same screen included PDGFB, PRKCQ and IL23R, but their four-guide effects
were neutral for TNF selection (P = 0.57, 0.28 and 0.68, respectively). These
results test a cytokine-output phenotype and do not directly test the predicted
target-mRNA effects. They are retained as tested-neutral rather than converted
to positive evidence.

## Source files

- `results/tables/external_ibd_study_summary.csv`
- `results/tables/external_gwas_gene_rankings.csv`
- `results/tables/external_gwas_block_bootstrap_pathways.csv`
- `results/tables/external_gwas_block_bootstrap_summary.json`
- `results/tables/external_experimental_evidence.csv`
- `results/tables/experimental_primary_macrophage_crispr_genes.csv`
