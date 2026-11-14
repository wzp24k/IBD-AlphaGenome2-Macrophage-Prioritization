# PRKCB causal-mechanism audit

## Scope

The analysis adds four computational mechanism layers for PRKCB rs7404003
(chr16:23853326:C>A):

- LD-aware two-sample MR and an SMR/HEIDI-style summary-statistic analysis.
- SuSiE fine-mapping and colocalization with a 1000 Genomes LD reference.
- Exact-variant coverage audit for eQTL, sQTL, pQTL and caQTL layers.
- Public overexpression-signature reversal and state-specific virtual knockdown.

These analyses do not include a newly generated allele-specific MPRA,
CRISPRi or primary macrophage experiment.

## LD reference

- 1000 Genomes high-coverage chr16 VCF.
- 633 European samples.
- GRCh38 positions.
- 877 variants matched across BLUEPRINT monocyte eQTL, IBD GWAS and LD.
- LD clumping r2 thresholds: 0.001, 0.01 and 0.10.

## Two-sample MR

At r2 < 0.01, three approximately independent instruments were retained.

- Mean F statistic: 418.8.
- IVW beta: 0.0517 (SE 0.00595; P = 3.25e-18).
- Weighted median beta: 0.0514 (SE 0.00728; P = 1.63e-12).
- MR-Egger beta: -0.178 (P = 0.032).
- MR-Egger intercept P = 0.025.
- Correct Steiger-style direction: 3/3 instruments.
- Leave-one-instrument-out IVW and weighted median: positive for all omissions.

The IVW and weighted-median directions agree. The Egger result and nonzero
intercept show that a single pleiotropy-free causal model is not established.

## SMR/HEIDI-style analysis

- Single-variant SMR beta: 0.0514.
- SMR P = 3.43e-9.
- Full-region LD-aware heterogeneity P < 1e-300.
- LD-pruned three-instrument heterogeneity P = 0.368.
- Local HEIDI windows of 50, 100, 150 and 250 kb remained heterogeneous.

The heterogeneity pattern supports a directional signal in a restricted
subset and argues against applying one effect size uniformly across the whole
regional signal. This is a summary-statistic implementation and is not output
from the official SMR executable.

## SuSiE-colocalization

- SuSiE L values: 1, 3 and 5.
- Colocalization priors p12: 1e-6, 5e-6 and 1e-5.
- Primary setting: L = 5, p12 = 5e-6.
- Primary maximum shared-credible-set H4: 1.0.
- Maximum H4 across sensitivity settings: 1.0.

This supports a shared causal-variant hypothesis for matched credible sets.
It does not resolve all regional signals into one causal variant.

## Molecular QTL coverage

- Exact-variant mRNA eQTL: covered in GTEx colon and eight myeloid datasets.
- Exact-variant sQTL: no record returned.
- Exact-variant pQTL: no integrated record returned.
- Exact-variant caQTL: no integrated record returned.

Missing sQTL, pQTL and caQTL records are resource non-coverage. They are not
treated as negative associations.

## Public perturbation and virtual knockdown

- Public overexpression signature reversed disease direction in 2/3 cohorts.
- Median correlation: 0.012.
- Median permutation P: 0.096.
- State-specific virtual network: 50 shared genes.
- Donor-blocked bootstrap: 200 replicates.
- Extended donor-blocked bootstrap: 1,000 replicates.
- Active-macrophage JAK/STAT virtual effect: -0.0246
  (95% CI -0.0355 to -0.0149).
- Extended active-macrophage JAK/STAT bootstrap: -0.0241
  (95% CI -0.0356 to -0.0126).
- JAK/STAT direction remained negative across shrinkage 0.02-0.30.
- Active-macrophage TNF/NF-kB virtual effect: -0.0104
  (95% CI -0.0258 to 0.0018).

The virtual model is a conditional Gaussian co-expression prediction. It is
not CellOracle and does not establish target engagement.

## Evidence boundary

Supported:

- The PRKCB risk allele is directionally associated with higher monocyte
  expression.
- Multiple computational causal frameworks prioritize PRKCB.
- Matched fine-mapped signals support a shared credible set.
- The virtual active-macrophage model predicts lower JAK/STAT activity after
  PRKCB knockdown.

Not supported:

- One uniform causal effect across every variant in the region.
- Allele-specific causality from rs7404003.
- Protein-level or chromatin-level mediation.
- A clinical treatment recommendation.
