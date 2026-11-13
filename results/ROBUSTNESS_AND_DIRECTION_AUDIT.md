# Robustness and direction audit

Audit date: 2026-11-13

## Corrected direction artifact

The first threshold-sensitivity implementation normalized pseudobulk profiles
using only a selected panel of 28 response genes. Because this panel is small
and biologically selected, its row sums are not valid estimates of library
size. The resulting composition shift could change the sign of a guide effect.

All threshold, leave-one-guide, leave-one-donor and donor-aware mixed-model
analyses now normalize using counts from all 18,129 gene features. The
robustness script also fails if leave-one-out or mixed-model input contains a
selected-gene panel rather than full-gene pseudobulk.

## Perturb-seq numerical checks

PDGFB self-knockdown at dominant-guide UMI thresholds 3, 5 and 10:

- -0.635 log2
- -0.636 log2
- -0.642 log2

IL23R self-effect at the same thresholds:

- +0.1070 log2
- +0.1068 log2
- +0.1070 log2

These values are positive and nonsignificant. The CRISPR perturbation did not
produce reproducible IL23R transcript reduction and is therefore classified as
inconclusive, not negative.

PRKCQ LPS IL-23 pathway effects:

- -0.217 log2 at threshold 3
- -0.215 log2 at threshold 5
- -0.216 log2 at threshold 10

PRKCQ LPS IL-6/JAK/STAT3 pathway effects:

- -0.104 log2 at threshold 3
- -0.102 log2 at threshold 5
- -0.100 log2 at threshold 10

## Model specification

A donor random intercept was used when estimable. When random-intercept
variance was zero, non-finite or singular, the model was explicitly refitted
with donor as a fixed effect. The output records model type and random-effect
variance so a fixed-effect fallback cannot be mistaken for a mixed model.

The PDGFB self-knockdown coefficient remained negative under donor blocking.
The PRKCQ LPS IL-23 coefficient was -0.274 with nominal P = 0.014. This is
pathway evidence, not allele-specific evidence, and the small PRKCQ transcript
effect is retained as an inconsistency rather than promoted to strong target
knockdown.

## Orthogonal secondary analyses

CAMK2A formal colocalization used 3,395 harmonized eQTL-GWAS variants. The
primary single-causal-variant approximation favored H2 (0.91) over H4 (0.09).
Across sensitivity priors, H4 never exceeded approximately 0.33 or the
prespecified 0.80 threshold. The manuscript therefore reports absent
colocalization support.

PRKCB overexpression showed positive JAK-STAT GSEA (NES 1.71, FDR 0.028),
IL-2/STAT5 (NES 1.65, FDR 0.038), TNF/NF-kB (NES 1.66, FDR 0.046) and
IL-6/JAK/STAT3 (NES 1.70, FDR 0.050). CHSY3 overexpression showed positive
TNF/NF-kB (NES 2.39, FDR < 0.001) and IL-6/JAK/STAT3 (NES 1.87, FDR 0.0033)
enrichment.

Ribosome, translation, E2F and MYC programs were also prominent. A sample-level
proliferation principal component from 495 E2F, G2-M and MYC genes explained
23.2% of their variance. After gene-wise residualization, all original
directions remained positive, but only CHSY3 TNF/NF-kB retained FDR < 0.05
(residualized NES 1.89, FDR 0.0125). PRKCB JAK-STAT residualized FDR was 0.098;
PRKCB IL-6/JAK/STAT3 and TNF/NF-kB and CHSY3 IL-6/JAK/STAT3 and IL-2/STAT5 no
longer passed FDR < 0.05. GSEA is therefore reported as hypothesis-supporting
and confounded, not as variant-level causality.

The PDGFB variant rs2143178 is in a Crohn DC-macrophage cCRE. Its Topic15
percentile is 98.4%, and Topic15 ranks first for macrophage and DC marker
enrichment. This is exact-variant chromatin annotation, not an allele-specific
perturbation.

## Evidence labels retained

- PRKCB has exact-variant colon eQTL support and gene-level perturbation data.
- CHSY3 has gene-level target engagement and pathway support.
- PDGFB has a cCRE/topic annotation and effective Perturb-seq knockdown.
- PRKCQ has LPS pathway evidence with inconsistent transcript-level knockdown.
- IL23R Perturb-seq is inconclusive because target reduction failed.
- ERGIC1, GCKR, ABCB5, IRF2BP2 and CAMK2A remain experimentally uncovered at
  the variant level, except for CAMK2A's non-colocalizing nominal eQTL.

No non-significant or uncovered result is relabeled as positive.
