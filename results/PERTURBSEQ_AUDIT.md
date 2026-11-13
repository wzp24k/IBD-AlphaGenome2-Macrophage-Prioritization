# GSE327644 Perturb-seq audit

Audit date: 2026-11-13

## Data and guide coverage

The GSE327644 archive contains eight filtered Cell Ranger H5 matrices from four
donors in resting and LPS-stimulated macrophages. Each matrix has 18129 gene
features and 2925 guide-capture features. The library contains four guides for
PDGFB, PRKCQ, IL23R, TNFSF15, MAP3K8, TNF and TNFAIP3.

## Assignment and pseudobulk

Cells were assigned only when the dominant guide had at least five UMIs and the
second guide had at most half the dominant count. Control guides were combined
into a matched NTC pool for each donor and state. A total of 124959 cells passed
assignment. Gene counts were aggregated by guide and sample, converted to log2
CPM using all 18,129 gene features, and compared with matched NTC profiles.

## Results

| Target | Combined self log2FC | FDR | FDR < 0.05, combined |
|---|---:|---:|---:|
| PDGFB | -0.636 | 0.0088 | 460 |
| PRKCQ | -0.104 | 0.83 | 548 |
| IL23R | +0.107 | 0.85 | 575 |
| TNFSF15 | -0.685 | 0.15 | 835 |
| MAP3K8 | -0.187 | 0.099 | 556 |
| TNF | -0.317 | 0.67 | 762 |
| TNFAIP3 | +0.319 | 0.078 | 2291 |

All four PRKCQ guides reduced the LPS IL-23 program (empirical FDR = 0.0021)
and the IL-6/JAK/STAT3 program (FDR = 0.020). PDGFB showed effective transcript
knockdown and a measurable transcriptome response. IL23R editing was
inconclusive because the target transcript was not reduced.

## Robustness and corrected direction artifact

An initial threshold-sensitivity implementation normalized pseudobulk counts
using only a selected 28-gene response panel. That introduced a compositional
direction artifact. The analysis now uses total counts across all gene features
for library-size normalization, including the threshold, leave-one-out and
mixed-model analyses.

Across dominant-guide thresholds 3, 5 and 10, PDGFB self-knockdown was -0.635,
-0.636 and -0.642 log2. The IL23R self-effect remained positive and
nonsignificant at approximately +0.107. PRKCQ LPS IL-23 effects were -0.217,
-0.215 and -0.216 log2. Leave-one-guide and leave-one-donor direction
concordance was 92-100% for selected PDGFB and PRKCQ metrics.

For 56 target/condition/pathway combinations, donor-aware models were fitted.
When the donor random-intercept variance was zero, non-finite or singular, the
model was refitted with donor as a fixed effect. PDGFB self-knockdown remained
negative and significant. PRKCQ LPS IL-23 had coefficient -0.274 and nominal
P = 0.014 under the donor-blocked specification.

## Interpretation

This resource adds functional transcriptome evidence for PDGFB and pathway-level
LPS evidence for PRKCQ. It does not provide allele-specific validation because
the perturbations target genes rather than patient variants. IL23R must not be
interpreted as negative.

## Outputs

- `results/tables/perturbseq_gse327644_pseudobulk.csv.gz`
- `results/tables/perturbseq_gse327644_cell_assignment.csv`
- `results/tables/perturbseq_gse327644_target_knockdown.csv`
- `results/tables/perturbseq_gse327644_differential_summary.csv`
- `results/tables/perturbseq_gse327644_guide_effects.csv`
- `results/tables/perturbseq_gse327644_pathway_effects.csv`
- `results/tables/perturbseq_gse327644_differential_effects.csv.gz`
- `results/tables/perturbseq_umi_threshold_sensitivity.csv`
- `results/tables/perturbseq_leave_one_guide_out.csv`
- `results/tables/perturbseq_leave_one_donor_out.csv`
- `results/tables/perturbseq_mixed_model_effects.csv`
- `results/tables/perturbseq_robustness_summary.json`
- `results/figures/figure9_perturbseq_validation.png`
- `results/figures/figure10_evidence_robustness.png`
