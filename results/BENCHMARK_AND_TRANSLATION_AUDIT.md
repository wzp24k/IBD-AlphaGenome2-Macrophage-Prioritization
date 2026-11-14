# External benchmark and translational annotation audit

## Exact-variant eQTL benchmark

- Candidate-target pairs queried: 384.
- Core macrophage or monocyte datasets: 4.
- Exact-variant eQTL in at least one dataset: 140/384.
- Replicated nominal eQTL in at least two datasets: 84/384.
- Top 10 matched-target pairs queried across 8 datasets: 7/10 nominal.
- Top 10 mismatched-target negative controls: 0/10 nominal.
- Fisher exact test: P = 0.003.
- Bonferroni threshold 0.05/8: 95 positive candidate-target pairs.
- Top 10 under the strict threshold: still 7/10 matched and 0/10 controls.

## Formal colocalization

- Top 10 eQTL-positive pairs evaluated by coloc.abf: 7.
- Pairs with primary H4 >= 0.8: 5/7.
- Supported: DENND1B, MAP3K8, PDGFB, PRKCB and TNFSF15.
- Not supported: CXCR2 (primary H4 = 0.746) and CDYL2
  (primary H4 = 0.231).

## Model comparison

External genetic support in both cohorts:

- Combined AUROC 0.673 (95% CI 0.579-0.775); AUPRC 0.057.
- PIP-only AUROC 0.652; AUPRC 0.051.
- L2G-only AUROC 0.434; AUPRC 0.042.
- Atlas-only AUROC 0.584; AUPRC 0.044.
- Live-only AUROC 0.537; AUPRC 0.041.

The combined score improved over L2G-only by Delta AUROC 0.239
(95% CI 0.094-0.396). The difference from PIP-only was 0.022
(-0.117 to 0.144), so universal superiority over simpler PIP ranking
was not established.

Replicated exact-variant eQTL:

- Combined AUROC 0.559 (95% CI 0.453-0.654); AUPRC 0.277.
- L2G-only AUROC 0.617; PIP-only AUROC 0.536.
- Superiority over PIP-only or L2G-only was not interval-confirmed.

Macrophage TNF CRISPR:

- Evaluable matched pairs: 77, including only 7 positive pairs.
- Combined AUROC 0.553; L2G-only 0.606.
- The sparse endpoint did not establish combined-score superiority.

## Calibration and transport

- External support, both cohorts: Brier score 0.035 and expected
  calibration error 0.0018.
- Replicated eQTL: Brier score 0.171 and expected calibration error 0.0104.
- UKB-to-MVP transport AUROC: 0.633.
- MVP-to-UKB transport AUROC: 0.614.
- Held-out credible-set grouped AUROC: 0.684.
- Held-out target-gene grouped AUROC: 0.715.

UKB and MVP transport is a cohort generalisation proxy. It is not a
formal ancestry-specific validation.

## Clinical translation annotation

- Top 10 targets with tractable modalities: 9/10.
- Top 10 targets with clinical candidates: 5/10.
- Clinical-stage IBD/UC/Crohn-indication signals: CXCR2, PRKCB and
  MAP3K8, each at phase 2.
- PDGFB and PLAU had non-IBD clinical candidates.
- MAP3K8 had the cleanest combination of variant-to-gene evidence,
  tractability and a target-specific inhibitor.
- PRKCB clinical pharmacology was dominated by pan-PKC inhibitors.
- CXCR2 did not reach the primary colocalization threshold.
- Clinical treatment recommendations supported: 0/10.

Clinical-stage IBD-indication signals are drug-repurposing hypotheses.
Treatment recommendations require established therapeutic direction,
target selectivity, prospective efficacy and safety evidence.

## Score-construction sensitivity

- Spearman correlation between the original priority score and the
  equal-weight benchmark score: 0.99998.
- External support in both cohorts: AUROC 0.67297 for the original score
  versus 0.67336 for the benchmark score.
- Replicated exact-variant eQTL: AUROC 0.55956 versus 0.55921.

The benchmark conclusions therefore do not depend on constructing a new
post hoc score.

## Evidence boundary

The exact-variant eQTL and colocalization results strengthen
variant-to-target evidence but remain observational. They do not establish
variant causality. The MPRA and Perturb-seq resources had too few
variant-level positive labels for stable AUROC estimation and are reported
as coverage rather than used to create sparse or misleading metrics.

Conclusions supported strongly:

- Matched Top 10 target assignments are enriched for exact-variant eQTL
  signal compared with mismatched targets.
- Five Top 10 candidate genes show convergent exact-variant eQTL and formal
  colocalization evidence across sensitivity settings.
- The original priority score reproduces the benchmark-score result.

Conclusions not supported:

- Universal superiority of the combined framework over every single-score
  baseline.
- Variant causality from eQTL or colocalization alone.
- IBD efficacy from clinical-stage drug annotations.

Clinical-stage hypotheses supported for further research:

- MAP3K8: strongest combined target-mapping, tractability and
  target-specific phase 2 inhibitor signal.
- PRKCB: strong target mapping with phase 2 pan-PKC pharmacology.
- CXCR2: phase 2 antagonist/inhibitor signal without formal
  colocalization support.
