# Clinical translation extension

Generated: 19 November 2026

## 1. CD/UC and cross-ancestry validation

- Top 10 targets supported by Crohn credible-set L2G at score >= 0.5: 5.
- Top 10 targets supported by UC credible-set L2G at score >= 0.5: 6.
- Targets supported in both diseases: 4.
- PRKCB cis-MR was positive in all four outcome cohorts:
  - Combined IBD: beta = 0.0517, P = 3.25e-18.
  - Crohn disease: beta = 0.0351, P = 1.68e-5.
  - Ulcerative colitis: beta = 0.0672, P = 1.38e-19.
  - East Asian IBD: beta = 0.0518, P = 1.20e-4.

## 2. pQTL, colocalization and PheW-MR

- pQTL MR tests: 12.
- Positive PRKCB direction: 9/12 tests.
- Nominal positive pQTL-MR tests: 3/12.
- SomaScan PKCB2: combined IBD beta = 0.1365, P = 0.00699; UC beta = 0.2914, P = 3.0e-6.
- Curated deCODE pQTL: East Asian IBD beta = 0.8679, P = 0.00672.
- pQTL regional colocalization was evaluable but did not produce H4 >= 0.8; eQTL-based SuSiE-colocalization remains the stronger colocalization layer.
- Targeted PheW-MR queried three OpenGWAS traits, with two having a matching instrument. No evaluable phenotype reached FDR < 0.05.

## 3. Public treatment-response programs

- Discovery/public datasets: GSE14580 (24 UC samples before infliximab) and GSE111761 (6 CD samples during anti-TNF therapy).
- PRKCB network program was lower in responders:
  - GSE14580 SMD = -1.456.
  - GSE111761 SMD = -3.351.
  - Random-effects meta-analysis SMD = -2.042 (95% CI -3.758 to -0.325; P = 0.0197).
- Macrophage marker program meta-analysis SMD = -1.128 (95% CI -2.114 to -0.142; P = 0.0249).
- Direct PRKCB expression meta-analysis SMD = -1.716 (95% CI -3.772 to 0.341; P = 0.102).

## 4. Clinical development readiness

- PRKCB: 6.0/10, high readiness class and highest-ranked target.
- PDGFB: 5.0/10, high.
- MAP3K8: 4.0/10, moderate.
- PRKCB combined cross-disease genetics, cross-ancestry MR, exact eQTL, eQTL colocalization, positive exploratory pQTL direction, treatment-response programs and direction-compatible clinical-stage annotation. The pQTL contribution is weighted at 0.5 because available pQTL instruments are not genome-wide significant.
- Across four weighting schemes, PRKCB ranked first under evidence-only and clinical-weighted scoring and second under treatment-excluded scoring.
- In 100,000 randomized component-weight draws, PRKCB ranked first in 70.1% and top two in 90.4%.

## Parallel robustness checks

- Parallel workers used 44 available CPU threads with 16 permutation chunks per feature.
- PRKCB cross-ancestry leave-one-instrument-out estimates remained positive; beta range 0.0349-0.0687.
- Treatment-response label permutations used 20,000 permutations per feature.
- PRKCB network program permutation P = 0.0002; macrophage marker program P = 0.0036.
- pQTL MR sensitivity included 20 IVW/weighted-median/MR-Egger/leave-one-out estimates; 70% retained a positive direction.
- pQTL colocalization used 72 prior combinations; maximum H4 was 0.185, so pQTL colocalization remains exploratory.

## Interpretation boundary

The extension is an evidence-integration analysis. pQTL MR and PheW-MR are exploratory and rely on available public instruments and phenotypes. The treatment-response analysis is observational and should be interpreted as program-level patient-enrichment evidence.
