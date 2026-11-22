# Clinical translation extension

Generated: 2026-11-22.

## 1. CD/UC and cross-ancestry validation

- Crohn targets with L2G >= 0.5: 5.
- UC targets with L2G >= 0.5: 6.
- Targets supported in both diseases: 4.
- PRKCB cis-MR was positive in all four outcome cohorts.

## 2. pQTL, colocalization and PheW-MR

- pQTL MR tests: 12.
- Positive PRKCB direction: 9/12.
- Nominal positive pQTL-MR tests: 3/12.
- pQTL colocalization reached H4 >= 0.8 in 0 primary tests.
- Available pQTL instruments remain exploratory because they are not genome-wide significant.
- Exploratory pQTL evidence contributed at most 0.25 of 10 points.
- PheW-MR queried 3 traits; 2 were evaluable.

## 3. Public treatment-response programs

- Datasets: GSE14580, GSE111761.
- The PRKCB network program was lower in responders across both cohorts.
- The macrophage marker program showed the same direction.
- Direct PRKCB expression had a wider interval that included zero.

## 4. Clinical development readiness

- PRKCB: 5.75/10 and the highest-ranked target.
- IBD-indication clinical stage, not maximum stage across all indications, contributes to scoring.
- PRKCB remained at rank 1-4 across four prespecified weighting schemes.
- PRKCB ranked first in 77.1% and top two in 90.4% of 100,000 randomized weight draws.

## 5. Held-out incremental benchmark

- The prespecified PIP+context model improved AUROC by 0.072 for external support in any cohort (95% CI 0.003 to 0.148).
- The advantage was not stable for support in both cohorts or replicated exact-variant eQTL.

## Parallel robustness checks

- The robustness run used 44 logical CPUs and 44 permutation chunks.
- Treatment-response permutations: 20000 per feature.
- pQTL sensitivity retained a positive direction in 70% of tests.
- Maximum pQTL-colocalization H4 across prior combinations: 0.185.

## Interpretation boundary

The extension is an evidence-integration analysis. pQTL MR and PheW-MR are exploratory. Treatment-response results are observational program-level evidence.
