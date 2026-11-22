# Clinical translation extension

Generated: 2026-11-22.

## 1. CD/UC and cross-ancestry validation

- Crohn targets with L2G >= 0.5: 5.
- UC targets with L2G >= 0.5: 6.
- Targets supported in both diseases: 4.
- PRKCB cis-MR was positive in all four outcome cohorts.

## 2. pQTL, colocalization and PheW-MR

- Sun 2018 plasma pQTL variants: 2,701.
- Independent pQTL instruments: 1.
- Positive pQTL-MR outcomes: 4/4.
- Colocalization H4 >= 0.8: 2/4 outcomes.
- Maximum H4: 0.998.
- Legacy SomaScan and curated pQTL results remain sensitivity analyses.
- PheW-MR queried 3 traits; 2 were evaluable.

## 3. Public treatment-response programs

- Datasets: GSE14580, GSE111761.
- The PRKCB network program was lower in responders across both cohorts.
- The macrophage marker program showed the same direction.
- Direct PRKCB expression had a wider interval that included zero.

## 4. Clinical development readiness

- PRKCB: 7.5/10 and the highest-ranked target.
- IBD-indication clinical stage, not maximum stage across all indications, contributes to scoring.
- PRKCB ranked 1 across four prespecified weighting schemes.
- PRKCB ranked first in 92.5% and top two in 97.6% of 100,000 randomized weight draws.

## 5. Held-out incremental benchmark

- The prespecified PIP+context model improved AUROC by 0.072 for external support in any cohort (95% CI 0.003 to 0.148).
- The advantage was not stable for support in both cohorts or replicated exact-variant eQTL.

## Parallel robustness checks

- The robustness run used 44 logical CPUs and 12 permutation chunks.
- Treatment-response permutations: 20000 per feature.
- pQTL sensitivity retained a positive direction in 70% of tests.
- Maximum pQTL-colocalization H4 across prior combinations: 0.185.

## Interpretation boundary

The extension is an evidence-integration analysis. The primary Sun 2018 pQTL is a one-instrument protein layer. Legacy pQTL and PheW-MR are exploratory. Treatment-response results are observational program-level evidence.
