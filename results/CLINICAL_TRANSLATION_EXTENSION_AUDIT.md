# Clinical translation extension

Generated: 19 November 2026

## 1. CD/UC and cross-ancestry validation

- Top 10 targets supported by Crohn credible-set L2G at score >= 0.5: 5.
- Top 10 targets supported by UC credible-set L2G at score >= 0.5: 6.
- Targets supported in both diseases: 4.
- PRKCB cis-MR was positive in all four outcome cohorts:
  - Combined IBD: beta = 0.0511, P = 1.07e-17.
  - Crohn disease: beta = 0.0351, P = 1.91e-5.
  - Ulcerative colitis: beta = 0.0661, P = 7.83e-19.
  - East Asian IBD: beta = 0.0484, P = 4.09e-4.

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

- PRKCB: 7.5/10, highest readiness class.
- PDGFB: 5.0/10, high.
- MAP3K8: 4.0/10, moderate.
- PRKCB combined cross-disease genetics, cross-ancestry MR, exact eQTL, eQTL colocalization, public perturbation, positive pQTL direction, treatment-response programs and direction-compatible clinical-stage annotation.

## Interpretation boundary

The extension is an evidence-integration analysis. pQTL MR and PheW-MR are exploratory and rely on available public instruments and phenotypes. The treatment-response analysis is observational and should be interpreted as program-level patient-enrichment evidence.
