# Sun 2018 PRKCB cis-pQTL

## Source

- Resource: eQTL Catalogue QTS000035 / QTD000584.
- Study: Sun 2018 plasma pQTL.
- Molecular trait: PRKCB.5475.10.3..1.
- Samples: 3,301.
- Assembly: GRCh38.
- Requested region: chr16:23,300,000-23,950,000.
- Observed region: chr16:23300993-23949568.

## Files

- `prkcb_sun2018_pqtl_cis_grch38.csv`: standardized PRKCB pQTL rows.
- `prkcb_sun2018_pqtl_ld_matched.csv`: pQTL rows with regional LD.
- `prkcb_sun2018_pqtl_ld_matrix.csv`: Pearson r matrix for matched variants.
- `prkcb_sun2018_pqtl_ld_variants.csv`: variant annotation for the LD matrix.
- `prkcb_sun2018_pqtl_mr_ibd.csv`: one-instrument cis-MR estimates.
- `prkcb_sun2018_pqtl_coloc.csv`: ABF colocalization estimates.

## Required Columns

`variant_id`, `rsid`, `chromosome`, `position_grch38`, `other_allele`, `effect_allele`, `beta`, `se`, `maf`, and `p`.

## LD

The LD matrix uses the existing 1000 Genomes high-coverage European reference. 865 pQTL variants matched the regional reference. Matrix values are Pearson r, not r-squared.

## Lead Signal

- Variant: chr16:23847471:T:C (rs28474709).
- Beta: 0.082620.
- SE: 0.011807.
- P: 3.141e-12.
