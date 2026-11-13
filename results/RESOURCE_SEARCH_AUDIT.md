# External resource search audit

Audit date: 2026-11-13

## MPRAbase and MPRAvar

Public MPRA portal endpoints were queried on 2026-11-13. MPRAbase returned an
SSL/502 gateway failure and MPRAvar did not complete an SSL handshake. A second
search used PubMed and GEO exact-variant queries for rs564109, rs780094,
rs76380568, rs10910476 and rs1549913.

- rs564109, rs76380568, rs10910476 and rs1549913 returned no MPRA, STARR-seq or
  CRISPR exact-variant records.
- rs780094 returned two publications but no retrievable exact-variant assay
  table that could be integrated under the same REF/ALT and cell-context rules.
- No result was promoted to experimental evidence without a retrievable assay
  record.

## Gut scATAC

The 2026 Crohn intestinal scATAC atlas was resolved to GSE273194 and Zenodo
record 21508422. Its exact-variant table was queried for all 10 Tier 1 variants.

- rs2143178 was present in a differentially accessible DC-macrophage cCRE
  (chr22:39264648-39265148); the published ABC maximum peak linked the cCRE to
  PDGFB.
- rs780094 was present but did not overlap a cCRE.
- The other eight variants were absent.

The scATAC result is retained as exact-variant chromatin annotation and is not
described as an allele-specific perturbation.

## eQTL

GTEx v10 colon and eight eQTL Catalogue monocyte/macrophage datasets were
queried by exact variant and nominated target. Results are in
`experimental_eqtl_catalogue_coverage.csv` and
`experimental_gtex_prkcb_eqtl.csv`.
