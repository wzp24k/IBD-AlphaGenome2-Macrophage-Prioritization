# Top-three variant primary macrophage MPRA audit

Audit date: 2026-11-13

## Public-data search

PubMed exact-variant searches returned zero records for rs2143178, rs7404003
and rs971891 on 2026-11-13. GEO searches combining each rsID with MPRA,
massively parallel reporter or reporter assay also returned zero records.

No existing exact-variant MPRA result was promoted to evidence for these three
variants. This is a public-coverage result, not a biological negative.

## Experiment selection

The first follow-up should be allele-specific lentiMPRA rather than CRISPRi.
MPRA directly tests whether the reference and alternate alleles drive different
reporter activity. CRISPRi at PRKCB, CHSY3 and PDGFB tests gene function but
cannot establish that a specific patient variant causes the observed effect.

## Candidate set

The selected experimental panel is PDGFB rs2143178, PRKCB rs7404003 and CHSY3
rs971891. These are not the literal first three global priority ranks. They are
the three strongest candidates for an interpretable variant-to-gene experiment:

- PDGFB: exact cCRE/topic annotation and effective Perturb-seq knockdown.
- PRKCB: exact colon eQTL and gene-level perturbation support.
- CHSY3: gene-level target engagement and inflammatory pathway evidence.

The literal global priority ranks 1 and 2 are TNFSF15 rs10117785 and CXCR2
rs4672872. They should not be conflated with the selected experimental panel.

## Generated design

- Genome build: hg38.
- Three 251 bp overlapping sequence contexts per variant.
- Both alleles for every sequence context.
- ETS2 positive control, IL2RA negative control and two scrambled sequences.
- Eight proposed barcodes per construct.
- Alternative/reference contrast rather than unverified risk/protective labels.

All reference alleles were checked against the UCSC hg38 sequence. The design is
ready for internal scientific and vector-compatibility review, but no synthetic
library, transduction, sequencing or biological result exists yet.
