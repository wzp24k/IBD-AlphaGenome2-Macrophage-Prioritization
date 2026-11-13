# External experimental evidence audit

Audit date: 2026-11-12

## Scope

Public experimental data were integrated only when the assay measured an exact
Tier 1 variant, a nominated target gene, or a specified pathway node.
Association statistics, L2G predictions and single-cell-conditioned model
outputs were not relabeled as experimental evidence.

## Variant-level evidence

GTEx v10 single-tissue eQTL queries were performed for all 10 Tier 1 variants
and their nominated targets in sigmoid and transverse colon.

- rs7404003 (chr16:23853326:C>A) was associated with PRKCB expression in
  sigmoid colon: P = 3.02 x 10^-5, normalized effect size = -0.171.
- No other Tier 1 variant-target pair passed the same exact-variant query.
- The PRKCB eQTL direction is concordant with the negative AlphaGenome effect.

Eight eQTL Catalogue monocyte or macrophage datasets were also queried for the
five remaining variants. ERGIC1, GCKR, ABCB5 and IRF2BP2 had no nominal
variant-target eQTL. CAMK2A rs1549913 had a nominal IFN-gamma macrophage
association (P = 0.00775, beta = -0.370), but the direction was discordant and
the Bonferroni value across 40 tests was approximately 0.31. It is retained as
a nominal discordant signal, not as validated colocalization.

GSE229472 was used as an external positive-control MPRA, not as validation of a
Tier 1 variant. The resource contains eight RNA pools and three sequence
contexts for rs2836882. The mean risk-allele effect was +0.169 log2 units, with
positive effects in two of three contexts.

## Gene-level perturbation evidence

GSE327124 contains primary human myeloid CRISPR screens with four sgRNAs per
target and TNF or CD80 high/low selection.

- TNF knockout: -1.447 log2 high/low.
- TNFAIP3 knockout: +0.711.
- MAP3K8 knockout: -0.385; permutation P < 0.0001; FDR = 0.0044.
- PDGFB knockout: +0.005; P = 0.57.
- PRKCQ knockout: +0.018; P = 0.28.
- IL23R knockout: -0.018; P = 0.68 in TNF selection and -0.012; P = 0.82 in
  CD80 selection.

The positive controls establish that the selection direction and analysis are
not reversed. Repeating the comparison against non-targeting controls with
Mann-Whitney tests produced the same neutral interpretation. The neutral
findings are therefore retained and are not described as failed assays. The
screen measures cytokine-output selection, not allele-specific target-mRNA
regulation.

GSE175685 contains lentiviral overexpression of IBD-associated genes in THP-1
myeloid cells.

- PRKCB expression delta versus empty vector: 2.68; P = 1.2 x 10^-15.
- CHSY3 expression delta versus empty vector: 6.82; P = 1.2 x 10^-9.
- PRKCB overexpression changed 15 genes and CHSY3 changed 43 genes at
  FDR < 0.05 and absolute expression delta > 1.

These experiments establish target engagement in a myeloid line. They do not
test the patient variants and are not called variant-level validation.

## Coverage

| Tier 1 target | Exact variant experiment | Gene-level perturbation | Interpretation |
|---|---|---|---|
| PDGFB | No | TNF CRISPR | Tested-neutral cytokine endpoint |
| PRKCB | GTEx colon eQTL | THP-1 overexpression | Directional variant support plus target engagement |
| CHSY3 | No | THP-1 overexpression | Target engagement |
| PRKCQ | No | TNF CRISPR | Tested-neutral cytokine endpoint |
| IL23R | No | TNF and CD80 CRISPR | Tested-neutral in two cytokine endpoints |
| ERGIC1, GCKR, ABCB5, IRF2BP2, CAMK2A | No | No | Experimentally uncovered |

Coverage summary: 1/10 exact-variant, 5/10 gene-level perturbation, 3/10
tested-neutral cytokine endpoints, and 5/10 experimentally uncovered.

## Reproducibility

- `scripts/31_integrate_external_experimental_data.py`
- `scripts/32_external_experimental_figure.py`
- `results/tables/external_experimental_evidence.csv`
- `results/tables/experimental_gtex_prkcb_eqtl.csv`
- `results/tables/experimental_ets2_mpra_context_effects.csv`
- `results/tables/experimental_primary_macrophage_crispr_guides.csv`
- `results/tables/experimental_primary_macrophage_crispr_genes.csv`
- `results/tables/experimental_thp1_overexpression_engagement.csv`
- `results/tables/experimental_evidence_summary.json`

## Perturb-seq

GSE327644 was analyzed as an independent transcriptome-level perturbation
resource. PDGFB editing reduced its own transcript by 0.636 log2 units
(FDR = 0.0088) and changed 460 genes at FDR < 0.05. All four PRKCQ guides reduced
LPS IL-23 and IL-6/JAK/STAT3 pathway expression. IL23R editing was inconclusive
because its transcript was not reduced. Complete methods and results are in
`PERTURBSEQ_AUDIT.md`.

## Gut scATAC

GSE273194 / Zenodo 21508422 was queried for exact Tier 1 variants. rs2143178 and
rs780094 were present in the Crohn scATAC variant table. Only rs2143178 overlapped
a differentially accessible DC-macrophage cCRE
(chr22:39264648-39265148), and the published ABC maximum peak linked that cCRE
to PDGFB. rs780094 did not overlap a cCRE. The other eight Tier 1 variants were
absent from the resource. This is variant-level chromatin evidence, not an
allele-specific perturbation.
