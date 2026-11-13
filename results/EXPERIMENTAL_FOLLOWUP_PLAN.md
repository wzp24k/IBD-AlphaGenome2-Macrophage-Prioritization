# Follow-up experiments for the five uncovered Tier 1 candidates

Date: 2026-11-12

## Candidate panel

| Target | Variant | Reference | Alternate |
|---|---|---|---|
| ERGIC1 | rs564109 | A | T |
| GCKR | rs780094 | T | C |
| ABCB5 | rs76380568 | T | C |
| IRF2BP2 | rs10910476 | C | T |
| CAMK2A | rs1549913 | G | A |

## 1. Allele-specific MPRA in primary macrophages

### Design

- Test both alleles for all five variants.
- Use 200-250 bp sequences centered on each variant.
- Include at least three sequence contexts per variant by shifting the tested
  interval by -75, 0 and +75 bp.
- Assign at least five random barcodes per allele-context construct.
- Include rs2836882/ETS2 as the published macrophage MPRA positive control.
- Include inert scramble constructs and empty reporter controls.

### Cells and conditions

- Primary monocyte-derived macrophages from at least three independent donors.
- Resting macrophages.
- LPS stimulation, 100 ng/mL for 6 h.
- IFN-gamma stimulation, 20 ng/mL for 18 h.
- IL-23 stimulation, 50 ng/mL for 24 h.
- A combined LPS + IFN-gamma condition is recommended as a fourth state.

### Readout and statistics

- Sequence DNA and reporter RNA barcodes.
- Normalize RNA barcodes to input DNA barcodes.
- Primary effect: alternate-versus-reference log2 reporter activity.
- Use a donor-aware mixed model with context and barcode as random effects.
- Use QuASAR-MPRA or MPRAnalyze as the primary allele-specific model.
- Report per-donor and pooled effects, 95% confidence intervals and FDR.
- Replicate a nominal variant in an independent macrophage donor panel.

## 2. Promoter CRISPRi in primary macrophages

### Guides and cells

- Design four to six sgRNAs targeting -250 to +150 bp around each predicted
  target TSS: ERGIC1, GCKR, ABCB5, IRF2BP2 and CAMK2A.
- Include non-targeting sgRNAs and positive-control guides against TNFAIP3,
  MAP3K8 and TNF.
- Use dCas9-KRAB primary human monocyte-derived macrophages or iPSC-derived
  macrophages from at least three donors.

### Validation and readout

- Confirm target-gene knockdown by qPCR and, where possible, RNA-seq.
- Require at least 50% reduction in the target transcript.
- Measure target expression and secreted TNF, IL-23 and CXCL10.
- Measure conditions in resting, LPS, IFN-gamma and IL-23-treated macrophages.
- Use donor as a biological replicate and model guide effects nested in donor.
- Do not interpret promoter CRISPRi as allele-specific variant validation.

## 3. eQTL and colocalization

### Completed exact-variant screen

Eight public monocyte or macrophage eQTL datasets were queried:

- Alasoo 2018 macrophage naive and IFN-gamma.
- BLUEPRINT monocyte.
- Nedelec 2016 macrophage naive.
- Quach 2016 monocyte naive and LPS.
- Schmiedel 2018 monocyte naive and CD16-positive monocyte.

Results:

| Target | Minimum exact-variant P | Interpretation |
|---|---:|---|
| ERGIC1 | 0.078 | Public myeloid eQTL not covered |
| GCKR | 0.663 | Public myeloid eQTL not covered |
| ABCB5 | 0.751 | Public myeloid eQTL not covered |
| IRF2BP2 | 0.073 | Public myeloid eQTL not covered |
| CAMK2A | 0.00775 | Nominal IFN-gamma macrophage signal; direction discordant with AlphaGenome and FDR approximately 0.31 across the 40 tests |

CAMK2A should not be called a replicated eQTL. A formal colocalization is the
next step only after the GWAS and eQTL summary statistics are harmonized to the
same effect allele and LD reference.

### Follow-up

- Download full IBD and myeloid eQTL summary statistics for the CAMK2A locus.
- Run SuSiE-coloc or coloc-SuSiE with a prespecified prior and sensitivity
  analysis.
- Require a posterior probability of shared causal variant >= 0.80.
- If no signal is present, retain the result as "public experimental data not
  covered" rather than negative biological evidence.

## Submission interpretation

Until these experiments are completed, the five targets remain model-based
Tier 1 candidates. They should not be presented as experimentally validated or
as having failed validation.

## Completed perturbation extension

GSE327644 Perturb-seq was processed for genes represented in the primary
macrophage library. PDGFB and PRKCQ gained transcriptome or LPS pathway-level
support. IL23R editing did not reduce the target transcript and is classified as
inconclusive. The remaining candidates still require allele-specific MPRA or
promoter CRISPRi.
