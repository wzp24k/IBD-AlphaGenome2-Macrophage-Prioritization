# Peer Reviewer 4 full audit

## Manuscript information

- File: `1114manuscript.docx`
- Review date: 2026-11-14
- Review focus: evidence taxonomy, declaration consistency, terminology drift,
  overstatement and whole-manuscript numerical coherence
- Calibration status: `NOT_CALIBRATED`

## Recommendation

**Minor revision for a computational genetics or translational-prioritization
journal.** A causal-mechanism venue would require major revision because the
study has no new allele-specific experimental assay.

## Corrections made

### W1: GTEx eQTL was misclassified as experimental evidence

**Finding**: The abstract, Results heading, Table 8, Figure 8 and Figure 11
described public evidence as "experimental." GTEx eQTL is an observational
molecular association, whereas the MPRA, CRISPR and Perturb-seq resources are
perturbational experiments.

**Correction**: The manuscript now uses "public molecular and perturbation
evidence" at section, table, figure and benchmark levels. The abstract separately
states "one exact-variant eQTL" and five gene-level perturbation targets.

**Severity**: Major. **Evidence Anchor**: table:
`Table 8` GTEx eQTL listed under the former experimental heading.
**Confidence**: 5.

### W2: Clinical trial declaration contradicted the Results

**Finding**: The manuscript discussed NCT00572585 but the declaration said
"Clinical trial number: Not applicable."

**Correction**: The declaration now states that the current work is a secondary
analysis and identifies NCT00572585 as the publicly reported trial cited.

**Severity**: Major. **Evidence Anchor**: text:
"The publicly reported trial cited in this study is NCT00572585."
**Confidence**: 5.

### W3: Planned lentiMPRA was described as prespecified

**Finding**: No evidence in the project artifacts establishes preregistration
of the lentiMPRA design.

**Correction**: "Prespecified" was replaced with "planned" in Methods and
Discussion.

**Severity**: Minor. **Evidence Anchor**: absence:
"Methods and Discussion" - expected a preregistration record before prespecified
language; checked project audit artifacts. **Confidence**: 5.

### W4: PRKCB overexpression was overstated as inflammatory validation

**Finding**: The THP-1 overexpression GSEA included growth and metabolic
programs. Only CHSY3 TNF/NF-kB survived growth-PC1 residualization.

**Correction**: The Discussion now says THP-1 overexpression changes
gene-expression programs, including inflammatory terms sensitive to growth
state.

**Severity**: Minor. **Evidence Anchor**: table:
`Table 10` growth-PC1 residualization row. **Confidence**: 5.

### W5: Discovery ranking language lacked scope

**Finding**: "The leading candidate" could be read as a field-wide claim.

**Correction**: The variant is now described as the "top-ranked candidate in
this discovery set." Archetype claims are restricted to "among the Top 10."

**Severity**: Minor. **Evidence Anchor**: text:
"The top-ranked candidate in this discovery set was rs2143178."
**Confidence**: 4.

## Whole-manuscript checks

- 395 primary variants and 281 independent credible sets.
- Tier counts 10, 57, 138 and 76 sum to 281.
- Study-level counts sum to 1,490.
- External analysis contains 510 genes, 73 supported genes and 19 genes
  supported by both cohort clusters.
- Four pathway rows match source effect sizes, confidence intervals and FDR.
- PRKCB MR, SMR, SuSiE, HEIDI and virtual-knockdown values match source tables.
- Tables 1-11 and Figures 1-12 all have body citations and matching captions.
- References: 42, each cited once and first used in order 1-42.
- Discussion uses 15 references; earlier sections use 27.
- No prose sentence exceeds 30 words.

## Remaining limitations

- No new allele-specific MPRA or primary macrophage CRISPRi result is reported.
- Five Tier 1 candidates remain experimentally uncovered.
- PRKCB direction differs between GTEx colon and myeloid datasets.
- Public resources vary in ancestry, medication exposure and disease duration.
