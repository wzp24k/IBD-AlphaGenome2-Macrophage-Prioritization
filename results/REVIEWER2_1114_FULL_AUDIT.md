# Peer Reviewer 2 full audit

## Manuscript information

- File: `1114manuscript.docx`
- Review date: 2026-11-14
- Review scope: full manuscript, references, numbers, logic, figures, tables and captions
- Calibration status: `NOT_CALIBRATED`

## Reviewer role and focus

Peer Reviewer 2: human and IBD regulatory genomics. The review focused on
domain accuracy, evidence hierarchy, direction consistency, target interpretation
and whether the conclusions match the presented data.

## Recommendation

**Minor revision after the corrections recorded below.**

The corrected manuscript has no unresolved internal numerical contradiction.
If the target journal expects a new experimental causal mechanism, the study
would instead require major revision because no new allele-specific experiment
is reported. The paper is internally consistent when framed as computational
prioritization and evidence triangulation.

## Strengths

### S1: The evidence ladder is explicit

The manuscript separates fine-mapping, calibrated prediction, macrophage
context, disease-state stability, independent genetics, public perturbation
and LD-aware mechanism analysis.

**Evidence Anchor**: text: "The analysis follows an evidence ladder."

### S2: Candidate interpretation distinguishes evidence classes

PDGFB is presented as the leading macrophage-regulatory candidate. PRKCB is
presented as the mechanistic-triangulation candidate. Gene-level perturbation,
chromatin annotation and allele-specific evidence are not treated as equivalent.

**Evidence Anchor**: text: "Among the Top 10, PDGFB is the leading macrophage-regulatory archetype."

### S3: Reproducibility and reporting controls are extensive

The package includes scripts, source tables, figure source data, 42 linked
references, nine declaration categories and audited figure/table captions.

**Evidence Anchor**: dataset: `IBD_AlphaGenome2_submission_package.zip`.

## Corrections made

### W1: PRKCB direction was described inconsistently

**Finding**: GTEx sigmoid-colon PRKCB NES is -0.171. The monocyte eQTL used for
MR is +1.138, and the queried myeloid datasets are also positive. The original
text stated that GTEx agreed with AlphaGenome and later described the risk
allele as increasing expression without identifying the tissue-context change.

**Correction**: The manuscript now states that the colon and myeloid directions
are discordant. Inhibitor relevance is restricted to the myeloid evidence.
Table 8 is labeled context-dependent. Table 11 includes a direction-boundary row.

**Severity**: Major. **Evidence Anchor**: table:
`Table 8` PRKCB GTEx NES = -0.17. **Confidence**: 5, core expertise:
IBD regulatory genomics.

### W2: Top-28 active-disease concordance wording was imprecise

**Finding**: The source table reports 21 loci supported by the composite
active-disease definition and 17 strictly replicated across both independent
cohorts. The prior sentence could be read as 21 loci concordant in at least one
independent cohort, which is not the definition used.

**Correction**: The text now states that 21 were concordant with the original
active pool plus at least one independent cohort. Seventeen were concordant in
both Crohn and Garrido cohorts.

**Severity**: Major. **Evidence Anchor**: table:
`ranked_lead_variants_active_ibd.csv` top-28 support counts. **Confidence**: 5,
core expertise: IBD single-cell analysis.

### W3: Claim calibration omitted the PRKCB direction boundary

**Finding**: The seven-row calibration table did not capture the newly identified
colon-versus-myeloid direction difference.

**Correction**: The table now has eight boundaries and explicitly permits only a
context-dependent PRKCB interpretation.

**Severity**: Minor. **Evidence Anchor**: table:
`Table 11` PRKCB eQTL direction row. **Confidence**: 5.

### W4: MR-Egger and heterogeneity wording was too categorical

**Finding**: "Unstable" was less precise than the observed evidence. Egger was
imprecise with three instruments and discordant with IVW. Local HEIDI
heterogeneity did not prove an artifact-free local model.

**Correction**: The manuscript now states that Egger was imprecise and discordant.
The local-window result is described as persistent heterogeneity rather than proof
that distant variants caused the full-region result.

**Severity**: Minor. **Evidence Anchor**: text:
"MR-Egger was imprecise with three instruments and discordant with IVW."
**Confidence**: 5.

## Full checks passed

- 395 primary variants and 281 independent credible sets.
- Study-level variant counts sum to 1,490.
- Tier counts sum to 281: 10, 57, 138 and 76.
- External analysis includes 510 genes; 73 are externally supported and 19 are
  supported by both cohort clusters.
- Four external pathway rows match the source confidence intervals and FDR.
- PRKCB IVW, SMR, SuSiE H4, leave-one-out, HEIDI and virtual-knockdown values
  match source tables.
- PDGFB, PRKCQ and IL23R Perturb-seq values match source tables.
- Tables 1-11 all have in-text citations and matching captions.
- Figures 1-12 all have in-text citations and matching captions.
- References: 42 entries, each cited once, first use ordered 1-42.
- Discussion uses 15 references; earlier sections use 27.
- No prose sentence exceeds 30 words.

## Remaining scientific limitations

- No new allele-specific MPRA or primary macrophage CRISPRi experiment is included.
- Five Tier 1 targets lack supportive experimental evidence.
- pQTL, sQTL and caQTL layers are not covered for rs7404003.
- PRKCB MR uses three approximately independent instruments.
- Public perturbation resources differ in cell type, assay and endpoint.

## Audit artifacts

- `results/1114_reference_audit.json`
- `results/1114_figure_table_audit.json`
- `results/1114_numeric_consistency_audit.json`
- `results/tables/reference_citation_audit_summary.json`
- `results/tables/reference_doi_verification_summary.json`
