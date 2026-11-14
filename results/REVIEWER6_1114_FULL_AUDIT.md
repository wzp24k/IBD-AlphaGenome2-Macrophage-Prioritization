# Peer Reviewer 6 full audit

## Manuscript information

- File: `1114manuscript.docx`
- Review date: 2026-11-14
- Review focus: final numerical provenance, bootstrap denominators, cross-section
  consistency and production-ready manuscript integrity
- Calibration status: `NOT_CALIBRATED`

## Recommendation

**Minor revision for a computational genetics or translational-prioritization
journal.** No unresolved numerical contradiction remains in the corrected
version. A venue requiring experimental causal closure would still require
major revision.

## Corrections made

### W1: Benchmark bootstrap count was wrong

**Finding**: Methods reported 500 locus-block bootstrap replicates. Main benchmark
AUROC and AUPRC metrics contain 1,000 valid replicates. The paired TNF CRISPR
differences contain 999 valid replicates.

**Correction**: Methods now report 1,000 replicates for AUROC and AUPRC and
explicitly state that 999 valid paired replicates were retained for the TNF
CRISPR endpoint.

**Severity**: Major. **Evidence Anchor**: table:
`benchmark_model_metrics.csv` and `benchmark_model_deltas.csv`
valid_bootstrap_replicates columns.
**Confidence**: 5.

### W2: Top-28 correlation could be overinterpreted

**Finding**: The top 28 loci were selected using related evidence axes. Their
high healthy-versus-active correlation cannot establish that selection effects
were absent globally.

**Correction**: The text now reports the preselected subset result without an
absence-of-selection claim.

**Severity**: Minor. **Evidence Anchor**: text:
"Among the preselected top 28 loci, correlation was stronger."
**Confidence**: 5.

### W3: CEDAR2 exact-variant and locus-level evidence were mixed

**Finding**: CEDAR2 showed positive PRKCB lead effects in monocytes, but the
exact rs7404003 association was significant only in granulocytes.

**Correction**: The Results now distinguish locus-level direction from
exact-variant evidence.

**Severity**: Minor. **Evidence Anchor**: table:
`prkcb_cedar2_independent_eqtl.csv`.
**Confidence**: 5.

### W4: Public evidence taxonomy was inconsistent

**Finding**: GTEx eQTL was grouped with experimental perturbation evidence.

**Correction**: The manuscript now uses "public molecular and perturbation
evidence." The abstract separately identifies exact-variant eQTL and gene-level
perturbation coverage.

**Severity**: Minor. **Evidence Anchor**: table:
`Table 8` public molecular and perturbation evidence.
**Confidence**: 5.

### W5: Donor-derived percentile omitted the absolute-value definition

**Finding**: The ranking uses the percentile of the absolute donor-derived
macrophage effect, not the signed effect percentile. The prior wording did not
distinguish these quantities.

**Correction**: Methods, Results and Table 5 now refer to the donor-derived
macrophage absolute-effect percentile.

**Severity**: Minor. **Evidence Anchor**: table:
`Table 5` donor-derived macrophage absolute-effect percentile.
**Confidence**: 5.

## Final numerical and display audit

- Benchmark replicate provenance: 1,000 valid replicates for primary benchmark
  endpoints and 999 for TNF CRISPR.
- Core cohort counts: 395 variants and 281 credible sets.
- Tier counts: 10, 57, 138 and 76, totaling 281.
- Study-level counts total 1,490.
- External analysis: 510 genes, 73 supported genes and 19 genes supported by
  both cohort clusters.
- Four pathway rows match source confidence intervals and FDR.
- PRKCB MR, SMR, SuSiE, HEIDI, leave-one-out and virtual-knockdown values match
  source tables.
- Tables 1-11: body references, captions and first-use order all match.
- Figures 1-12: body references, captions and first-use order all match.
- References: 42, each cited once and first used in order 1-42.
- Discussion uses 15 references; earlier sections use 27.
- No prose sentence exceeds 30 words.

## Remaining scientific limitations

- No new allele-specific MPRA or primary macrophage CRISPRi is included.
- PRKCB direction differs between GTEx colon and myeloid datasets.
- Five Tier 1 targets remain experimentally uncovered.
- Public resources vary in ancestry, treatment exposure and disease duration.
