# Peer Reviewer 5 full audit

## Manuscript information

- File: `1114manuscript.docx`
- Review date: 2026-11-14
- Review focus: final semantic consistency, selection effects, numeric traceability,
  display-item ordering and declaration consistency
- Calibration status: `NOT_CALIBRATED`

## Recommendation

**Minor revision for a computational genetics or translational-prioritization
journal.** No unresolved numerical contradiction remains. A causal experimental
journal would still require major revision because no new allele-specific assay
is reported.

## Corrections made

### W1: Top-28 correlation was overinterpreted

**Finding**: The top 28 loci were selected using the same evidence axes later
correlated across healthy and active macrophages. A high correlation in this
selected subset does not prove that the global prioritization was unaffected by
absolute-effect selection.

**Correction**: The manuscript now reports the stronger correlation among the
preselected top 28 without using it to claim absence of selection effects.

**Severity**: Minor. **Evidence Anchor**: text:
"Among the preselected top 28 loci, correlation was stronger."
**Confidence**: 5.

### W2: CEDAR2 locus-level and exact-variant evidence were conflated

**Finding**: The prior sentence could imply that the exact PRKCB variant was
replicated in monocytes. CEDAR2 showed positive lead effects in monocytes, but
the exact variant was significant only in granulocytes.

**Correction**: The text now distinguishes locus-level direction from exact-variant
evidence and states that the exact variant was not supported in a colon
macrophage node.

**Severity**: Minor. **Evidence Anchor**: table:
`prkcb_cedar2_independent_eqtl.csv` exact-variant and lead-effect rows.
**Confidence**: 5.

### W3: "Model novelty" was the wrong evidence concept

**Finding**: Model novelty is not an evidence layer. The relevant dimension is
model prediction and its calibration.

**Correction**: The Abstract and Discussion now use "model prediction."

**Severity**: Minor. **Evidence Anchor**: text:
"distinguish model prediction, cell-context specificity and causal evidence."
**Confidence**: 5.

### W4: Neutral CRISPR wording was defensively framed

**Finding**: The phrase "not failed assays or evidence against target-mRNA
effects" mixed interpretation of assay validity with the measured endpoint.

**Correction**: The manuscript now states directly that the results are neutral
for cytokine selection and do not measure target-mRNA effects.

**Severity**: Minor. **Evidence Anchor**: text:
"These results are neutral for cytokine selection."
**Confidence**: 5.

## Final integrity checks

- 42 references, each cited once and ordered by first appearance from 1 to 42.
- Discussion: 15 distinct references; earlier sections: 27.
- Citation groups: at most two references, sorted and unique.
- Word citation hyperlinks: all 42 occurrences point to the correct bookmarks.
- Tables 1-11 and Figures 1-12: body references, captions and first-use order
  all match.
- Core cohort, pathway, perturbation and PRKCB mechanism values match source
  tables.
- Clinical trial number NCT00572585 matches the Results and Discussion.
- No prose sentence exceeds 30 words.

## Remaining scientific limitations

- No new allele-specific MPRA or primary macrophage CRISPRi result is included.
- PRKCB direction differs between GTEx colon and myeloid datasets.
- Five Tier 1 targets remain experimentally uncovered.
- Public data vary in ancestry, treatment exposure and disease duration.
