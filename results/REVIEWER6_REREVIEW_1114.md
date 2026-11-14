# Peer Reviewer 6 re-review

## Manuscript

- File: `1114manuscript.docx`
- Review round: 2
- Review date: 2026-11-14
- Calibration status: `NOT_CALIBRATED`

## Decision

**Minor revision for a computational genetics or translational-prioritization
journal.** No unresolved numerical or citation inconsistency remains.

## New corrections verified

### 1. Benchmark bootstrap denominators

The Methods now distinguish:

- 1,000 valid bootstrap replicates for benchmark AUROC and AUPRC metrics.
- 999 valid paired resampling replicates for the TNF CRISPR comparison.

This matches both `benchmark_model_metrics.csv` and
`benchmark_model_deltas.csv`.

### 2. Donor-derived percentile definition

The ranking and narrative now specify donor-derived macrophage
absolute-effect percentile. This matches the computation in
`10_rank_candidates.py`.

### 3. Selection-bias wording

The Top 28 correlation is reported as a descriptive result within a preselected
subset. The text no longer uses it to claim absence of selection effects.

### 4. CEDAR2 precision

The text distinguishes positive locus-level effects in monocytes from the
exact-variant result available in granulocytes. It does not call the exact
variant replicated in colon macrophages.

## Final delivery verification

- Word paragraphs: 265.
- Tables: 11.
- Embedded figures: 12.
- Excel sheets: 81.
- Submission ZIP entries: 305, with no duplicate archive paths.
- References: 42.
- Each reference is cited only once.
- First-use citation order: 1-42.
- Discussion references: 15.
- Pre-Discussion references: 27.
- Maximum references in one citation group: 2.
- All 42 citation occurrences link to reference bookmarks.
- Table first-use order: 1-11.
- Figure first-use order: 1-12.
- Sentences longer than 30 words: 0.
- Final QC: 158/158 passed.
- Numeric consistency audit: all checks passed.

## Remaining limitation

The paper remains a computational prioritization study. It does not contain a
new allele-specific MPRA or primary macrophage CRISPRi experiment.
