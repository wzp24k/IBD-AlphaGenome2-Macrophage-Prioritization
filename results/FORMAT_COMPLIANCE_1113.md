# Manuscript format compliance

Date: 2026-11-13

## Document format

- US Letter page size.
- One-inch margins on all four sides.
- Times New Roman, 12 pt body text.
- 1.5 line spacing.
- 18 pt bold title and 14 pt first-level headings.
- No header or footer.
- Embedded figures with captions below.
- Three-line tables with notes below.

## Text and statistics

- Structured abstract: 249 words.
- Abstract labels: Background, Methods, Results and Conclusion.
- Keywords: five title-case terms.
- No prose sentence exceeds 30 words.
- P values use three decimal places when P >= 0.001.
- P values below 0.001 are reported as `P < 0.001`.
- Strong causal wording was reduced where the evidence is predictive,
  observational or pathway-level.

## Figures and tables

- Figure 1 is a four-stage technical roadmap.
- Figures 2-12 use a unified six-color palette.
- Multi-panel figures include A-F panel labels and captions below.
- Figure types include bubble, heatmap, bar, interval, ECDF, strip, line and
  concordance displays.
- Tables 1-11 are sequential, three-line tables with interpretation notes.
- The main figure-generation pipeline remains Python. R 4.5.3 is used only for
  SuSiE fine-mapping and colocalization; Figure 12 is rendered in Python from
  the exported R results.

## References

- 42 references, all verified as `journal-article` records in Crossref.
- Every reference has a valid DOI.
- Citation numbers follow first appearance from 1 to 42.
- Every reference is cited in the manuscript.
- No preprint, dataset record, website or non-journal source is included.

## Declarations

The manuscript includes Data availability, Code availability, Ethics approval
and consent to participate, Patient and public involvement, Authors'
contributions, Funding, Generative AI declaration, Consent for publication,
Clinical trial number and Competing interests.

The public code repository is:
https://github.com/wzp24k/IBD-AlphaGenome2-Macrophage-Prioritization

It excludes the manuscript, raw data, large single-cell objects, API credentials
and unverified files.

## Final checks

The three review passes covered content and logic, numbers and figure citations,
reference numbering and DOI integrity, document formatting, and final package
integrity. All automated checks in `QC_REPORT.json` pass.
