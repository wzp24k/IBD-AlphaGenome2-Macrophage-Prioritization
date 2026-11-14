# Peer Reviewer 3 full audit

## Manuscript information

- File: `1114manuscript.docx`
- Review date: 2026-11-14
- Review focus: cross-disciplinary translation, practical implementation,
  overinterpretation risk and broader validity
- Calibration status: `NOT_CALIBRATED`

## Reviewer identity

Cross-disciplinary reviewer in translational genomics and health-data
implementation. The review asked whether the framework would help research
teams choose experiments and whether clinical readers could reasonably
interpret the candidate labels.

## Recommendation

**Minor revision for a computational genomics or translational-prioritization
journal.** A journal that expects a completed wet-laboratory causal mechanism
would still require major revision.

## Strengths

### S1: The output is organized around decisions

The evidence ladder separates target prioritization, donor-state stability,
independent validation, public perturbation and mechanism analysis. This is more
useful than a single opaque score for experimental triage.

**Evidence Anchor**: text: "The analysis follows an evidence ladder."

### S2: The manuscript distinguishes experimental scales

Chromatin annotation, reporter assays, gene-level perturbation and virtual
prediction are not merged into one validation label.

**Evidence Anchor**: text: "Gene-level experiments do not establish that a patient variant causes the observed effect."

### S3: The major clinical inference is bounded by cell context

The revised text states that PRKCB direction differs between colon tissue and
myeloid cells. The inhibitor rationale is restricted to myeloid evidence.

**Evidence Anchor**: table: `Table 11` PRKCB eQTL direction row.

## Corrections made

### W1: Benchmark language implied stronger evidence than the data support

**Finding**: The phrase "excluding dependence on post hoc score construction"
was stronger than the near-identical benchmark version could establish. The
Top-10 eQTL comparison also used a small number of pairs.

**Correction**: The text now states that the conclusions were not driven by the
post hoc construction and that the eQTL result "supports enrichment." The word
"strong" was removed from both claims.

**Severity**: Major. **Evidence Anchor**: text:
"indicating that conclusions were not driven by the post hoc score construction."
**Confidence**: 4, cross-disciplinary benchmarking perspective.

### W2: Practical implementation barriers were underreported

**Finding**: The manuscript discussed biological limitations but not donor
consent, primary macrophage availability, stimulation standardization, ancestry,
medication exposure or transfer across health systems.

**Correction**: A practical implementation paragraph was added. It identifies
these constraints without changing the performance claims.

**Severity**: Minor. **Evidence Anchor**: absence:
"Discussion limitations" - expected implementation barriers; checked Discussion
and Limitations. **Confidence**: 4.

### W3: Candidate labels lacked scope

**Finding**: "Leading candidate" can be read as a field-wide claim rather than a
statement about the prioritized Top 10.

**Correction**: PDGFB and PRKCB archetype statements and the conclusion now say
"among the Top 10."

**Severity**: Minor. **Evidence Anchor**: text:
"Among the Top 10, PDGFB is the leading macrophage-regulatory archetype."
**Confidence**: 4.

### W4: Virtual knockdown was described as defining a mechanism

**Finding**: A computational perturbation should generate a hypothesis, not
define an established mechanism.

**Correction**: The text now says that the predictions support a
context-dependent JAK/STAT hypothesis for functional follow-up.

**Severity**: Minor. **Evidence Anchor**: text:
"These predictions support a context-dependent JAK/STAT hypothesis."
**Confidence**: 5.

## Cross-disciplinary assessment

### Practical impact

The framework can help research teams decide which candidates receive
allele-specific reporter assays and which receive gene-level perturbation.
It does not replace assay design, donor selection or clinical evaluation.

### Stakeholders

The principal users are IBD genetics and macrophage-biology teams. Clinical
readers can interpret the drug annotations as translational leads, while the
manuscript keeps prospective efficacy and target-selective safety separate.

### Equity and generalizability

The added limitations note donor-ancestry, medication-exposure and
disease-duration variation. These are material because public single-cell
resources do not represent all IBD populations or treatment settings.

### Unintended use

The evidence labels reduce the risk that a reader will treat a model score as a
diagnostic or treatment-selection result. The manuscript keeps the output tied
to experimental triage.

## Full checks passed

- 42 references, each cited once and ordered by first use from 1 to 42.
- Discussion uses 15 references; earlier sections use 27.
- Tables 1-11 and Figures 1-12 all have body references and matching captions.
- Core cohort, pathway, Perturb-seq and PRKCB mechanism statistics match sources.
- No prose sentence exceeds 30 words.
- Final Word/format QC and numeric consistency audit both pass.

## Remaining limitations

- No new allele-specific MPRA or primary macrophage CRISPRi is included.
- PRKCB has context-dependent eQTL direction.
- Five Tier 1 targets remain experimentally uncovered.
- Public resources underrepresent some ancestries and treatment states.
