# Reviewer 1 methodology audit

Date: 2026-11-13

Scope: logic, numerical consistency, sample sizes, statistical interpretation,
figure and table citations, and claim calibration. Reference numbering and DOI
validation were excluded at the author's request.

## Corrections made

1. The external pathway bootstrap was described as resampling 293 loci. The
   upstream table contains 293 credible-set identifiers, but only 276 had
   usable target evidence and entered the bootstrap. The text and Table 3 now
   state 276 loci and explain the 293-to-276 reduction.

2. Atlas target-level scores were missing for 22 of 395 variants, not 13.
   Thirteen were indels and nine were SNVs without an interpretable
   target-level summary. The manuscript now reports all 22 and retains them in
   live-model analysis.

3. CAMK2A had been grouped with "uncovered" targets despite a nominal eQTL and
   formal colocalization analysis. The manuscript now states that five targets
   lacked supportive experimental evidence; CAMK2A is explicitly described as
   tested but discordant and non-colocalized.

4. TNFAIP3 had been called a positive control without verified self-transcript
   knockdown. The self-effect was +0.319 log2 with FDR 0.078. It is now labeled
   a downstream inflammatory-response control, not evidence of successful
   transcript knockdown.

5. Table 2, Table 7, Figure 1, Figure 2, Figure 5 and Figure 6 lacked explicit
   in-text first references. These were added so every displayed item is cited.

6. Table 1 priority scores were reported to one decimal place; they are now
   reported to two decimal places, matching the manuscript rule.

## Numerical checks passed

- 395 primary non-coding variants and 281 independent SuSiE credible sets.
- Study table sums to 1,490 eligible rows across six studies.
- 373 Atlas-scored variants, 386 live-model variants and 384 donor-conditioned
  variants.
- 10 Tier 1 and 57 Tier 2 loci; 138 Tier 3 and 76 Tier 4 loci sum to 281.
- 276 informative loci, 510 target genes, 73 externally supported genes and 19
  genes supported by both external cohort clusters.
- Pathway directions, confidence intervals and FDR values match the source
  table.
- PRKCQ attenuated classification and repeated-pool interpretation match the
  source tables.
- Perturb-seq self-effects and responsive-gene counts match the combined tables.
- Donor-aware modelling used 56 target/condition/pathway combinations: 54
  donor-fixed-effect fallbacks and two estimable donor-random-intercept models.
- CAMK2A primary H4 = 0.09; sensitivity H4 remains below 0.80.
- Growth-sensitivity results match the residualized GSEA output.

## Logic and claim checks

- Fine-mapping, L2G, Atlas calibration, live scoring and experimental evidence
  are described as distinct evidence layers.
- Model predictions are not presented as experimental measurements.
- Gene-level perturbation is not presented as allele-specific validation.
- PDGFB is not described as externally supported at the target level.
- PRKCQ is not described as strong transcript knockdown.
- IL23R is retained as inconclusive, not negative.
- CAMK2A is retained as non-positive with no colocalization support.
- JAK-STAT and IL-6/JAK/STAT3 retain FDR-supported one-sided tails but
  intervals crossing zero.
- The GSEA sensitivity is not presented as proof that all inflammatory signal
  is absent.

## Remaining limitations

- No new wet-laboratory variant-specific experiment was performed.
- Five targets lack supportive experimental evidence.
- Activity-state pools contain 256 cells and aggregate heterogeneous disease
  and treatment states.
- The frequency-matched background contains 71 variants.
- L2G target assignment remains statistical rather than causal.

After correction, the manuscript's methods, results, tables and interpretive
claims are internally consistent. The principal remaining limitation is
external experimental validation, not an unresolved numerical contradiction.
