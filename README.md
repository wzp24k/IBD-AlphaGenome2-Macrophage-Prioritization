# IBD non-coding variant prioritization with AlphaGenome 2

This repository contains the analysis code, configuration, figure files and
aggregate result tables for a study of intestinal macrophage-focused
prioritization of non-coding inflammatory bowel disease variants.

## Scope

- Primary discovery study: GCST90292538.
- Primary candidate set: 395 non-coding variants from 281 independent credible sets.
- Regulatory inputs: Open Targets L2G, AlphaGenome Atlas and live AlphaGenome 2.
- External genetic clusters: UK Biobank WGS 2025 and MVP 2024.
- Public experimental evidence: GTEx, eQTL Catalogue, GSE229472, GSE327124,
  GSE327644, GSE175685 and GSE273194.
- Output: 10 Tier 1 and 57 Tier 2 loci with explicit evidence boundaries.

## Repository contents

- `scripts/`: numbered analysis, figure and audit scripts.
- `config/`: analysis configuration.
- `results/figures/`: final Figure 1-10 files, including the editable
  `figure1_analysis_flow.drawio` source.
- `results/tables/`: aggregate results smaller than 20 MB.
- `results/*.md`: evidence and quality-control audits.
- `run_core_pipeline.ps1`: pipeline entry point for a prepared local data tree.

The manuscript citation audit verifies sequential use of references 1-42,
one citation per reference, ascending citation groups and Word bookmarks for
every entry.

The 37 MB Perturb-seq pseudobulk matrix is excluded from version control.
Target-knockdown summaries, differential-expression outputs and robustness
tables remain available in `results/tables/`.

## Data boundary

This repository does not redistribute raw credentialed data, large single-cell
objects, raw GEO archives, the manuscript or API credentials. Public source data
must be obtained from Open Targets, Ensembl, AlphaGenome, CELLxGENE, NCBI GEO,
GTEx and the cited GWAS resources under their access terms.

## Reproduction

Create `API.txt` locally with a valid AlphaGenome API key or set
`ALPHAGENOME_API_KEY`. Prepare the external and processed data paths described
in `config/analysis.json`, then run:

```powershell
powershell -ExecutionPolicy Bypass -File .\run_core_pipeline.ps1
```

The pipeline retains intermediate summaries and excludes API credentials from
all outputs. Large raw datasets are intentionally absent from version control.

## Citation

The manuscript citing this repository is not included here. Please cite the
associated article when available.
