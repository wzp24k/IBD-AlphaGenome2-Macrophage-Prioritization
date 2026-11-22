from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT
from remote_tabix import query_region


DATA_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/"
    "sumstats/QTS000035/QTD000584/QTD000584.all.tsv.gz"
)
INDEX_URL = DATA_URL + ".tbi"
CACHE_DIR = ROOT / "work" / "tabix_cache" / "pqtl_sun2018"
TABLE_DIR = ROOT / "results" / "tables"

COLUMNS = [
    "molecular_trait_id",
    "chromosome",
    "position",
    "ref",
    "alt",
    "variant",
    "ma_samples",
    "maf",
    "pvalue",
    "beta",
    "se",
    "type",
    "ac",
    "an",
    "r2",
    "molecular_trait_object_id",
    "gene_id",
    "median_tpm",
    "rsid",
]


def main() -> None:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    lines = query_region(
        DATA_URL,
        INDEX_URL,
        "16",
        23_300_000,
        23_950_000,
        CACHE_DIR,
    )
    rows = [
        line.split("\t")
        for line in lines
        if line.startswith("PRKCB.")
    ]
    if not rows:
        raise RuntimeError("No PRKCB pQTL rows were returned.")

    raw = pd.DataFrame(rows, columns=COLUMNS)
    numeric = ["position", "maf", "pvalue", "beta", "se", "ac", "an", "r2"]
    for column in numeric:
        raw[column] = pd.to_numeric(raw[column], errors="coerce")
    raw["rsid"] = raw["rsid"].astype(str).str.strip()
    raw["variant_id"] = (
        "chr"
        + raw["chromosome"].astype(str)
        + ":"
        + raw["position"].astype(str)
        + ":"
        + raw["ref"]
        + ":"
        + raw["alt"]
    )
    raw = raw.sort_values(["variant_id", "rsid"]).drop_duplicates(
        "variant_id",
        keep="first",
    )
    standardized = pd.DataFrame(
        {
            "variant_id": raw["variant_id"],
            "rsid": raw["rsid"],
            "chromosome": raw["chromosome"].astype(str).radd("chr"),
            "position_grch38": raw["position"].astype(int),
            "other_allele": raw["ref"].astype(str).str.upper(),
            "effect_allele": raw["alt"].astype(str).str.upper(),
            "beta": raw["beta"],
            "se": raw["se"],
            "maf": raw["maf"],
            "p": raw["pvalue"],
            "n_alt_alleles": raw["ac"].astype("Int64"),
            "n_total_alleles": raw["an"].astype("Int64"),
            "imputation_r2": raw["r2"],
            "molecular_trait_id": raw["molecular_trait_id"],
            "gene_id": raw["gene_id"],
            "source": "Sun_2018_eQTL_Catalogue_QTD000584",
            "assembly": "GRCh38",
        }
    ).sort_values(["position_grch38", "variant_id"])
    standardized = standardized[
        standardized["p"].gt(0)
        & standardized["se"].gt(0)
        & standardized["beta"].notna()
    ].reset_index(drop=True)

    output = TABLE_DIR / "prkcb_sun2018_pqtl_cis_grch38.csv"
    standardized.to_csv(output, index=False)
    key_variants = {
        "rs7404003",
        "rs9806836",
        "rs144833389",
        "rs28474709",
        "rs12598807",
    }
    standardized[standardized["rsid"].isin(key_variants)].to_csv(
        TABLE_DIR / "prkcb_sun2018_pqtl_key_variants.csv",
        index=False,
    )

    ld_variants = pd.read_csv(TABLE_DIR / "prkcb_ld_analysis_variants.csv")
    ld = pd.read_csv(TABLE_DIR / "prkcb_ld_matrix.csv", index_col=0)
    ld_rsids = set(ld.index.astype(str))
    matched = standardized[standardized["rsid"].isin(ld_rsids)].copy()
    matched = matched.drop_duplicates("rsid").sort_values("position_grch38")
    matched.to_csv(
        TABLE_DIR / "prkcb_sun2018_pqtl_ld_matched.csv",
        index=False,
    )
    matched_ids = [rsid for rsid in matched["rsid"] if rsid in ld.index]
    ld.loc[matched_ids, matched_ids].rename_axis("rsid").to_csv(
        TABLE_DIR / "prkcb_sun2018_pqtl_ld_matrix.csv"
    )
    ld_variants[ld_variants["rsid"].isin(matched_ids)].to_csv(
        TABLE_DIR / "prkcb_sun2018_pqtl_ld_variants.csv",
        index=False,
    )

    lead = standardized.loc[standardized["p"].idxmin()]
    metadata = {
        "source": "eQTL Catalogue QTS000035 / QTD000584",
        "study": "Sun 2018 plasma pQTL",
        "trait": "PRKCB.5475.10.3..1",
        "trait_object": "ENSG00000166501",
        "assembly": "GRCh38",
        "requested_region": "chr16:23,300,000-23,950,000",
        "observed_region": {
            "start": int(standardized["position_grch38"].min()),
            "end": int(standardized["position_grch38"].max()),
        },
        "samples": 3301,
        "variants": len(standardized),
        "unique_rsids": int(standardized["rsid"].nunique()),
        "ld_matched_variants": len(matched_ids),
        "lead_variant": {
            "variant_id": lead["variant_id"],
            "rsid": lead["rsid"],
            "position_grch38": int(lead["position_grch38"]),
            "other_allele": lead["other_allele"],
            "effect_allele": lead["effect_allele"],
            "beta": float(lead["beta"]),
            "se": float(lead["se"]),
            "maf": float(lead["maf"]),
            "p": float(lead["p"]),
        },
        "output": str(output),
        "ld_matrix": str(TABLE_DIR / "prkcb_sun2018_pqtl_ld_matrix.csv"),
    }
    (TABLE_DIR / "prkcb_sun2018_pqtl_metadata.json").write_text(
        json.dumps(metadata, indent=2),
        encoding="utf-8",
    )
    readme = (
        "# Sun 2018 PRKCB cis-pQTL\n\n"
        "## Source\n\n"
        "- Resource: eQTL Catalogue QTS000035 / QTD000584.\n"
        "- Study: Sun 2018 plasma pQTL.\n"
        "- Molecular trait: PRKCB.5475.10.3..1.\n"
        "- Samples: 3,301.\n"
        "- Assembly: GRCh38.\n"
        "- Requested region: chr16:23,300,000-23,950,000.\n"
        f"- Observed region: chr16:{metadata['observed_region']['start']}-"
        f"{metadata['observed_region']['end']}.\n\n"
        "## Files\n\n"
        "- `prkcb_sun2018_pqtl_cis_grch38.csv`: standardized PRKCB pQTL rows.\n"
        "- `prkcb_sun2018_pqtl_ld_matched.csv`: pQTL rows with regional LD.\n"
        "- `prkcb_sun2018_pqtl_ld_matrix.csv`: Pearson r matrix for matched variants.\n"
        "- `prkcb_sun2018_pqtl_ld_variants.csv`: variant annotation for the LD matrix.\n"
        "- `prkcb_sun2018_pqtl_mr_ibd.csv`: one-instrument cis-MR estimates.\n"
        "- `prkcb_sun2018_pqtl_coloc.csv`: ABF colocalization estimates.\n\n"
        "## Required Columns\n\n"
        "`variant_id`, `rsid`, `chromosome`, `position_grch38`, "
        "`other_allele`, `effect_allele`, `beta`, `se`, `maf`, and `p`.\n\n"
        "## LD\n\n"
        "The LD matrix uses the existing 1000 Genomes high-coverage European "
        f"reference. {metadata['ld_matched_variants']} pQTL variants matched "
        "the regional reference. Matrix values are Pearson r, not r-squared.\n\n"
        "## Lead Signal\n\n"
        f"- Variant: {metadata['lead_variant']['variant_id']} "
        f"({metadata['lead_variant']['rsid']}).\n"
        f"- Beta: {metadata['lead_variant']['beta']:.6f}.\n"
        f"- SE: {metadata['lead_variant']['se']:.6f}.\n"
        f"- P: {metadata['lead_variant']['p']:.3e}.\n"
    )
    (TABLE_DIR / "prkcb_sun2018_pqtl_README.md").write_text(
        readme,
        encoding="utf-8",
    )
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
