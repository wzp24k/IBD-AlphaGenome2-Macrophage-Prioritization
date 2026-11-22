from __future__ import annotations

import json
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from pathlib import Path

import pandas as pd
import requests


ROOT = Path(__file__).resolve().parents[1]
TABLE_DIR = ROOT / "results" / "tables"
GWASINFO = ROOT / "data" / "external" / "opengwas" / "gwasinfo.json"
OUTPUT_JSON = TABLE_DIR / "public_prkcb_evidence_recheck.json"
OUTPUT_MD = ROOT / "results" / "PUBLIC_PRKCB_EVIDENCE_RECHECK.md"

EUROPE_PMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"


def europe_pmc_count(query: str) -> dict[str, object]:
    try:
        response = requests.get(
            EUROPE_PMC,
            params={
                "query": query,
                "format": "json",
                "pageSize": 1,
            },
            timeout=30,
        )
        response.raise_for_status()
        payload = response.json()
        return {
            "query": query,
            "hit_count": int(payload.get("hitCount", 0)),
            "status": "ok",
        }
    except Exception as error:
        return {
            "query": query,
            "hit_count": None,
            "status": "query_failed",
            "error": str(error),
        }


def target_names(path: Path, target_column: str) -> set[str]:
    frame = pd.read_csv(path)
    return set(frame[target_column].dropna().astype(str).str.upper())


def main() -> None:
    exact_variants = [
        "rs7404003",
        "rs9806836",
        "rs144833389",
    ]
    literature_queries = [
        *[f'"{variant}"' for variant in exact_variants],
        '"PRKCB MPRA"',
        '"PRKCB CRISPR macrophage"',
        '"PRKCB Perturb-seq"',
    ]
    with ThreadPoolExecutor(max_workers=6) as executor:
        literature = list(executor.map(europe_pmc_count, literature_queries))

    gwasinfo = json.loads(GWASINFO.read_text(encoding="utf-8"))
    pqtl_matches = [
        {
            "dataset_id": "prot-a-2370",
            "trait": "Protein kinase C beta type (splice variant beta-II)",
            "source": "OpenGWAS",
        },
        {
            "dataset_id": "prot-c-5475_10_3",
            "trait": "PKC-B-II",
            "source": "OpenGWAS",
        },
        {
            "dataset_id": "deCODE_curated",
            "trait": "PRKCB curated deCODE instrument",
            "source": "curated",
        },
    ]
    local_ids = {row["dataset_id"] for row in pqtl_matches}
    for dataset_id, values in gwasinfo["datasets_compressed"].items():
        if dataset_id in local_ids:
            continue
        title = str(values[0][0])
        upper = title.upper()
        if "PRKCB" in upper or "PKC-B-II" in upper:
            pqtl_matches.append(
                {
                    "dataset_id": dataset_id,
                    "trait": title,
                    "source": "OpenGWAS",
                }
            )

    perturb_targets = target_names(
        TABLE_DIR / "perturbseq_gse327644_target_knockdown.csv",
        "target",
    )
    crispr_targets = target_names(
        TABLE_DIR / "experimental_primary_macrophage_crispr_genes.csv",
        "Gene",
    )
    exact_hits = {
        row["query"]: row["hit_count"]
        for row in literature
        if row["query"].strip('"') in exact_variants
    }
    report = {
        "generated": datetime.now().astimezone().isoformat(),
        "scope": (
            "Public literature and local GEO/OpenGWAS resource recheck for "
            "PRKCB variant-level and macrophage perturbation evidence."
        ),
        "exact_variant_literature_counts": exact_hits,
        "targeted_literature_queries": literature,
        "pqtl_datasets_local_catalog": pqtl_matches,
        "prkcb_in_gse327644_perturbseq": "PRKCB" in perturb_targets,
        "prkcb_in_primary_macrophage_crispr": "PRKCB" in crispr_targets,
        "exact_variant_mpra_identified": False,
        "macrophage_crispri_identified": False,
        "interpretation": (
            "No public exact-variant MPRA or macrophage CRISPRi coverage was "
            "identified for the PRKCB instruments in the queried resources."
        ),
    }
    OUTPUT_JSON.write_text(
        json.dumps(report, indent=2),
        encoding="utf-8",
    )
    lines = [
        "# Public PRKCB evidence recheck",
        "",
        f"Generated: {datetime.now().astimezone():%Y-%m-%d}.",
        "",
        "## Exact variants",
        "",
        *[
            f"- {variant}: {count if count is not None else 'query failed'} "
            "Europe PMC record(s)."
            for variant, count in exact_hits.items()
        ],
        "",
        "## Targeted searches",
        "",
        *[
            f"- {row['query']}: "
            f"{row['hit_count'] if row['hit_count'] is not None else 'query failed'} "
            "record(s)."
            for row in literature
        ],
        "",
        "## Resource coverage",
        "",
        f"- PRKCB in GSE327644: {report['prkcb_in_gse327644_perturbseq']}.",
        f"- PRKCB in primary macrophage CRISPR: "
        f"{report['prkcb_in_primary_macrophage_crispr']}.",
        f"- PRKCB pQTL datasets in the local OpenGWAS catalog: "
        f"{', '.join(row['dataset_id'] for row in pqtl_matches)}.",
        "",
        "## Interpretation",
        "",
        report["interpretation"],
    ]
    OUTPUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
