from __future__ import annotations

import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT, configure_proxy, load_config, post_graphql


OPEN_TARGETS_URL = "https://api.platform.opentargets.org/api/v4/graphql"
QUERY = """
query Variant($id: String!) {
  variant(variantId: $id) {
    id
    alleleFrequencies {
      populationName
      alleleFrequency
    }
  }
}
"""


def fetch_variant(
    session: requests.Session,
    canonical_variant: str,
    open_targets_id: str,
) -> dict:
    result = post_graphql(
        session,
        OPEN_TARGETS_URL,
        QUERY,
        {"id": open_targets_id},
        attempts=3,
    )
    payload = result["data"]["variant"] or {}
    frequencies = {
        row["populationName"]: row["alleleFrequency"]
        for row in payload.get("alleleFrequencies") or []
        if row.get("alleleFrequency") is not None
    }
    values = list(frequencies.values())
    af_max = max(values) if values else None
    return {
        "canonical_variant": canonical_variant,
        "open_targets_variant_id": open_targets_id,
        "frequency_found": bool(payload),
        "frequency_source": "Open Targets harmonized population AF",
        "af_max": af_max,
        "af_matched": min(af_max, 1 - af_max) if af_max is not None else None,
        "af_nfe": frequencies.get("nfe_adj"),
        "af_eas": frequencies.get("eas_adj"),
        "af_afr": frequencies.get("afr_adj"),
        "af_sas": frequencies.get("sas_adj"),
        "population_count": len(frequencies),
    }


def main() -> None:
    configure_proxy()
    config = load_config()
    path = ROOT / "data" / "processed" / "ibd_noncoding_credible_variants.csv"
    candidates = pd.read_csv(path)
    primary = candidates[candidates["study_role"] == "primary_ibd"][
        ["canonical_variant", "variant_id"]
    ]
    positive = pd.DataFrame(config["positive_controls"])[
        ["variant", "id"]
    ].rename(
        columns={"variant": "canonical_variant", "id": "variant_id"}
    )
    variants = (
        pd.concat([primary, positive], ignore_index=True)
        .drop_duplicates("canonical_variant")
        .to_dict("records")
    )

    session = requests.Session()
    records = []
    with ThreadPoolExecutor(max_workers=8) as executor:
        futures = {
            executor.submit(
                fetch_variant,
                session,
                variant["canonical_variant"],
                variant["variant_id"],
            ): variant
            for variant in variants
        }
        for completed, future in enumerate(as_completed(futures), start=1):
            variant = futures[future]
            try:
                records.append(future.result())
            except Exception as exc:
                records.append(
                    {
                        "canonical_variant": variant["canonical_variant"],
                        "open_targets_variant_id": variant["variant_id"],
                        "frequency_found": False,
                        "frequency_error": str(exc),
                    }
                )
            if completed % 50 == 0 or completed == len(variants):
                print(f"Frequency annotation: {completed}/{len(variants)}")

    frame = pd.DataFrame(records).sort_values("canonical_variant")
    frame["af_common"] = frame["af_matched"].fillna(0) >= 0.01
    output_dir = ROOT / "results" / "tables"
    output_dir.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output_dir / "primary_variants_frequency.csv", index=False)
    print(frame["frequency_found"].value_counts(dropna=False).to_string())
    print(frame["af_common"].value_counts(dropna=False).to_string())


if __name__ == "__main__":
    main()
