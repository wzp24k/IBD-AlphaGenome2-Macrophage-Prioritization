from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT, configure_proxy, load_config, post_graphql


OPEN_TARGETS_URL = "https://api.platform.opentargets.org/api/v4/graphql"
PAGE_SIZE = 250

QUERY = """
query CredibleSets($studyIds: [String!]!, $page: Pagination) {
  credibleSets(studyIds: $studyIds, page: $page) {
    count
    rows {
      studyLocusId
      studyId
      credibleSetIndex
      confidence
      finemappingMethod
      beta
      standardError
      pValueMantissa
      pValueExponent
      l2GPredictions {
        count
        rows {
          score
          target {
            id
            approvedSymbol
            approvedName
          }
        }
      }
      locus {
        count
        rows {
          posteriorProbability
          is95CredibleSet
          is99CredibleSet
          r2Overall
          variant {
            id
            chromosome
            position
            referenceAllele
            alternateAllele
            rsIds
            mostSevereConsequence {
              id
              label
            }
          }
        }
      }
    }
  }
}
"""


def flatten_rows(page_rows: list[dict]) -> list[dict]:
    records: list[dict] = []
    for outer in page_rows:
        locus = outer.get("locus") or {}
        l2g_rows = ((outer.get("l2GPredictions") or {}).get("rows")) or []
        top_l2g = max(l2g_rows, key=lambda row: row.get("score", 0), default={})
        l2g_target = top_l2g.get("target") or {}
        for inner in locus.get("rows") or []:
            variant = inner.get("variant") or {}
            consequence = variant.get("mostSevereConsequence") or {}
            records.append(
                {
                    "studyLocusId": outer.get("studyLocusId"),
                    "studyId": outer.get("studyId"),
                    "credibleSetIndex": outer.get("credibleSetIndex"),
                    "confidence": outer.get("confidence"),
                    "finemappingMethod": outer.get("finemappingMethod"),
                    "study_beta": outer.get("beta"),
                    "study_standardError": outer.get("standardError"),
                    "pValueMantissa": outer.get("pValueMantissa"),
                    "pValueExponent": outer.get("pValueExponent"),
                    "l2g_score": top_l2g.get("score"),
                    "l2g_target_id": l2g_target.get("id"),
                    "l2g_target_symbol": l2g_target.get("approvedSymbol"),
                    "l2g_target_name": l2g_target.get("approvedName"),
                    "pip": inner.get("posteriorProbability"),
                    "is95CredibleSet": inner.get("is95CredibleSet"),
                    "is99CredibleSet": inner.get("is99CredibleSet"),
                    "r2Overall": inner.get("r2Overall"),
                    "variant_id": variant.get("id"),
                    "chromosome": variant.get("chromosome"),
                    "position": variant.get("position"),
                    "referenceAllele": variant.get("referenceAllele"),
                    "alternateAllele": variant.get("alternateAllele"),
                    "rsIds": ",".join(variant.get("rsIds") or []),
                    "consequence_id": consequence.get("id"),
                    "most_severe_consequence": consequence.get("label"),
                }
            )
    return records


def main() -> None:
    configure_proxy()
    config = load_config()
    raw_dir = ROOT / "data" / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)

    session = requests.Session()
    all_records: list[dict] = []
    summary: list[dict] = []

    for study in config["studies"]:
        study_id = study["id"]
        raw_path = raw_dir / f"open_targets_v2_{study_id}.json"
        records_for_study: list[dict] = []
        total = None
        page_index = 0

        if raw_path.exists():
            payload = json.loads(raw_path.read_text(encoding="utf-8"))
            pages = payload["pages"]
        else:
            pages = []
            while total is None or page_index * PAGE_SIZE < total:
                data = post_graphql(
                    session,
                    OPEN_TARGETS_URL,
                    QUERY,
                    {
                        "studyIds": [study_id],
                        "page": {"index": page_index, "size": PAGE_SIZE},
                    },
                )
                result = data["data"]["credibleSets"]
                total = int(result["count"])
                rows = result.get("rows") or []
                pages.append(rows)
                if not rows:
                    break
                page_index += 1
            payload = {"study": study, "count": total, "pages": pages}
            raw_path.write_text(
                json.dumps(payload, indent=2, ensure_ascii=True),
                encoding="utf-8",
            )

        for page_rows in pages:
            records_for_study.extend(flatten_rows(page_rows))

        if not records_for_study:
            summary.append(
                {
                    "study_id": study_id,
                    "role": study["role"],
                    "credible_set_count": total or 0,
                    "variant_row_count": 0,
                    "locus_count": 0,
                }
            )
            continue

        frame = pd.DataFrame(records_for_study)
        frame["study_role"] = study["role"]
        frame["study_label"] = study["label"]
        all_records.extend(frame.to_dict("records"))
        summary.append(
            {
                "study_id": study_id,
                "role": study["role"],
                "credible_set_count": total or 0,
                "variant_row_count": len(frame),
                "locus_count": frame["studyLocusId"].nunique(),
            }
        )
        print(
            f"{study_id}: {total} credible sets, "
            f"{len(frame)} variant rows, {frame['studyLocusId'].nunique()} loci"
        )

    output = pd.DataFrame(all_records)
    output.to_csv(raw_dir / "open_targets_credible_sets_v2.csv", index=False)
    pd.DataFrame(summary).to_csv(raw_dir / "open_targets_study_summary.csv", index=False)
    print(f"Wrote {len(output)} variant rows to data/raw.")


if __name__ == "__main__":
    main()
