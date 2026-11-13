from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT, configure_proxy, post_graphql


URL = "https://api.platform.opentargets.org/api/v4/graphql"
STUDY_CLUSTERS = {
    "UKB_WGS_2025": [
        "GCST90473822",
        "GCST90473823",
        "GCST90473826",
    ],
    "MVP_2024": [
        "GCST90480318",
        "GCST90480319",
        "GCST90480320",
    ],
}
PAGE_SIZE = 250
QUERY = """
query ExternalL2G($studyIds: [String!]!, $page: Pagination) {
  credibleSets(studyIds: $studyIds, page: $page) {
    count
    rows {
      studyLocusId
      studyId
      l2GPredictions {
        rows {
          score
          target {
            id
            approvedSymbol
            approvedName
          }
        }
      }
    }
  }
}
"""


def main() -> None:
    configure_proxy()
    raw_dir = ROOT / "data" / "raw" / "external_ibd"
    raw_dir.mkdir(parents=True, exist_ok=True)
    session = requests.Session()
    rows = []
    metadata = []
    for cluster, studies in STUDY_CLUSTERS.items():
        for study_id in studies:
            pages = []
            total = None
            page_index = 0
            while total is None or page_index * PAGE_SIZE < total:
                data = post_graphql(
                    session,
                    URL,
                    QUERY,
                    {
                        "studyIds": [study_id],
                        "page": {"index": page_index, "size": PAGE_SIZE},
                    },
                )
                result = data["data"]["credibleSets"]
                total = int(result["count"])
                page_rows = result.get("rows") or []
                pages.append(page_rows)
                if not page_rows:
                    break
                page_index += 1
            payload = {
                "cluster": cluster,
                "studyId": study_id,
                "count": total or 0,
                "pages": pages,
            }
            (raw_dir / f"{study_id}.json").write_text(
                json.dumps(payload, indent=2),
                encoding="utf-8",
            )
            target_count = 0
            for page_rows in pages:
                for locus in page_rows:
                    targets = (locus.get("l2GPredictions") or {}).get("rows") or []
                    for target_prediction in targets:
                        target = target_prediction.get("target") or {}
                        rows.append(
                            {
                                "cluster": cluster,
                                "studyId": study_id,
                                "studyLocusId": locus.get("studyLocusId"),
                                "l2g_score": target_prediction.get("score"),
                                "target_id": target.get("id"),
                                "target_symbol": target.get("approvedSymbol"),
                                "target_name": target.get("approvedName"),
                            }
                        )
                        target_count += 1
            metadata.append(
                {
                    "cluster": cluster,
                    "studyId": study_id,
                    "credible_sets": total or 0,
                    "l2g_rows": target_count,
                }
            )
    output = pd.DataFrame(rows)
    output.to_csv(
        ROOT / "results" / "tables" / "external_ibd_l2g_predictions.csv",
        index=False,
    )
    pd.DataFrame(metadata).to_csv(
        ROOT / "results" / "tables" / "external_ibd_study_summary.csv",
        index=False,
    )
    print(pd.DataFrame(metadata).to_string(index=False))
    print(f"External L2G rows: {len(output)}")


if __name__ == "__main__":
    main()
