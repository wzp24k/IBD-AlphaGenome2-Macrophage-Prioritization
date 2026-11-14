from __future__ import annotations

import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT


MANUSCRIPT = ROOT / "manuscript" / "manuscript.md"
BACKUP = (
    ROOT
    / "manuscript"
    / "archive"
    / "manuscript_pre_reviewer5_edits.md"
)


def main() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")
    replacements = [
        (
            "distinguish model novelty, cell-context specificity and causal evidence",
            "distinguish model prediction, cell-context specificity and causal evidence",
        ),
        (
            "separates genetic plausibility from model novelty, cell-context specificity, independent validation and causal inference",
            "separates genetic plausibility from model prediction, cell-context specificity, independent validation and causal inference",
        ),
        (
            "The top 28 loci showed much stronger stability (rho = 0.83, P < 0.001), indicating that the prioritization was not driven only by absolute-effect outliers.",
            "Among the preselected top 28 loci, correlation was stronger (rho = 0.83, P < 0.001).",
        ),
        (
            "These are neutral results for a cytokine-selection readout, not failed assays or evidence against target-mRNA effects.",
            "These results are neutral for cytokine selection and do not measure target-mRNA effects.",
        ),
        (
            "Classical and intermediate monocytes showed positive PRKCB eQTL lead effects, but exact rs7404003 was significant only in granulocytes (beta = +0.618). The direction was therefore concordant in independent blood myeloid or granulocyte data but not significantly replicated in an independent colon macrophage node.",
            "Classical and intermediate monocytes showed positive PRKCB eQTL lead effects. Exact rs7404003 was significant only in granulocytes (beta = +0.618). CEDAR2 therefore supports the locus-level direction in blood myeloid and granulocyte populations, but not the exact variant in a colon macrophage node.",
        ),
    ]
    updated = text
    counts = []
    for old, new in replacements:
        count = updated.count(old)
        if count == 0:
            raise SystemExit(f"Required reviewer 5 edit not found: {old}")
        counts.append(count)
        updated = updated.replace(old, new)
    BACKUP.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(MANUSCRIPT, BACKUP)
    MANUSCRIPT.write_text(updated, encoding="utf-8")
    print(f"Applied {len(replacements)} reviewer 5 edits: {counts}")


if __name__ == "__main__":
    main()
