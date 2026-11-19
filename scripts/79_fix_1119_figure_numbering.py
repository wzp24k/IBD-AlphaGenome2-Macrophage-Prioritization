from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

from docx import Document
from docx.oxml.ns import qn
from docx.shared import Pt

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT


DOCX = ROOT / "1119manuscript.docx"
COPY = ROOT / "manuscript" / "Manuscript_AlphaGenome2_IBD_clinical_translation.docx"
BACKUP = ROOT / "_archive" / "1119manuscript_before_figure_number_fix.docx"
AUDIT = ROOT / "results" / "1119_figure_number_fix_audit.json"


def replace_paragraph(paragraph, text: str) -> None:
    first_run = paragraph.runs[0] if paragraph.runs else None
    bold = first_run.bold if first_run else None
    size = first_run.font.size if first_run else Pt(11)
    for child in list(paragraph._p):
        if child.tag != qn("w:pPr"):
            paragraph._p.remove(child)
    run = paragraph.add_run(text)
    run.font.name = "Times New Roman"
    run.font.size = size
    if bold is not None:
        run.bold = bold


def main() -> None:
    BACKUP.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(DOCX, BACKUP)
    document = Document(DOCX)
    replacements = {
        "PRKCQ repeated-pool effects": (
            "Among the top 28 loci, 21 were concordant with the original active "
            "pool and at least one independent active-disease cohort. Seventeen "
            "were concordant in both Crohn and Garrido cohorts. Nine Tier 1 loci "
            "were persistent in the initial analysis and one was attenuated. "
            "PRKCQ repeated-pool effects were near zero in active colitis (-0.00), "
            "active Crohn macrophages (-0.00) and active IBD myeloid cells (+0.00). "
            "These values were below the 0.005 interpretability threshold, while "
            "healthy Oliver macrophages remained negative (median -0.01; Table 5; "
            "Figure 5). The repeated pools quantify sampling stability within "
            "available donors. PRKCQ was therefore classified as attenuated. "
            "ABCB5, IRF2BP2 and IL23R were additionally marked target-uncertain "
            "because L2G was below 0.50."
        ),
        "The external background contained": (
            "The external background contained 71 frequency-matched common SNVs "
            "from GRCh38 regions at least 1 Mb from primary IBD loci (Figure 6). "
            "Empirical tail probabilities provided calibration statistics, and the "
            "minimum attainable one-sided probability was approximately 0.014. "
            "PDGFB had macrophage-effect and specificity tail probabilities of "
            "0.014 and 0.056. PRKCB, CHSY3, PRKCQ, ABCB5 and IRF2BP2 also had "
            "macrophage-effect tail probabilities of 0.014, while other Tier 1 "
            "candidates had weaker background-tail support. Evidence components "
            "were retained separately to provide component-level ranking."
        ),
        "Figure 6. Activity-state evaluation": (
            "Figure 5. Activity-state evaluation of prioritized loci. (A) Predicted "
            "target effects across healthy and active intestinal macrophage "
            "conditions. (B) Healthy-macrophage effect versus active-colitis effect "
            "for the top 28 loci. (C) Six resampled donor-balanced PRKCQ pools per "
            "condition. The pools quantify sampling stability across donor-balanced "
            "resamples."
        ),
    }
    changed = []
    for paragraph in document.paragraphs:
        text = paragraph.text.strip()
        for prefix, replacement in replacements.items():
            if prefix in text:
                replace_paragraph(paragraph, replacement)
                changed.append(prefix)
                break
    document.save(DOCX)
    COPY.write_bytes(DOCX.read_bytes())
    final_text = "\n".join(p.text for p in document.paragraphs)
    checks = {
        "prkcq_reference_fixed": "Table 5; Figure 5" in final_text,
        "background_reference_fixed": "primary IBD loci (Figure 6)" in final_text,
        "activity_caption_fixed": "Figure 5. Activity-state evaluation" in final_text,
        "figure_caption_sequence_has_five_six": "Figure 5. Activity-state evaluation" in final_text
        and "Figure 6. Empirical background distribution" in final_text,
    }
    audit = {
        "document": str(DOCX),
        "backup": str(BACKUP),
        "changed": changed,
        "all_checks_pass": all(checks.values()),
        "checks": checks,
    }
    AUDIT.write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))
    if not audit["all_checks_pass"]:
        raise SystemExit("Figure numbering fix failed.")


if __name__ == "__main__":
    main()
