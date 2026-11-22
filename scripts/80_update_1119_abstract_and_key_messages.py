from __future__ import annotations

import json
import re
import sys
from pathlib import Path

from docx import Document
from docx.oxml.ns import qn
from docx.shared import Pt

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT


DOCX = ROOT / "1119manuscript.docx"
COPY = ROOT / "manuscript" / "Manuscript_AlphaGenome2_IBD_clinical_translation.docx"
AUDIT = ROOT / "results" / "1119_abstract_key_message_audit.json"


def replace_text(paragraph, text: str) -> None:
    first = paragraph.runs[0] if paragraph.runs else None
    bold = first.bold if first else None
    size = first.font.size if first else Pt(11)
    for child in list(paragraph._p):
        if child.tag != qn("w:pPr"):
            paragraph._p.remove(child)
    run = paragraph.add_run(text)
    run.font.name = "Times New Roman"
    run.font.size = size
    if bold is not None:
        run.bold = bold


def main() -> None:
    document = Document(DOCX)
    changed = []
    repeated_sentences = [
        "Cross-ancestry, pQTL and treatment-response analyses extended the clinical translation layer.",
        "Cross-ancestry MR and treatment-response programs ranked PRKCB highest for translational readiness.",
        "A strong cis-pQTL supported PRKCB protein in combined IBD and UC; cross-ancestry MR and treatment-response programs ranked it highest.",
        "Cross-ancestry, pQTL and treatment-response analyses identify PRKCB as the leading translation candidate.",
    ]
    for paragraph in document.paragraphs:
        text = paragraph.text.strip()
        for sentence in repeated_sentences:
            while f"{sentence} {sentence}" in text:
                text = text.replace(f"{sentence} {sentence}", sentence)
        if text.startswith("Methods:"):
            replacement = text
            replacement = replacement.replace(
                "Cross-ancestry, pQTL and treatment-response analyses extended "
                "the clinical translation layer.",
                "Cross-ancestry MR, a strong cis-pQTL and treatment-response "
                "programs extended the translational layer.",
            )
            if "strong cis-pQTL" not in replacement:
                replacement += (
                    " Cross-ancestry MR, a strong cis-pQTL and treatment-response "
                    "programs extended the translational layer."
                )
            replace_text(paragraph, replacement)
            changed.append("abstract_methods")
        elif text.startswith("Background:"):
            replacement = text.replace(
                "A useful prioritization system integrates model prediction, "
                "cell-context specificity and causal evidence.",
                "Prioritization integrates model prediction, cell-context "
                "specificity and translational evidence.",
            )
            replace_text(paragraph, replacement)
            changed.append("abstract_background")
        elif text.startswith("Results:"):
            replacement = text.replace(
                "Virtual knockdown predicted lower active-macrophage JAK/STAT "
                "activity, and signature reversal was strongest in Crohn and "
                "Garrido macrophages.",
                "Virtual knockdown supported JAK/STAT modulation.",
            )
            replacement = replacement.replace(
                "Cross-ancestry MR and treatment-response programs ranked PRKCB "
                "highest for translational readiness.",
                "A strong cis-pQTL supported PRKCB protein in combined IBD and UC; "
                "cross-ancestry MR and treatment-response programs ranked it highest.",
            )
            replacement = replacement.replace(
                "Its three instruments and regional heterogeneity defined the scope "
                "of the mediation signal.",
                "Three instruments and regional heterogeneity defined the mediation signal.",
            )
            if "strong cis-pQTL" not in replacement:
                replacement += (
                    " A strong cis-pQTL supported PRKCB protein in combined IBD and "
                    "UC; cross-ancestry MR and treatment-response programs ranked it "
                    "highest."
                )
            replace_text(paragraph, replacement)
            changed.append("abstract_results")
        elif text.startswith("The framework converts genetic credible sets"):
            replacement = text
            replacement = replacement.replace(
                "Cross-ancestry, pQTL and treatment-response analyses identify "
                "PRKCB as the leading translation candidate.",
                "A strong cis-pQTL, cross-ancestry MR and treatment-response "
                "programs make PRKCB the leading candidate for experimental "
                "validation and drug development assessment.",
            )
            if "leading candidate for experimental validation" not in replacement:
                replacement += (
                    " A strong cis-pQTL, cross-ancestry MR and treatment-response "
                    "programs make PRKCB the leading candidate for experimental "
                    "validation and drug development assessment."
                )
            replace_text(paragraph, replacement)
            changed.append("key_message_3")
        elif text.startswith("Conclusion:"):
            replace_text(
                paragraph,
                "Conclusion: The evidence ladder converts a large credible-set "
                "inventory into mechanistically stratified macrophage candidates. "
                "PDGFB leads macrophage-regulatory follow-up, and PRKCB meets "
                "criteria for experimental validation and drug development assessment.",
            )
            changed.append("abstract_conclusion")

    document.save(DOCX)
    COPY.write_bytes(DOCX.read_bytes())
    revised = Document(DOCX)
    paragraphs = [paragraph.text.strip() for paragraph in revised.paragraphs]
    abstract_start = paragraphs.index("Abstract")
    keywords_index = next(
        index for index, text in enumerate(paragraphs)
        if text.startswith("Keywords:")
    )
    abstract_text = " ".join(paragraphs[abstract_start + 1 : keywords_index])
    abstract_words = len(
        re.findall(r"[A-Za-z0-9][A-Za-z0-9'./-]*", abstract_text)
    )
    checks = {
        "methods_updated": "strong cis-pQTL" in "\n".join(paragraphs),
        "results_updated": "supported PRKCB protein in combined IBD and UC" in "\n".join(paragraphs),
        "key_message_updated": "leading candidate for experimental validation" in "\n".join(paragraphs),
        "conclusion_updated": "PRKCB meets criteria for experimental validation" in "\n".join(paragraphs),
        "abstract_250_or_less": abstract_words <= 250,
    }
    audit = {
        "document": str(DOCX),
        "changed": changed,
        "abstract_words": abstract_words,
        "all_checks_pass": all(checks.values()),
        "checks": checks,
    }
    AUDIT.write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))
    if not audit["all_checks_pass"]:
        raise SystemExit("Abstract or Key Messages update failed.")


if __name__ == "__main__":
    main()
