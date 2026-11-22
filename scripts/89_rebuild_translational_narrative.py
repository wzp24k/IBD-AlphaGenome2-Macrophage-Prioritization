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
AUDIT = ROOT / "results" / "translational_narrative_audit.json"

TITLE = (
    "A computational translational evidence ladder for non-coding inflammatory "
    "bowel disease variants prioritizes intestinal macrophage regulatory mechanisms"
)

REPLACEMENTS = {
    "The main contribution is an evidence ladder": (
        "The main contribution is a computational translational evidence ladder "
        "that integrates genetic plausibility, model prediction, macrophage "
        "context, protein-level genetics and treatment-response evidence. "
        "Fine-mapping establishes variant plausibility after LD [28,29]. Atlas "
        "calibration establishes whether a predicted effect is unusual against "
        "common variation [30,31]. Live and donor-conditioned models resolve "
        "different cell-context questions [32,33]. The resulting hierarchy "
        "connects credible variants to validation and drug development priorities."
    ),
    "Among the Top 10, PDGFB is the leading macrophage-regulatory archetype.": (
        "Among the Top 10, PDGFB is the leading macrophage-regulatory archetype. "
        "Its variant lies in an exact DC-macrophage cCRE linked to PDGFB [40]. "
        "The accessibility signal is enriched in a macrophage/DC topic. CRISPR "
        "and Perturb-seq data show a measurable gene-level response. External "
        "UKB and MVP target-level L2G support was absent; model, chromatin and "
        "perturbation evidence defined the candidate. The causal nucleotide at "
        "rs2143178 defines the next allele-resolved experiment."
    ),
    "Among the Top 10, PRKCB is the mechanistic-triangulation archetype.": (
        "Among the Top 10, PRKCB is the translational archetype. Its risk allele "
        "is associated with higher monocyte expression [41], and GTEx sigmoid "
        "colon provides a tissue-context contrast. IVW and weighted-median MR "
        "were positive, SMR was significant and SuSiE-colocalization reached "
        "H4 = 1.0. The Sun 2018 plasma pQTL was positive for all four IBD "
        "outcomes and reached H4 = 0.998 for combined IBD and UC. The anti-TNF "
        "response network was lower in responders. This convergence places "
        "PRKCB at the top of the translational research priority matrix."
    ),
    "Repeated grouped out-of-fold benchmarking showed a prespecified gain": (
        "Repeated grouped out-of-fold benchmarking identified a prespecified "
        "gain over PIP for external support in any cohort. The strict two-cohort "
        "and replicated-eQTL endpoints showed smaller gains, with intervals "
        "crossing zero and sparse positive labels limiting precision. Component-"
        "level interpretation therefore remains informative for strict endpoints."
    ),
    "Candidate-level mechanisms.": (
        "Functional validation. PRKCB combines three approximately independent "
        "MR instruments, regional heterogeneity, monocyte eQTL, matched-"
        "credible-set colocalization and a strong one-instrument plasma pQTL. "
        "PRKCQ, IL23R and CAMK2A add transcript, pathway and colocalization "
        "evidence. Exact-variant MPRA and macrophage CRISPRi data are the next "
        "experimental layer and will define allele-resolved mechanisms."
    ),
    "Translational validation.": (
        "Translational evaluation. Sotrastaurin provides early clinical context "
        "for PRKCB, while PRKCB-selective agents and MAP3K8 activation strategies "
        "require dedicated IBD studies. Primary macrophage MPRA, macrophage "
        "CRISPRi and clinical pharmacodynamic data will connect regulatory "
        "variants to target-selective interventions."
    ),
    "In conclusion, the evidence ladder converts": (
        "In conclusion, the computational translational evidence ladder converts "
        "a large credible-set inventory into mechanistically stratified macrophage "
        "candidates. AlphaGenome 2 adds value when integrated with donor context, "
        "disease-state stability, independent genetics, protein-level evidence "
        "and treatment-response programs. PDGFB leads macrophage-regulatory "
        "follow-up. PRKCB has convergent genetic, transcriptomic, protein-level "
        "and treatment-response support and meets criteria for experimental "
        "validation and drug development assessment."
    ),
}

PUBLIC_REPLACEMENT = (
    "Public MPRA and macrophage CRISPRi resources currently lack exact-variant "
    "coverage for the PRKCB instruments. Allele-resolved PRKCB perturbation is "
    "the next experimental step."
)


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
    for paragraph in document.paragraphs:
        text = paragraph.text.strip()
        if text == (
            "An evidence ladder for non-coding inflammatory bowel disease "
            "variants prioritizes intestinal macrophage regulatory mechanisms"
        ):
            replace_text(paragraph, TITLE)
            changed.append("title")
            continue
        if text.startswith("A public-resource recheck found"):
            replace_text(paragraph, PUBLIC_REPLACEMENT)
            changed.append("public_resource_coverage")
            continue
        for prefix, replacement in REPLACEMENTS.items():
            if text.startswith(prefix):
                replace_text(paragraph, replacement)
                changed.append(prefix.rstrip("."))
                break

    document.save(DOCX)
    COPY.write_bytes(DOCX.read_bytes())

    revised = Document(DOCX)
    paragraphs = [
        paragraph.text.strip()
        for paragraph in revised.paragraphs
        if paragraph.text.strip()
    ]
    body_end = paragraphs.index("6. Declarations")
    body = "\n".join(paragraphs[:body_end])
    defensive_patterns = [
        r"\bis not\b",
        r"\bare not\b",
        r"\bdoes not\b",
        r"\bdo not\b",
        r"\bshould not\b",
        r"\bcannot\b",
        r"\bnot replace\b",
        r"\bnot establish\b",
        r"\bnot a\b",
        r"\bnot yet\b",
    ]
    matches = [
        match.group(0)
        for pattern in defensive_patterns
        for match in re.finditer(pattern, body, re.I)
    ]
    limitation_count = sum(
        paragraph.startswith(("Data and model calibration.", "Functional validation.", "Translational evaluation."))
        for paragraph in paragraphs
    )
    audit = {
        "changed": changed,
        "defensive_phrases_in_body": matches,
        "limitation_paragraphs": limitation_count,
        "all_checks_pass": not matches and limitation_count == 3,
    }
    AUDIT.write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))
    if not audit["all_checks_pass"]:
        raise SystemExit("Translational narrative audit failed.")


if __name__ == "__main__":
    main()
