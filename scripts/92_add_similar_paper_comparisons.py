from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

from docx import Document
from docx.oxml.ns import qn
from docx.shared import Pt

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT


DOCX = ROOT / "1119manuscript.docx"
COPY = ROOT / "manuscript" / "Manuscript_AlphaGenome2_IBD_clinical_translation.docx"
HELPER_PATH = ROOT / "scripts" / "25_build_template_docx.py"

COMPARISONS = [
    (
        "Cell-type-resolved genetic variation shapes inflammatory bowel disease risk",
        "Alegbe et al. mapped cis-eQTLs across intestinal and blood cell types and "
        "showed that cell-type-resolved regulation increases colocalization with "
        "IBD loci and nominates myeloid effector genes [43]. Our framework "
        "addresses the same problem at macrophage-credible-set resolution and "
        "adds protein-level and treatment-response evidence, which places PRKCB "
        "in a development-assessment tier.",
    ),
    (
        "Prioritizing Crohn's disease genes by integrating association signals with gene expression implicates monocyte subsets",
        "Gettler et al. integrated Crohn's association signals with bulk ileal and "
        "stimulated-monocyte expression to prioritize genes and monocyte subsets "
        "[44]. Their approach operates mainly at gene level. The present analysis "
        "resolves the PRKCB locus through eQTL/pQTL MR, SuSiE colocalization and a "
        "clinical-readiness matrix, which supports experimental validation and "
        "drug development assessment.",
    ),
]

REFERENCE_43 = (
    "43. Alegbe T, Harris BT, Fachal L, et al. Cell-type-resolved genetic "
    "variation shapes inflammatory bowel disease risk. Nature. 2026;656:129-139. "
    "doi:10.1038/s41586-026-10627-z."
)
REFERENCE_44 = (
    "44. Gettler K, Giri M, Kenigsberg E, et al. Prioritizing Crohn's disease "
    "genes by integrating association signals with gene expression implicates "
    "monocyte subsets. Genes & Immunity. 2019;20:577-588. "
    "doi:10.1038/s41435-019-0059-y."
)


def load_helper():
    spec = importlib.util.spec_from_file_location("template_builder", HELPER_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def set_paragraph_format(paragraph, size: float = 11) -> None:
    paragraph.paragraph_format.space_before = Pt(0)
    paragraph.paragraph_format.space_after = Pt(0)
    paragraph.paragraph_format.line_spacing = 1.5
    for run in paragraph.runs:
        run.font.name = "Times New Roman"
        run.font.size = Pt(size)


def main() -> None:
    document = Document(DOCX)
    existing = "\n".join(p.text for p in document.paragraphs)
    if "ref_43" in {
        bookmark.get(qn("w:name"))
        for bookmark in document.element.body.iter(qn("w:bookmarkStart"))
    } or REFERENCE_43 in existing:
        print("Similar-paper comparisons already present.")
        return

    helper = load_helper()
    limitations = next(
        paragraph
        for paragraph in document.paragraphs
        if paragraph.text.strip() == "4.4 Limitations"
    )
    for _, text in COMPARISONS:
        paragraph = limitations.insert_paragraph_before()
        helper.add_inline(paragraph, text)
        set_paragraph_format(paragraph)

    reference_numbers = {
        int(match.group(1))
        for paragraph in document.paragraphs
        if (match := __import__("re").match(r"^(\d+)\.", paragraph.text))
    }
    next_bookmark_id = max(
        int(bookmark.get(qn("w:id")))
        for bookmark in document.element.body.iter(qn("w:bookmarkStart"))
        if bookmark.get(qn("w:id")) is not None
    ) + 1
    for number, text in [(43, REFERENCE_43), (44, REFERENCE_44)]:
        paragraph = document.add_paragraph()
        run = paragraph.add_run(text)
        run.font.name = "Times New Roman"
        run.font.size = Pt(10)
        helper.add_reference_bookmark(paragraph, number, next_bookmark_id)
        next_bookmark_id += 1

    document.save(DOCX)
    COPY.write_bytes(DOCX.read_bytes())
    revised = Document(DOCX)
    text = "\n".join(p.text for p in revised.paragraphs)
    checks = {
        "alegbe_comparison": "Alegbe et al." in text,
        "gettler_comparison": "Gettler et al." in text,
        "reference_43": "doi:10.1038/s41586-026-10627-z" in text,
        "reference_44": "doi:10.1038/s41435-019-0059-y" in text,
    }
    print(checks)
    if not all(checks.values()):
        raise SystemExit("Similar-paper comparison insertion failed.")


if __name__ == "__main__":
    main()
