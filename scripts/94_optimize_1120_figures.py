from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
import sys
import zipfile
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont, ImageOps

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT


DOCX = ROOT / "1120manuscript.docx"
BACKUP = ROOT / "_archive" / "1120manuscript_before_figure_optimization.docx"
FIGURE_DIR = ROOT / "results" / "figures"
WORK_DIR = ROOT / "work" / "1120_figure_optimization"
AUDIT = ROOT / "results" / "1120_figure_optimization_audit.json"

MEDIA_TO_FIGURE = {
    "image1.png": "figure1_analysis_flow.png",
    "image2.png": "figure4_ets2_positive_control.png",
    "image3.png": "figure2_quantile_vs_specificity.png",
    "image4.png": "figure3_top_candidate_heatmap.png",
    "image5.png": "figure6_active_ibd_validation.png",
    "image6.png": "figure5_ets2_background_null.png",
    "image7.png": "figure7_all_l2g_pathways.png",
    "image8.png": "figure8_external_experimental_evidence.png",
    "image9.png": "figure9_perturbseq_validation.png",
    "image10.png": "figure10_evidence_robustness.png",
    "image11.png": "figure11_benchmark_translation.png",
    "image12.png": "figure12_prkcb_mechanism.png",
    "image13.png": "figure13_clinical_translation.png",
}

FIGURE_TO_MEDIA = {
    "figure1_analysis_flow.png": "image1.png",
    "figure2_quantile_vs_specificity.png": "image3.png",
    "figure3_top_candidate_heatmap.png": "image4.png",
    "figure4_ets2_positive_control.png": "image2.png",
    "figure5_ets2_background_null.png": "image6.png",
    "figure6_active_ibd_validation.png": "image5.png",
    "figure7_all_l2g_pathways.png": "image7.png",
    "figure8_external_experimental_evidence.png": "image8.png",
    "figure9_perturbseq_validation.png": "image9.png",
    "figure10_evidence_robustness.png": "image10.png",
    "figure11_benchmark_translation.png": "image11.png",
    "figure12_prkcb_mechanism.png": "image12.png",
    "figure13_clinical_translation.png": "image13.png",
}


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def fit_to_slot(source: Path, size: tuple[int, int]) -> Image.Image:
    image = Image.open(source).convert("RGB")
    contained = ImageOps.contain(image, size, method=Image.Resampling.LANCZOS)
    canvas = Image.new("RGB", size, "white")
    left = (size[0] - contained.width) // 2
    top = (size[1] - contained.height) // 2
    canvas.paste(contained, (left, top))
    return canvas


def prepare_images() -> dict[str, Path]:
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    outputs: dict[str, Path] = {}
    with zipfile.ZipFile(DOCX) as archive:
        for figure_name, media_name in FIGURE_TO_MEDIA.items():
            source = FIGURE_DIR / figure_name
            if not source.exists():
                raise FileNotFoundError(source)
            old = Image.open(io.BytesIO(archive.read(f"word/media/{media_name}")))
            optimized = fit_to_slot(source, old.size)
            output = WORK_DIR / media_name
            optimized.save(output, dpi=(300, 300), optimize=True)
            outputs[media_name] = output
    return outputs


def make_contact_sheet(outputs: dict[str, Path]) -> Path:
    figure_order = [
        (index, figure, media)
        for figure, media in FIGURE_TO_MEDIA.items()
        for index in range(1, 14)
        if media == f"image{index}.png"
    ]
    figure_order = sorted(figure_order)
    thumb_width = 380
    cells = []
    for index, source, media in figure_order:
        image = Image.open(outputs[media]).convert("RGB")
        height = round(image.height * thumb_width / image.width)
        cells.append(
            (
                index,
                image.resize((thumb_width, height), Image.Resampling.LANCZOS),
                height,
            )
        )
    columns = 3
    gutter = 24
    label_height = 30
    rows = (len(cells) + columns - 1) // columns
    row_heights = [
        max(cells[row * columns + column][2] for column in range(columns) if row * columns + column < len(cells))
        for row in range(rows)
    ]
    width = columns * thumb_width + (columns - 1) * gutter
    height = sum(row_heights) + rows * (gutter + label_height) + gutter
    sheet = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(sheet)
    font = ImageFont.load_default()
    y = gutter
    for row in range(rows):
        for column in range(columns):
            index = row * columns + column
            if index >= len(cells):
                break
            figure_number, image, image_height = cells[index]
            x = gutter + column * (thumb_width + gutter)
            draw.text(
                (x, y),
                f"Figure {figure_number}",
                fill="black",
                font=font,
            )
            sheet.paste(image, (x, y + label_height))
        y += label_height + row_heights[row] + gutter
    output = WORK_DIR / "1120_optimized_contact_sheet.png"
    sheet.save(output, dpi=(300, 300))
    return output


def embed_images(outputs: dict[str, Path]) -> dict:
    original_bytes = DOCX.read_bytes()
    if not BACKUP.exists():
        BACKUP.parent.mkdir(parents=True, exist_ok=True)
        BACKUP.write_bytes(original_bytes)
    temp = DOCX.with_suffix(".figure_optimization.tmp.docx")
    changed = []
    with zipfile.ZipFile(io.BytesIO(original_bytes), "r") as source_zip:
        with zipfile.ZipFile(
            temp,
            "w",
            compression=zipfile.ZIP_DEFLATED,
        ) as target_zip:
            for info in source_zip.infolist():
                data = source_zip.read(info.filename)
                media_name = info.filename.removeprefix("word/media/")
                if info.filename.startswith("word/media/") and media_name in outputs:
                    data = outputs[media_name].read_bytes()
                    changed.append(info.filename)
                target_zip.writestr(info, data)
    with zipfile.ZipFile(DOCX) as old_zip, zipfile.ZipFile(temp) as new_zip:
        old_names = set(old_zip.namelist())
        new_names = set(new_zip.namelist())
        document_unchanged = old_zip.read("word/document.xml") == new_zip.read(
            "word/document.xml"
        )
        nonmedia_unchanged = all(
            old_zip.read(name) == new_zip.read(name)
            for name in old_names
            if not name.startswith("word/media/")
        )
    os.replace(temp, DOCX)
    audit = {
        "document": str(DOCX),
        "backup": str(BACKUP),
        "figures_replaced": len(changed),
        "media_changed": changed,
        "document_xml_unchanged": document_unchanged,
        "nonmedia_entries_unchanged": nonmedia_unchanged,
        "all_checks_pass": (
            len(changed) == 13
            and document_unchanged
            and nonmedia_unchanged
            and old_names == new_names
        ),
        "optimized_sha256": sha256(DOCX.read_bytes()),
    }
    AUDIT.write_text(json.dumps(audit, indent=2), encoding="utf-8")
    return audit


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--embed", action="store_true")
    args = parser.parse_args()

    outputs = prepare_images()
    contact_sheet = make_contact_sheet(outputs)
    print(f"Prepared {len(outputs)} optimized figures.")
    print(contact_sheet)
    if args.embed:
        audit = embed_images(outputs)
        print(json.dumps(audit, indent=2))
        if not audit["all_checks_pass"]:
            raise SystemExit("1120 figure embedding audit failed.")


if __name__ == "__main__":
    main()
