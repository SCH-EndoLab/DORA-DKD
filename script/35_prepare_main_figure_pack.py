from __future__ import annotations

import os
import subprocess
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[3]
EXPORT = ROOT / "figure" / "export"
PNG_SRC = ROOT / "figure" / "final_main_png_src"
FINAL = ROOT / "figure" / "final_main"

PNG_SRC.mkdir(parents=True, exist_ok=True)
FINAL.mkdir(parents=True, exist_ok=True)


def font(size: int, bold: bool = False):
    candidates = []
    if bold:
        candidates.extend(
            [
                "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
                "/System/Library/Fonts/Supplemental/Helvetica.ttc",
            ]
        )
    else:
        candidates.extend(
            [
                "/System/Library/Fonts/Supplemental/Arial.ttf",
                "/System/Library/Fonts/Supplemental/Helvetica.ttc",
            ]
        )
    for candidate in candidates:
        if os.path.exists(candidate):
            return ImageFont.truetype(candidate, size=size)
    return ImageFont.load_default()


FONT_TITLE = font(50, bold=True)
FONT_SUBTITLE = font(28, bold=True)
FONT_BODY = font(24)
FONT_LABEL = font(34, bold=True)


def wrap_text(text: str, max_chars: int) -> list[str]:
    words = text.split()
    lines = []
    current = []
    for word in words:
        tentative = " ".join(current + [word])
        if len(tentative) <= max_chars:
            current.append(word)
        else:
            if current:
                lines.append(" ".join(current))
            current = [word]
    if current:
        lines.append(" ".join(current))
    return lines


def ql_png(pdf_path: Path) -> Path:
    out_path = PNG_SRC / f"{pdf_path.name}.png"
    if out_path.exists():
        return out_path
    subprocess.run(
        [
            "/usr/bin/qlmanage",
            "-t",
            "-s",
            "2400",
            "-o",
            str(PNG_SRC),
            str(pdf_path),
        ],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    generated = PNG_SRC / f"{pdf_path.name}.png"
    if not generated.exists():
        raise FileNotFoundError(f"Failed to render preview for {pdf_path}")
    return generated


def load_image(path: Path) -> Image.Image:
    if path.suffix.lower() == ".pdf":
        path = ql_png(path)
    img = Image.open(path)
    if img.mode != "RGB":
        img = img.convert("RGB")
    return img


def fit_image(img: Image.Image, target_w: int, target_h: int, pad: int = 12) -> Image.Image:
    working_w = max(20, target_w - 2 * pad)
    working_h = max(20, target_h - 2 * pad)
    scale = min(working_w / img.width, working_h / img.height)
    new_size = (max(1, int(img.width * scale)), max(1, int(img.height * scale)))
    resized = img.resize(new_size, Image.Resampling.LANCZOS)
    canvas = Image.new("RGB", (target_w, target_h), "white")
    x = (target_w - resized.width) // 2
    y = (target_h - resized.height) // 2
    canvas.paste(resized, (x, y))
    return canvas


def annotate_panel(draw: ImageDraw.ImageDraw, x: int, y: int, label: str) -> None:
    draw.rounded_rectangle((x, y, x + 52, y + 44), radius=8, fill="#f0f0f0", outline="#707070", width=2)
    draw.text((x + 14, y + 5), label, font=FONT_LABEL, fill="#1a1a1a")


def make_composite(title: str, panels: list[dict], out_stem: str, canvas_size: tuple[int, int] = (3200, 2200)) -> None:
    canvas = Image.new("RGB", canvas_size, "white")
    draw = ImageDraw.Draw(canvas)
    draw.text((70, 45), title, font=FONT_TITLE, fill="#111111")

    for panel in panels:
        panel_img = load_image(panel["path"])
        slot = fit_image(panel_img, panel["w"], panel["h"])
        canvas.paste(slot, (panel["x"], panel["y"]))
        annotate_panel(draw, panel["x"] + 12, panel["y"] + 12, panel["label"])

    png_path = FINAL / f"{out_stem}.png"
    tiff_path = FINAL / f"{out_stem}.tiff"
    canvas.save(png_path, format="PNG", optimize=True)
    canvas.save(tiff_path, format="TIFF", compression="tiff_lzw")


def make_study_design() -> None:
    width, height = 3200, 2200
    img = Image.new("RGB", (width, height), "#fcfcfa")
    draw = ImageDraw.Draw(img)

    draw.text((80, 50), "Figure 1. Multi-omics study design", font=FONT_TITLE, fill="#111111")
    draw.text((80, 115), "Diabetes-centered DKD project structure for direct manuscript insertion", font=FONT_SUBTITLE, fill="#35505b")

    boxes = [
        (120, 250, 1000, 530, "#e7f2ef", "Disease frame", [
            "Type 2 diabetes",
            "DKD progression",
            "Renal injury and failure risk",
        ]),
        (1120, 250, 3080, 530, "#eef3f8", "Data layers", [
            "Bulk kidney transcriptomes: GSE30528, GSE96804, GSE104954",
            "Proximal tubule methylation: GSE233758",
            "Kidney and urine single-cell data: GSE209781, GSE279086, GSE266146",
            "Blood diabetes-axis support cohort: GSE189007",
        ]),
        (120, 680, 960, 1140, "#fff3e6", "Bulk discovery and mechanism layer", [
            "Gene-level processing and limma",
            "Strong glomerular consensus construction",
            "Curated oxidative stress-associated mechanism scoring",
            "Direction-stratified functional enrichment",
        ]),
        (1030, 680, 1880, 1140, "#f7ecff", "Epigenetic support", [
            "Differential methylation in proximal tubules",
            "Mechanism-gene and candidate-gene CpG support",
        ]),
        (1950, 680, 3080, 1140, "#edf7ea", "Cell-state localization", [
            "Coarse cell-type inference",
            "PT and PT_INJURED focus",
            "Urine-side injury support",
        ]),
        (120, 1300, 1450, 1820, "#f4f6fa", "Candidate integration", [
            "Bulk priority",
            "Methylation support",
            "Tubular localization",
            "External kidney support",
            "Urine detectability",
        ]),
        (1560, 1300, 3080, 1820, "#fef0f0", "Main manuscript outputs", [
            "Loss of intracellular stress-handling and mitochondrial homeostasis",
            "Gain of ECM and immune remodeling",
            "Proximal tubular injury axis",
            "Main candidates: DCXR, ETFB, MACROD1, ACOT7",
        ]),
    ]

    for x1, y1, x2, y2, color, title, bullets in boxes:
        draw.rounded_rectangle((x1, y1, x2, y2), radius=28, fill=color, outline="#60717b", width=4)
        draw.text((x1 + 30, y1 + 24), title, font=FONT_SUBTITLE, fill="#162228")
        cursor_y = y1 + 90
        for bullet in bullets:
            wrapped = wrap_text(bullet, 42 if (x2 - x1) < 1000 else 58)
            for line_index, line in enumerate(wrapped):
                prefix = "- " if line_index == 0 else "  "
                draw.text((x1 + 38, cursor_y), prefix + line, font=FONT_BODY, fill="#26363d")
                cursor_y += 34
            cursor_y += 16

    arrows = [
        ((560, 530), (560, 680)),
        ((2100, 530), (2100, 680)),
        ((540, 1140), (540, 1300)),
        ((1455, 1560), (1560, 1560)),
        ((2450, 1140), (2450, 1300)),
    ]
    for (x1, y1), (x2, y2) in arrows:
        draw.line((x1, y1, x2, y2), fill="#7b8a91", width=10)
        draw.polygon([(x2, y2), (x2 - 20, y2 - 35), (x2 + 20, y2 - 35)], fill="#7b8a91")

    out_png = FINAL / "Figure_1_Study_Design.png"
    out_tiff = FINAL / "Figure_1_Study_Design.tiff"
    img.save(out_png, format="PNG", optimize=True)
    img.save(out_tiff, format="TIFF", compression="tiff_lzw")


def main() -> None:
    make_study_design()

    make_composite(
        title="Figure 2. Bulk kidney oxidative stress-associated remodeling in DKD",
        out_stem="Figure_2_Bulk_Mechanism",
        panels=[
            {"label": "A", "path": EXPORT / "mechanism" / "bulk_mechanism_scores_boxplot.pdf", "x": 70, "y": 170, "w": 1540, "h": 1900},
            {"label": "B", "path": EXPORT / "mechanism" / "bulk_oxeiptosis_core_gene_expression.pdf", "x": 1620, "y": 170, "w": 1510, "h": 1900},
        ],
    )

    make_composite(
        title="Figure 3. Functional enrichment of direction-stratified glomerular consensus genes",
        out_stem="Figure_3_Functional_Enrichment",
        panels=[
            {"label": "A", "path": EXPORT / "functional" / "glomerular_go_enrichment_dotplot.pdf", "x": 70, "y": 170, "w": 1540, "h": 1900},
            {"label": "B", "path": EXPORT / "functional" / "glomerular_pathway_enrichment_dotplot.pdf", "x": 1620, "y": 170, "w": 1510, "h": 1900},
        ],
    )

    make_composite(
        title="Figure 4. Proximal tubule methylation support in DKD",
        out_stem="Figure_4_Methylation",
        panels=[
            {"label": "A", "path": EXPORT / "mechanism" / "GSE233758_methylation_volcano.pdf", "x": 70, "y": 170, "w": 1020, "h": 1900},
            {"label": "B", "path": EXPORT / "mechanism" / "GSE233758_mechanism_gene_methylation_summary.pdf", "x": 1100, "y": 170, "w": 1020, "h": 1900},
            {"label": "C", "path": EXPORT / "mechanism" / "GSE233758_top_candidate_methylation_summary.pdf", "x": 2130, "y": 170, "w": 1000, "h": 1900},
        ],
    )

    make_composite(
        title="Figure 5. Single-cell localization of the DKD tubular injury axis",
        out_stem="Figure_5_Single_Cell_Localization",
        panels=[
            {"label": "A", "path": EXPORT / "single_cell" / "single_cell_typing_sample_composition.pdf", "x": 70, "y": 170, "w": 1500, "h": 890},
            {"label": "B", "path": EXPORT / "single_cell" / "single_cell_typing_group_fraction_heatmap.pdf", "x": 1600, "y": 170, "w": 1530, "h": 890},
            {"label": "C", "path": EXPORT / "single_cell" / "single_cell_tubular_axis_fraction_boxplot.pdf", "x": 70, "y": 1080, "w": 1020, "h": 990},
            {"label": "D", "path": EXPORT / "single_cell" / "single_cell_tubular_axis_score_boxplot.pdf", "x": 1100, "y": 1080, "w": 1020, "h": 990},
            {"label": "E", "path": EXPORT / "single_cell" / "single_cell_typing_group_oxeiptosis_heatmap.pdf", "x": 2130, "y": 1080, "w": 1000, "h": 990},
        ],
    )

    make_composite(
        title="Figure 6. Candidate localization and integrative prioritization",
        out_stem="Figure_6_Candidate_Localization",
        panels=[
            {"label": "A", "path": EXPORT / "single_cell" / "single_cell_candidate_gene_expression_heatmap.pdf", "x": 70, "y": 170, "w": 1040, "h": 1900},
            {"label": "B", "path": EXPORT / "single_cell" / "single_cell_candidate_gene_detection_heatmap.pdf", "x": 1120, "y": 170, "w": 1040, "h": 1900},
            {"label": "C", "path": EXPORT / "single_cell" / "single_cell_celltype_score_focus_dotplot.pdf", "x": 2170, "y": 170, "w": 960, "h": 930},
            {"label": "D", "path": EXPORT / "mechanism" / "multiomics_candidate_ppi_network.png", "x": 2170, "y": 1130, "w": 960, "h": 940},
        ],
    )

    make_composite(
        title="Figure 7. Whole-blood diabetes-axis contextual support",
        out_stem="Figure_7_Blood_Diabetes_Axis",
        panels=[
            {"label": "A", "path": EXPORT / "mechanism" / "GSE189007_diabetes_axis_mechanism_scores.pdf", "x": 70, "y": 170, "w": 1540, "h": 1900},
            {"label": "B", "path": EXPORT / "mechanism" / "GSE189007_diabetes_axis_top_candidates.pdf", "x": 1620, "y": 170, "w": 1510, "h": 1900},
        ],
    )

    manifest = FINAL / "README_figure_pack.md"
    manifest.write_text(
        "\n".join(
            [
                "# Main Figure Pack",
                "",
                "These files are prepared for direct insertion into the Word manuscript.",
                "",
                "- Figure 1: `Figure_1_Study_Design.png`",
                "- Figure 2: `Figure_2_Bulk_Mechanism.png`",
                "- Figure 3: `Figure_3_Functional_Enrichment.png`",
                "- Figure 4: `Figure_4_Methylation.png`",
                "- Figure 5: `Figure_5_Single_Cell_Localization.png`",
                "- Figure 6: `Figure_6_Candidate_Localization.png`",
                "- Figure 7: `Figure_7_Blood_Diabetes_Axis.png`",
                "",
                "Each figure is also exported as a LZW-compressed TIFF with the same stem.",
            ]
        ),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
