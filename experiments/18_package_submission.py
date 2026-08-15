#!/usr/bin/env python3
"""Export primary figures as Frontiers-ready RGB 300-dpi TIFF files."""
from __future__ import annotations

from pathlib import Path

from PIL import Image

ROOT = Path(__file__).resolve().parent.parent
RD = ROOT / "data/results_v2"
OUT = ROOT / "submission/figures"
FIGURES = [
    ("figure_main.png", "Figure_1.tiff"),
    ("fig_feature_importance.png", "Figure_2.tiff"),
    ("fig_phist_compare.png", "Figure_3.tiff"),
    ("cocktail_coverage.png", "Figure_4.tiff"),
    ("temporal_dynamics.png", "Figure_5.tiff"),
]


def to_rgb(image: Image.Image) -> Image.Image:
    if image.mode == "RGB":
        return image
    if image.mode in {"RGBA", "LA"}:
        background = Image.new("RGBA", image.size, "white")
        background.alpha_composite(image.convert("RGBA"))
        return background.convert("RGB")
    return image.convert("RGB")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for source_name, output_name in FIGURES:
        source = RD / source_name
        if not source.is_file():
            raise FileNotFoundError(source)
        with Image.open(source) as image:
            rgb = to_rgb(image)
            output = OUT / output_name
            rgb.save(output, format="TIFF", compression="tiff_lzw",
                     dpi=(300, 300))
        with Image.open(output) as check:
            if check.mode != "RGB":
                raise AssertionError(f"{output_name}: expected RGB, got {check.mode}")
            dpi = check.info.get("dpi", (0, 0))
            if min(dpi) < 299:
                raise AssertionError(f"{output_name}: expected 300 dpi, got {dpi}")
            print(f"{output_name}: {check.size[0]}x{check.size[1]}, RGB, "
                  f"{float(dpi[0]):.0f} dpi")


if __name__ == "__main__":
    main()
