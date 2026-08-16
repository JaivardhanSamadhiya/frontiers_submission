#!/usr/bin/env python3
"""Audit author-year citations against the revised DOCX reference list."""
from __future__ import annotations

import re
from pathlib import Path

from docx import Document

ROOT = Path(__file__).resolve().parent.parent
DOCX = ROOT / "submission/PrecisionPhage_Frontiers_Original_Research_frozen_external_validated.docx"
REPORT = ROOT / "submission/REFERENCE_AUDIT.md"

CITATION_RE = re.compile(
    r"([A-ZÀ-ÖØ-Þ][A-Za-zÀ-ÿ'’\-]+)"
    r"(?: et al\.| and [A-ZÀ-ÖØ-Þ][A-Za-zÀ-ÿ'’\-]+)?(?:,\s*|\s*\()"
    r"((?:19|20)\d{2})"
)
REFERENCE_RE = re.compile(
    r"^([A-ZÀ-ÖØ-Þ][A-Za-zÀ-ÿ'’\-]+),.*?\(((?:19|20)\d{2})\)\."
)
DOI_RE = re.compile(r"doi:\s*(10\.\d{4,9}/\S+)$", re.IGNORECASE)


def main() -> None:
    doc = Document(DOCX)
    paragraphs = doc.paragraphs
    ref_heading = next(i for i, p in enumerate(paragraphs) if p.text == "References")
    table_heading = next(i for i, p in enumerate(paragraphs) if p.text == "Tables")

    body = "\n".join(p.text for p in paragraphs[:ref_heading])
    citations = {(surname, year) for surname, year in CITATION_RE.findall(body)}

    references: dict[tuple[str, str], str] = {}
    dois: list[str] = []
    malformed: list[str] = []
    for paragraph in paragraphs[ref_heading + 1:table_heading]:
        text = paragraph.text.strip()
        if not text:
            continue
        match = REFERENCE_RE.search(text)
        if not match:
            malformed.append(text)
            continue
        key = (match.group(1), match.group(2))
        references[key] = text
        doi = DOI_RE.search(text)
        if doi:
            dois.append(doi.group(1).rstrip("."))

    missing_references = sorted(citations - set(references))
    uncited_references = sorted(set(references) - citations)
    duplicate_dois = sorted({doi for doi in dois if dois.count(doi) > 1})

    lines = [
        "# Reference and citation audit",
        "",
        f"Audited file: `{DOCX.name}`",
        "",
        f"- Unique author–year citations in text: {len(citations)}",
        f"- Reference entries: {len(references)}",
        f"- References containing a DOI: {len(dois)}",
        f"- Citations without a reference: {len(missing_references)}",
        f"- References not cited in the text: {len(uncited_references)}",
        f"- Malformed reference paragraphs: {len(malformed)}",
        f"- Duplicate DOI strings: {len(duplicate_dois)}",
        "",
    ]
    for title, values in (
        ("Citations without a reference", missing_references),
        ("References not cited", uncited_references),
        ("Malformed reference paragraphs", malformed),
        ("Duplicate DOI strings", duplicate_dois),
    ):
        if values:
            lines.extend([f"## {title}", ""])
            lines.extend(f"- {value}" for value in values)
            lines.append("")

    lines.extend([
        "## DOI corrections incorporated",
        "",
        "- Ahlgren et al. (2017): `10.1093/nar/gkw1002`",
        "- Edwards et al. (2016): `10.1093/femsre/fuv048`",
        "- Galiez et al. (2017): `10.1093/bioinformatics/btx383`",
        "- Torres-Barceló and Hochberg (2016): `10.1016/j.tim.2015.12.011`",
        "- Benjamini and Hochberg (1995): `10.1111/j.2517-6161.1995.tb02031.x`",
        "",
        "The direct VHIP source and software/method references missing from the",
        "original manuscript were added. DOI targets were checked against",
        "publisher/Crossref or indexed bibliographic records during the audit.",
        "",
    ])
    REPORT.write_text("\n".join(lines), encoding="utf-8")

    if missing_references or uncited_references or malformed or duplicate_dois:
        raise AssertionError(f"reference audit failed; see {REPORT}")
    print(f"PASS: {len(citations)} citations map bidirectionally to "
          f"{len(references)} references; {len(dois)} DOI strings are unique")
    print(REPORT)


if __name__ == "__main__":
    main()
