"""Recover versioned training records without claiming raw-byte reconstruction.

Uses the header-derived metadata included with Supplementary Data S1. No network
requests occur unless --download is explicitly supplied. Existing files are
verified, never overwritten. New FASTAs preserve original order and header text,
but standardize uppercase sequence/70-column/LF layout. Sequence identity and
historical raw-byte identity are separate checks. Python 3.11 standard library.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import textwrap
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[2]
PROVENANCE = ROOT / "data/results_v2/frontiers_revision_2026_10_14/genome_provenance"


def sha256(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def fasta_records(raw: bytes) -> list[tuple[str, bytes]]:
    records, header, chunks = [], None, []
    for line in raw.splitlines():
        if line.startswith(b">"):
            if header is not None:
                records.append((header, b"".join(chunks).upper()))
            header, chunks = line[1:].decode("utf-8").strip(), []
        elif line.strip():
            if header is None:
                raise ValueError("Sequence data before FASTA header")
            chunks.append(b"".join(line.split()))
    if header is not None:
        records.append((header, b"".join(chunks).upper()))
    if not records:
        raise ValueError("No FASTA records")
    return records


def verify_sequence(sequence: bytes, record: dict[str, str]) -> None:
    if len(sequence) != int(record["length"]):
        raise ValueError(f"Sequence length mismatch for {record['leading_versioned_record_accession']}")
    if sha256(sequence) != record["sequence_sha256_uppercase"]:
        raise ValueError(f"Sequence SHA-256 mismatch for {record['leading_versioned_record_accession']}")


def verify_file(raw: bytes, expected: list[dict[str, str]]) -> bool:
    parsed = fasta_records(raw)
    if len(parsed) != len(expected):
        raise ValueError("FASTA record-count mismatch")
    for (header, sequence), record in zip(parsed, expected):
        if header != record["header"]:
            raise ValueError("Existing FASTA header/order differs from inventory; file is not overwritten")
        verify_sequence(sequence, record)
    return sha256(raw) == expected[0]["frozen_file_sha256"]


class Fetcher:
    def __init__(self, delay: float, attempts: int, email: str | None):
        self.delay, self.attempts, self.email = max(delay, 0.4), attempts, email
        self.last_request = 0.0
        self.requests = 0

    def fetch(self, record: dict[str, str]) -> bytes:
        accession = record["leading_versioned_record_accession"]
        if not re.fullmatch(r"[A-Z][A-Z0-9_]*\d\.\d+", accession):
            raise ValueError(f"Missing/invalid versioned accession: {accession!r}")
        parameters = {"db": "nuccore", "id": accession, "rettype": "fasta",
                      "retmode": "text", "tool": "PrecisionPhage_provenance_recovery"}
        if self.email:
            parameters["email"] = self.email
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urlencode(parameters)
        max_bytes = 2 * int(record["length"]) + 1_000_000
        for attempt in range(self.attempts):
            time.sleep(max(0.0, self.delay - (time.monotonic() - self.last_request)))
            self.last_request = time.monotonic()
            self.requests += 1
            try:
                request = Request(url, headers={"User-Agent": "PrecisionPhage-record-recovery/1.0"})
                with urlopen(request, timeout=60) as response:
                    raw = response.read(max_bytes + 1)
                if len(raw) > max_bytes:
                    raise ValueError(f"Unexpectedly large response for {accession}")
                parsed = fasta_records(raw)
                if len(parsed) != 1:
                    raise ValueError(f"NCBI returned {len(parsed)} records for {accession}, expected one")
                header, sequence = parsed[0]
                if not re.search(r"(?<![A-Z0-9_])" + re.escape(accession) + r"(?![A-Z0-9_.])", header):
                    raise ValueError(f"NCBI response does not identify exact version {accession}")
                verify_sequence(sequence, record)
                return sequence
            except HTTPError as exc:
                if exc.code not in {429, 500, 502, 503, 504} or attempt + 1 >= self.attempts:
                    raise
                retry_after = exc.headers.get("Retry-After", "")
                pause = float(retry_after) if retry_after.isdigit() else 2 ** attempt
                time.sleep(min(30.0, max(1.0, pause)))
            except (URLError, TimeoutError):
                if attempt + 1 >= self.attempts:
                    raise
                time.sleep(min(30.0, 2 ** attempt))
        raise RuntimeError("Unreachable fetch retry state")


def load_inventory(path: Path) -> list[tuple[str, list[dict[str, str]]]]:
    grouped = defaultdict(list)
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            relative = PurePosixPath(row["fasta_path"])
            if relative.is_absolute() or ".." in relative.parts or ":" in row["fasta_path"]:
                raise ValueError(f"Unsafe relative path in inventory: {relative}")
            grouped[relative.as_posix()].append(row)
    for path, records in grouped.items():
        records.sort(key=lambda record: int(record["record_index"]))
        if [int(row["record_index"]) for row in records] != list(range(1, len(records) + 1)):
            raise ValueError(f"Nonconsecutive record order: {path}")
        if len({row["frozen_file_sha256"] for row in records}) != 1:
            raise ValueError(f"Inconsistent file checksums: {path}")
    return sorted(grouped.items(), key=lambda item: (item[1][0]["axis"], item[1][0]["entity"]))


def reconstruct(records: list[dict[str, str]], fetcher: Fetcher) -> bytes:
    parts = []
    for record in records:
        sequence = fetcher.fetch(record).decode("ascii")
        parts.append(">" + record["header"] + "\n" + "\n".join(textwrap.wrap(sequence, 70)) + "\n")
    result = "".join(parts).encode("utf-8")
    verify_file(result, records)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", type=Path, default=PROVENANCE / "training_contig_accessions.csv")
    parser.add_argument("--target-root", type=Path, default=ROOT / "recovered_training_fastas")
    parser.add_argument("--axis", choices=["all", "phage", "host"], default="all")
    parser.add_argument("--limit", type=int, help="Maximum files, for a small test")
    parser.add_argument("--download", action="store_true", help="Explicitly allow NCBI requests for missing files")
    parser.add_argument("--delay", type=float, default=0.6, help="Seconds between requests (minimum 0.4)")
    parser.add_argument("--attempts", type=int, default=4)
    parser.add_argument("--email", help="Optional contact email sent to NCBI; no API key needed")
    parser.add_argument("--report", type=Path, help="Optional NEW JSON report; existing reports are never overwritten")
    args = parser.parse_args()
    if (args.limit is not None and args.limit < 1) or args.attempts < 1:
        parser.error("--limit and --attempts must be positive")
    root = args.target_root.resolve()
    selected = [(path, rows) for path, rows in load_inventory(args.inventory)
                if args.axis == "all" or rows[0]["axis"] == args.axis]
    if args.limit:
        selected = selected[:args.limit]
    if not selected:
        raise ValueError("No inventory records selected")
    fetcher, results = Fetcher(args.delay, args.attempts, args.email), []
    for relative, records in selected:
        target = (root / Path(relative)).resolve()
        if not target.is_relative_to(root):
            raise ValueError(f"Resolved output escapes target root: {target}")
        if target.exists():
            raw = target.read_bytes()
            raw_exact = verify_file(raw, records)
            status = "existing-sequence-verified-no-overwrite"
        elif not args.download:
            results.append({"path": relative, "status": "missing-network-not-authorized"})
            continue
        else:
            # Verify every record before a new output file is created.
            raw = reconstruct(records, fetcher)
            target.parent.mkdir(parents=True, exist_ok=True)
            with target.open("xb") as handle:
                handle.write(raw)
            raw_exact = verify_file(target.read_bytes(), records)
            status = "downloaded-sequence-verified-standardized-layout"
        results.append({"path": relative, "status": status, "records": len(records),
                        "file_sha256": sha256(raw), "historical_file_sha256": records[0]["frozen_file_sha256"],
                        "raw_bytes_identical_to_historical_file": raw_exact,
                        "all_record_sequences_length_and_sha256_verified": True})
        print(f"Verified {relative}: sequence match; historical raw-byte match={raw_exact}", flush=True)
    summary = {"created_utc": datetime.now(timezone.utc).isoformat(),
               "inventory_sha256": sha256(args.inventory.read_bytes()), "network_allowed": args.download,
               "http_requests": fetcher.requests, "target_root": str(root),
               "status_counts": dict(Counter(row["status"] for row in results)),
               "historical_raw_file_matches": sum(row.get("raw_bytes_identical_to_historical_file", False) for row in results),
               "limitations": ["Exact accession.version and sequence length/hash are checked; changed or unavailable records fail closed.",
                               "New files preserve original headers/order but use uppercase/70-column/LF layout.",
                               "Sequence identity is NOT proof of raw-byte recovery. The historical raw-file SHA-256 must separately match.",
                               "No existing input is overwritten; mismatching existing sequences or headers raise an error."],
               "files": results}
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        with args.report.open("x", encoding="utf-8", newline="\n") as handle:
            json.dump(summary, handle, indent=2)
            handle.write("\n")
    print(json.dumps({key: value for key, value in summary.items() if key != "files"}, indent=2))


if __name__ == "__main__":
    main()
