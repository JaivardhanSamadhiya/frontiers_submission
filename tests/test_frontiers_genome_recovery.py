"""Offline safety/identity checks for versioned training-FASTA recovery."""
from __future__ import annotations

import csv
import importlib.util
import io
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

MODULE = Path(__file__).resolve().parents[1] / "scripts/frontiers_revision/recover_training_fastas.py"
SPEC = importlib.util.spec_from_file_location("frontiers_genome_recovery", MODULE)
recovery = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(recovery)


class RecoveryChecks(unittest.TestCase):
    def setUp(self):
        self.sequence = b"ACGTN" * 23
        self.raw = b">NC_000866.4 preserved original header\r\n" + self.sequence + b"\r\n\r\n"
        self.record = {"axis": "phage", "entity": "example", "fasta_path": "original/example.fasta",
                       "record_index": "1", "header": "NC_000866.4 preserved original header",
                       "leading_versioned_record_accession": "NC_000866.4", "length": str(len(self.sequence)),
                       "sequence_sha256_uppercase": recovery.sha256(self.sequence),
                       "frozen_file_sha256": recovery.sha256(self.raw)}

    def test_original_bytes_and_normalized_sequence_are_distinct(self):
        self.assertTrue(recovery.verify_file(self.raw, [self.record]))
        standard = self.raw.replace(b"\r\n", b"\n")
        self.assertFalse(recovery.verify_file(standard, [self.record]))

    def test_changed_base_or_header_fails_closed(self):
        with self.assertRaises(ValueError):
            recovery.verify_file(self.raw.replace(b"ACGTN", b"TCGTN", 1), [self.record])
        with self.assertRaises(ValueError):
            recovery.verify_file(self.raw.replace(b"original header", b"new header"), [self.record])

    def test_versioned_fetch_checks_server_header_and_hash_offline(self):
        server = b">NC_000866.4 NCBI current description\n" + self.sequence.lower() + b"\n"
        with patch.object(recovery, "urlopen", return_value=io.BytesIO(server)), patch.object(recovery.time, "sleep"):
            fetcher = recovery.Fetcher(0.6, 1, None)
            self.assertEqual(fetcher.fetch(self.record), self.sequence)
            self.assertEqual(fetcher.requests, 1)
        bad_version = server.replace(b"NC_000866.4", b"NC_000866.5")
        with patch.object(recovery, "urlopen", return_value=io.BytesIO(bad_version)), patch.object(recovery.time, "sleep"):
            with self.assertRaises(ValueError):
                recovery.Fetcher(0.6, 1, None).fetch(self.record)

    def test_reconstruction_keeps_original_order_and_headers(self):
        second = dict(self.record, record_index="2", header="NZ_CP187563.1 second ordered header",
                      leading_versioned_record_accession="NZ_CP187563.1", length="3",
                      sequence_sha256_uppercase=recovery.sha256(b"GCA"))
        class CachedFetcher:
            def fetch(inner, record):
                return self.sequence if record["record_index"] == "1" else b"GCA"
        standardized = recovery.reconstruct([self.record, second], CachedFetcher())
        parsed = recovery.fasta_records(standardized)
        self.assertEqual([p[0] for p in parsed], [self.record["header"], second["header"]])
        self.assertEqual([p[1] for p in parsed], [self.sequence, b"GCA"])
        self.assertNotIn(b"\r", standardized)
        self.assertLessEqual(max(map(len, standardized.splitlines()[1:])), 70)

    def test_inventory_rejects_path_traversal_and_nonconsecutive_order(self):
        with tempfile.TemporaryDirectory(prefix="precisionphage-recovery-test-") as directory:
            path = Path(directory) / "inventory.csv"
            for bad in (dict(self.record, fasta_path="../bad.fasta"), dict(self.record, record_index="2")):
                with path.open("w", encoding="utf-8", newline="") as handle:
                    writer = csv.DictWriter(handle, fieldnames=list(bad))
                    writer.writeheader()
                    writer.writerow(bad)
                with self.assertRaises(ValueError):
                    recovery.load_inventory(path)


if __name__ == "__main__":
    unittest.main()
