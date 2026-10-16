import sys
from pathlib import Path
import tempfile
import unittest
import numpy as np
import sourmash
from Bio.Seq import Seq
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from precisionphage.splits.canonical import sketch_fasta, distance_matrix


class CanonicalTests(unittest.TestCase):
    def test_reverse_complement_invariance(self):
        sequence=''.join(np.random.default_rng(13).choice(list('ACGT'),10000))
        a,b=sourmash.MinHash(2048,21),sourmash.MinHash(2048,21)
        a.add_sequence(sequence); b.add_sequence(str(Seq(sequence).reverse_complement()))
        self.assertEqual(a.hashes,b.hashes)
        self.assertEqual(distance_matrix([a,b])[0,1],0)

    def test_records_and_ambiguous_bases_do_not_join(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'tiny.fa'
            path.write_text('>a\nAACGTTT\n>b\nCCGGANNTAC\n')
            actual=sketch_fasta(path,num=100,k=5)
            expected=sourmash.MinHash(100,5)
            expected.add_sequence('AACGTTT'); expected.add_sequence('CCGGANNTAC',force=True)
            self.assertEqual(actual.hashes,expected.hashes)

    def test_no_valid_windows_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'empty.fa'; path.write_text('>a\nNNNNN\n')
            with self.assertRaises(ValueError): sketch_fasta(path)


if __name__=='__main__': unittest.main()
