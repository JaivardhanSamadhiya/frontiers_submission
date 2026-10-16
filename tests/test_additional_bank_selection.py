import importlib.util
import io
from pathlib import Path
import unittest
import numpy as np
from Bio import Phylo
spec=importlib.util.spec_from_file_location('bank',Path(__file__).resolve().parents[1]/'experiments/41_additional_bank_selection.py')
bank=importlib.util.module_from_spec(spec); spec.loader.exec_module(bank)


class BankTests(unittest.TestCase):
    def test_only_assembly_versions_are_aliased(self):
        self.assertEqual(bank.accession_alias('GCA_123456.1'),'GCA_123456')
        self.assertEqual(bank.accession_alias('strain.1'),'strain.1')

    def test_duplicate_accession_alias_rejected(self):
        tree=Phylo.read(io.StringIO('(GCA_123456.1:0.1,GCA_123456.2:0.2);'),'newick')
        with self.assertRaises(ValueError): bank.tree_distances(tree,['GCA_123456'])

    def test_exact_patristic_distance(self):
        tree=Phylo.read(io.StringIO('((A:0.1,B:0.2):0.3,C:0.4);'),'newick')
        distance,_=bank.tree_distances(tree,['A','B','C'])
        np.testing.assert_allclose(distance,[[0,.3,.8],[.3,0,.9],[.8,.9,0]])

    def test_greedy_baseline_covers_complementary_hosts(self):
        y=np.array([[1,1,0,0],[1,1,0,0],[0,0,1,0],[0,0,0,1]])
        selected=bank.greedy_coverage(y,np.ones_like(y,dtype=bool))
        self.assertEqual(selected.tolist(),[0,2,3])


if __name__=='__main__': unittest.main()
