import importlib.util
from pathlib import Path, PurePosixPath, PureWindowsPath
import unittest
import numpy as np
import pandas as pd
import json
from sklearn.metrics import roc_auc_score

spec=importlib.util.spec_from_file_location('doublets',Path(__file__).with_name('doublet_validation.py'))
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

def recorded_relative_path(value):
    """Interpret preserved Windows/POSIX provenance paths on either platform."""
    normalized=value.replace('\\','/')
    path=PurePosixPath(normalized)
    if path.is_absolute() or PureWindowsPath(value).drive or '..' in path.parts:
        raise ValueError('Expected a repository-relative provenance path')
    return Path(*path.parts)

class TestMeasuredMenus(unittest.TestCase):
    def test_recorded_path_separators_are_portable(self):
        expected=Path('scripts')/'study'/'protocol.md'
        self.assertEqual(recorded_relative_path('scripts/study/protocol.md'),expected)
        self.assertEqual(recorded_relative_path(r'scripts\study\protocol.md'),expected)
    def test_recorded_path_cannot_escape_repository(self):
        for value in ['/tmp/input.csv',r'C:\study\input.csv',r'..\input.csv','../input.csv',r'\\server\share\input.csv']:
            with self.subTest(value=value),self.assertRaises(ValueError):
                recorded_relative_path(value)
    def test_randomized_ties(self):
        value,ties=m.choose_expected(np.array([[.8,.8,.1],[.2,.2,.2]]),np.array([[0,1,1],[0,1,1]]))
        np.testing.assert_allclose(value,[.5,2/3])
        np.testing.assert_array_equal(ties.sum(axis=1),[2,3])
    def test_ties_do_not_depend_on_row_order(self):
        s=np.array([[.9,.9,.3]]);y=np.array([[1,0,0]])
        self.assertEqual(m.choose_expected(s,y)[0][0],m.choose_expected(s[:,::-1],y[:,::-1])[0][0])
    def test_perfect_and_reverse(self):
        self.assertEqual(m.pair_concordance([1,0],[.9,.1]),(1.,1.))
        self.assertEqual(m.pair_concordance([1,0],[.1,.9]),(0.,1.))
    def test_constant_and_singleclass(self):
        self.assertEqual(m.pair_concordance([1,0,0],[.5,.5,.5]),(1.,2.))
        self.assertEqual(m.pair_concordance([1,1],[.1,.9]),(0.,0.))
        self.assertTrue(np.isnan(m.conditional_metrics(np.array([1,1]),np.array([.1,.9]),np.array(['x','x']),np.ones(2))).all())
    def test_weighted_concordance_sklearn(self):
        rng=np.random.default_rng(11)
        for _ in range(50):
            y=rng.integers(0,2,30);s=rng.integers(0,6,30)/5;w=rng.integers(0,5,30)
            a,b=m.pair_concordance(y,s,w)
            if b:self.assertAlmostEqual(a/b,roc_auc_score(y,s,sample_weight=w),places=12)
    def test_composition_control(self):
        v=m.conditional_metrics(np.array([1,0,1,0]),np.array([.1,.1,.8,.8]),np.array(['a','a','b','b']),np.ones(4))
        np.testing.assert_allclose(v,[.5,.5])
    def test_cluster_multiplicity(self):
        w=m.bootstrap_weights(np.array([0,0,1,2]),np.array(['a','a','b','b']))
        np.testing.assert_array_equal(w[:,0],w[:,1])
        np.testing.assert_array_equal(w[:,0]+w[:,2]+w[:,3],np.full(m.REPS,3))
    def test_scores_never_use_outcome(self):
        x=np.array([[[.2,.7],[.8,.3]]])
        np.testing.assert_allclose(m.aggregate(x,'product'),[[.76,.86]])
        with self.assertRaises(ValueError):m.aggregate(np.array([[np.nan,.3]]),'product')
    def test_saved_decisions_reconcile_to_source(self):
        if not (m.OUT/'summary.json').exists(): self.skipTest('Run analysis first')
        rows=pd.read_csv(m.OUT/'doublet_option_scores.csv')
        host=pd.read_csv(m.OUT/'doublet_host_decisions.csv').set_index('bacteria')
        source=pd.read_csv(m.INPUTS['outcomes'],sep=';').replace({'bacteria':{'BCH857':'BCH857-72'}}).set_index('bacteria')
        for (method,name),g in rows.groupby(['method','bacteria']):
            self.assertEqual(len(g),3)
            self.assertAlmostEqual(g.selection_probability.sum(),1.)
            self.assertAlmostEqual(float((g.active*g.selection_probability).sum()),host.loc[name,method])
            for row in g.itertuples():
                self.assertEqual(row.assay_value,source.loc[name,row.option])
                self.assertEqual(row.active,int(source.loc[name,row.option]>0))
        report=json.loads((m.OUT/'summary.json').read_text())
        self.assertEqual(report['informative_triplet_compositions'],3)
        self.assertEqual(report['doublet_positive'],224)
        for data in report['inputs'].values():
            self.assertEqual(data['sha256'],m.sha(m.ROOT/recorded_relative_path(data['path'])))

if __name__=='__main__':unittest.main()
