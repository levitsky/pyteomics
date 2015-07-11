import unittest
import string
from itertools import count
import operator as op
import numpy as np
import pandas as pd
from pyteomics import auxiliary as aux

class QvalueTest(unittest.TestCase):
    def setUp(self):
        psms = list(zip(count(), string.ascii_uppercase + string.ascii_lowercase))
        np.random.shuffle(psms)
        self.psms = iter(psms)
        self.key = op.itemgetter(0)
        self.is_decoy = lambda x: x[1].islower()

    def _run_check(self, q, formula):
        self.assertTrue(np.allclose(q['q'][:26], 0))
        if formula == 2:
            self.assertTrue(np.allclose(q['q'][26:], 2 * np.arange(1., 27.) / (26 + np.arange(1, 27))))
        else:
            self.assertTrue(np.allclose(q['q'][26:], np.arange(1., 27.) / 26))
        self.assertTrue(np.allclose(q['is decoy'][:26], 0))
        self.assertTrue(np.allclose(q['is decoy'][26:], 1))
        self.assertTrue(np.allclose(q['score'], np.arange(52)))
        self.setUp()

    def test_qvalues(self):
        q = aux.qvalues(self.psms, key=self.key, is_decoy=self.is_decoy, remove_decoy=True)
        self.assertTrue(np.allclose(q['q'], 0))
        self.assertTrue(np.allclose(q['is decoy'], 0))
        self.assertTrue(np.allclose(q['score'], np.arange(26)))

    def test_qvalues_with_decoy(self):
        q = aux.qvalues(self.psms, key=self.key, is_decoy=self.is_decoy, remove_decoy=False)
        self._run_check(q, 2)
        q = aux.qvalues(self.psms, key=self.key, is_decoy=self.is_decoy, remove_decoy=False, formula=1)
        self._run_check(q, 1)

    def test_qvalues_full_output(self):
        q = aux.qvalues(self.psms, key=self.key, is_decoy=self.is_decoy, remove_decoy=False, full_output=True)
        self._run_check(q, 2)

    def test_qvalues_from_numpy(self):
        dtype = [('score', np.int8), ('label', np.str_, 1)]
        psms = np.array(list(self.psms), dtype=dtype)
        q = aux.qvalues(psms, key=self.key, is_decoy=self.is_decoy, remove_decoy=False, formula=1)
        self._run_check(q, 1)
        q = aux.qvalues(psms, key=self.key, is_decoy=self.is_decoy, remove_decoy=False, formula=1,
            full_output=True)
        self._run_check(q, 1)
        self.assertTrue(q['psm'].dtype == dtype)

    def test_qvalues_from_dataframe(self):
        dtype = [('score', np.int8), ('label', np.str_, 1)]
        psms = pd.DataFrame(np.array(list(self.psms), dtype=dtype))
        q = aux.qvalues(psms, key=self.key, is_decoy=self.is_decoy, remove_decoy=False, formula=1)
        self._run_check(q, 1)
        q = aux.qvalues(psms, key=self.key, is_decoy=self.is_decoy, remove_decoy=False, formula=1,
            full_output=True)
        self._run_check(q, 1)

    def test_qvalues_from_numpy_string_key(self):
        dtype = [('score', np.int8), ('label', np.str_, 1)]
        psms = np.array(list(self.psms), dtype=dtype)
        q = aux.qvalues(psms, key='score', is_decoy=self.is_decoy, remove_decoy=False, formula=1)
        self._run_check(q, 1)
        q = aux.qvalues(psms, key='score', is_decoy=self.is_decoy, remove_decoy=False, formula=1,
            full_output=True)
        self._run_check(q, 1)
        self.assertTrue(q['psm'].dtype == dtype)

    def test_qvalues_from_dataframe_string_key(self):
        dtype = [('score', np.int8), ('label', np.str_, 1)]
        psms = pd.DataFrame(np.array(list(self.psms), dtype=dtype))
        q = aux.qvalues(psms, key='score', is_decoy=self.is_decoy, remove_decoy=False, formula=1)
        self._run_check(q, 1)
        q = aux.qvalues(psms, key='score', is_decoy=self.is_decoy, remove_decoy=False, formula=1,
            full_output=True)
        self._run_check(q, 1)

    def test_qvalues_from_dataframe_string_key_and_is_decoy(self):
        dtype = [('score', np.int8), ('label', np.str_, 1)]
        psms = pd.DataFrame(np.array(list(self.psms), dtype=dtype))
        psms['is decoy'] = [self.is_decoy(row) for _, row in psms.iterrows()]
        q = aux.qvalues(psms, key='score', is_decoy='is decoy', remove_decoy=False, formula=1)
        self._run_check(q, 1)
        q = aux.qvalues(psms, key='score', is_decoy='is decoy', remove_decoy=False, formula=1,
            full_output=True)
        self._run_check(q, 1)

class FilterTest(unittest.TestCase):
    def setUp(self):
        self.psms = list(zip(count(), string.ascii_uppercase + string.ascii_lowercase))
        np.random.shuffle(self.psms)
        self.key = op.itemgetter(0)
        self.is_decoy = lambda x: x[1].islower()

    def test_filter(self):
        f11 = aux.filter(self.psms, key=self.key, is_decoy=self.is_decoy, fdr=0.5)
        f12 = aux.filter(self.psms, key=self.key, is_decoy=self.is_decoy, fdr=0.5,
            formula=2)
        f21 = aux.filter(self.psms, key=self.key, is_decoy=self.is_decoy, fdr=0.5,
            remove_decoy=False, formula=1)
        f22 = aux.filter(iter(self.psms), key=self.key, is_decoy=self.is_decoy, fdr=0.5,
            remove_decoy=False)
        self.assertEqual(f11.size, 26)
        self.assertEqual(f12.size, 26)
        self.assertEqual(len(f21), 39)
        self.assertEqual(f22.shape[0], 34)

    def test_filter_with(self):
        with aux.filter(self.psms, key=self.key, is_decoy=self.is_decoy, fdr=0.5,
            full_output=False) as f:
            f11 = list(f)
        with aux.filter(self.psms, key=self.key, is_decoy=self.is_decoy, fdr=0.5,
            formula=2, full_output=False) as f:
            f12 = list(f)
        with aux.filter(self.psms, key=self.key, is_decoy=self.is_decoy, fdr=0.5,
            remove_decoy=False, formula=1, full_output=False) as f:
            f21 = list(f)
        with aux.filter(iter(self.psms), key=self.key, is_decoy=self.is_decoy, fdr=0.5,
            remove_decoy=False, full_output=False) as f:
            f22 = list(f)
        self.assertEqual(len(f11), 26)
        self.assertEqual(len(f12), 26)
        self.assertEqual(len(f21), 39)
        self.assertEqual(len(f22), 34)

    def test_filter_array(self):
        dtype = [('score', np.int8), ('label', np.str_, 1)]
        psms = np.array(self.psms, dtype=dtype)
        f11 = aux.filter(psms, key=self.key, is_decoy=self.is_decoy, fdr=0.5)
        f12 = aux.filter(psms, key=self.key, is_decoy=self.is_decoy, fdr=0.5, formula=2)
        f21 = aux.filter(psms, key=self.key, is_decoy=self.is_decoy, fdr=0.5, remove_decoy=False, formula=1)
        f22 = aux.filter(psms, key=self.key, is_decoy=self.is_decoy, fdr=0.5, remove_decoy=False)
        self.assertEqual(f11.size, 26)
        self.assertEqual(f12.size, 26)
        self.assertEqual(len(f21), 39)
        self.assertEqual(f22.shape[0], 34)

    def test_filter_array_with(self):
        dtype = [('score', np.int8), ('label', np.str_, 1)]
        psms = np.array(self.psms, dtype=dtype)
        with aux.filter(psms, key=self.key, is_decoy=self.is_decoy, fdr=0.5,
            full_output=False) as f:
            f11 = list(f)
        with aux.filter(psms, key=self.key, is_decoy=self.is_decoy, fdr=0.5,
            formula=2, full_output=False) as f:
            f12 = list(f)
        with aux.filter(psms, key=self.key, is_decoy=self.is_decoy, fdr=0.5,
            remove_decoy=False, formula=1, full_output=False) as f:
            f21 = list(f)
        with aux.filter(psms, key=self.key, is_decoy=self.is_decoy, fdr=0.5,
            remove_decoy=False, full_output=False) as f:
            f22 = list(f)
        self.assertEqual(len(f11), 26)
        self.assertEqual(len(f12), 26)
        self.assertEqual(len(f21), 39)
        self.assertEqual(len(f22), 34)

    def test_filter_array_str_is_decoy(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('is decoy', np.bool)]
        psms = np.array([(s, l, self.is_decoy((s, l))) for s, l in self.psms], dtype=dtype)
        f11 = aux.filter(psms, key=self.key, is_decoy='is decoy', fdr=0.5)
        f12 = aux.filter(psms, key=self.key, is_decoy='is decoy', fdr=0.5, formula=2)
        f21 = aux.filter(psms, key=self.key, is_decoy='is decoy', fdr=0.5, remove_decoy=False, formula=1)
        f22 = aux.filter(psms, key=self.key, is_decoy='is decoy', fdr=0.5, remove_decoy=False)
        self.assertEqual(f11.size, 26)
        self.assertEqual(f12.size, 26)
        self.assertEqual(len(f21), 39)
        self.assertEqual(f22.shape[0], 34)

    def test_filter_dataframe(self):
        dtype = [('score', np.int8), ('label', np.str_, 1)]
        psms = pd.DataFrame(np.array(self.psms, dtype=dtype))
        f11 = aux.filter(psms, key=self.key, is_decoy=self.is_decoy, fdr=0.5)
        f12 = aux.filter(psms, key=self.key, is_decoy=self.is_decoy, fdr=0.5, formula=2)
        f21 = aux.filter(psms, key=self.key, is_decoy=self.is_decoy, fdr=0.5, remove_decoy=False, formula=1)
        f22 = aux.filter(psms, key=self.key, is_decoy=self.is_decoy, fdr=0.5, remove_decoy=False)
        self.assertEqual(f11.size, 26)
        self.assertEqual(f12.size, 26)
        self.assertEqual(len(f21), 39)
        self.assertEqual(f22.shape[0], 34)

if __name__ == '__main__':
    unittest.main()