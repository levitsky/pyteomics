import unittest
import string
from itertools import count
import operator as op
import numpy as np
import pandas as pd
from os import path
import pyteomics
pyteomics.__path__ = [path.abspath(path.join(path.dirname(__file__), path.pardir, 'pyteomics'))]
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

    def test_qvalues_pep_exceptions(self):
        self.assertRaises(aux.PyteomicsError, aux.qvalues,
            self.psms, pep='score', is_decoy=self.is_decoy)
        self.assertRaises(aux.PyteomicsError, aux.qvalues,
            self.psms, pep='score', remove_decoy=False)
        self.assertRaises(aux.PyteomicsError, aux.qvalues,
            self.psms, pep='score', correction=0)

class FilterTest(unittest.TestCase):
    def setUp(self):
        self.psms = list(zip(count(), string.ascii_uppercase + string.ascii_lowercase))
        np.random.shuffle(self.psms)
        self.key = op.itemgetter(0)
        self.is_decoy = lambda x: x[1].islower()

    def _run_check(self, *args, **kwargs):
        key = kwargs.get('key', self.key)
        is_decoy = kwargs.get('is_decoy', self.is_decoy)
        f11 = aux.filter(*args, key=key, is_decoy=is_decoy, fdr=0.5)
        f12 = aux.filter(*args, key=key, is_decoy=is_decoy, fdr=0.5, formula=2)
        f21 = aux.filter(*args, key=key, is_decoy=is_decoy, fdr=0.5, remove_decoy=False, formula=1)
        f22 = aux.filter(*args, key=key, is_decoy=is_decoy, fdr=0.5, remove_decoy=False)

        self.assertEqual(f11.shape[0], 26)
        self.assertEqual(f12.shape[0], 26)
        self.assertEqual(f21.shape[0], 39)
        self.assertEqual(f22.shape[0], 34)

        with aux.filter(*args, key=key, is_decoy=is_decoy, fdr=0.5, full_output=False) as f:
            f11 = list(f)
        with aux.filter(*args, key=key, is_decoy=is_decoy, fdr=0.5, formula=2, full_output=False) as f:
            f12 = list(f)
        with aux.filter(*args, key=key, is_decoy=is_decoy, fdr=0.5, remove_decoy=False, formula=1, full_output=False) as f:
            f21 = list(f)
        with aux.filter(*args, key=key, is_decoy=is_decoy, fdr=0.5, remove_decoy=False, full_output=False) as f:
            f22 = list(f)

        self.assertEqual(len(f11), 26)
        self.assertEqual(len(f12), 26)
        self.assertEqual(len(f21), 39)
        self.assertEqual(len(f22), 34)

    def test_filter(self):
        self._run_check(self.psms)

    def test_filter_chain(self):
        f11 = aux.filter.chain(self.psms, self.psms, key=self.key, is_decoy=self.is_decoy, fdr=0.5)
        self.assertEqual(f11.shape[0], 52)

    def test_filter_chain_with(self):
        with aux.filter.chain(self.psms, self.psms, key=self.key, is_decoy=self.is_decoy,
            fdr=0.5, full_output=False) as f:
            f11 = list(f)
        self.assertEqual(len(f11), 52)

    def test_filter_chain_arr_str_key(self):
        dtype = [('score', np.int8), ('label', np.str_, 1)]
        psms = np.array(self.psms, dtype=dtype)
        f11 = aux.filter.chain(psms, psms, key='score', is_decoy=self.is_decoy, fdr=0.5)
        self.assertEqual(f11.shape[0], 52)

    def test_filter_chain_arr_str_key_with(self):
        dtype = [('score', np.int8), ('label', np.str_, 1)]
        psms = np.array(self.psms, dtype=dtype)
        with aux.filter.chain(psms, psms, key='score', is_decoy=self.is_decoy, fdr=0.5, full_output=False) as f:
            f11 = list(f)
        self.assertEqual(len(f11), 52)

    def test_filter_chain_from_iterable(self):
        f11 = aux.filter.chain.from_iterable([self.psms, self.psms], key=self.key, is_decoy=self.is_decoy, fdr=0.5)
        self.assertEqual(f11.shape[0], 52)

    def test_filter_chain_from_iterable_with(self):
        with aux.filter.chain.from_iterable([self.psms, self.psms],
            key=self.key, is_decoy=self.is_decoy, fdr=0.5, full_output=False) as f:
            f11 = list(f)
        self.assertEqual(len(f11), 52)

    def test_filter_array(self):
        dtype = [('score', np.int8), ('label', np.str_, 1)]
        psms = np.array(self.psms, dtype=dtype)
        self._run_check(psms)

    def test_filter_array_str_is_decoy(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('is decoy', np.bool)]
        psms = np.array([(s, l, self.is_decoy((s, l))) for s, l in self.psms], dtype=dtype)
        self._run_check(psms, is_decoy='is decoy')

    def test_filter_array_str_is_decoy_str_key(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('is decoy', np.bool)]
        psms = np.array([(s, l, self.is_decoy((s, l))) for s, l in self.psms], dtype=dtype)
        self._run_check(psms, is_decoy='is decoy', key='score')

    def test_filter_array_list_key(self):
        dtype = [('score', np.int8), ('label', np.str_, 1)]
        psms = np.array(self.psms, dtype=dtype)
        key = [self.key(psm) for psm in psms]
        self._run_check(psms, key=key)

    def test_filter_array_gen_key(self):
        dtype = [('score', np.int8), ('label', np.str_, 1)]
        psms = np.array(self.psms, dtype=dtype)
        key = (self.key(psm) for psm in psms)
        f11 = aux.filter(psms, key=key, is_decoy=self.is_decoy, fdr=0.5)
        self.assertEqual(f11.shape[0], 26)
        key = (self.key(psm) for psm in psms)
        with aux.filter(psms, key=key, is_decoy=self.is_decoy, fdr=0.5, full_output=False) as f:
            f11 = list(f)
        self.assertEqual(len(f11), 26)

    def test_filter_array_iter_key_str_is_decoy(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('is decoy', np.bool)]
        psms = np.array([(s, l, self.is_decoy((s, l))) for s, l in self.psms], dtype=dtype)
        key = iter([self.key(psm) for psm in psms])
        f11 = aux.filter(psms, key=key, is_decoy='is decoy', fdr=0.5)
        self.assertEqual(f11.shape[0], 26)
        key = iter(self.key(psm) for psm in psms)
        with aux.filter(psms, key=key, is_decoy='is decoy', fdr=0.5, full_output=False) as f:
            f11 = list(f)
        self.assertEqual(len(f11), 26)

    def test_filter_array_arr_is_decoy(self):
        dtype = [('score', np.int8), ('label', np.str_, 1)]
        psms = np.array(self.psms, dtype=dtype)
        is_decoy = np.array([self.is_decoy(psm) for psm in self.psms])
        self._run_check(psms, is_decoy=is_decoy)

    def test_filter_dataframe(self):
        dtype = [('score', np.int8), ('label', np.str_, 1)]
        psms = pd.DataFrame(np.array(self.psms, dtype=dtype))
        self._run_check(psms)

    def test_filter_dataframe_str_key(self):
        dtype = [('score', np.int8), ('label', np.str_, 1)]
        psms = pd.DataFrame(np.array(self.psms, dtype=dtype))
        self._run_check(psms, key='score')

    def test_filter_dataframe_str_is_decoy(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('is decoy', np.bool)]
        psms = np.array([(s, l, self.is_decoy((s, l))) for s, l in self.psms], dtype=dtype)
        psms = pd.DataFrame(psms)
        self._run_check(psms, is_decoy='is decoy')

    def test_filter_dataframe_str_key_str_is_decoy(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('is decoy', np.bool)]
        psms = np.array([(s, l, self.is_decoy((s, l))) for s, l in self.psms], dtype=dtype)
        psms = pd.DataFrame(psms)
        self._run_check(psms, key='score', is_decoy='is decoy')

    def test_filter_dataframe_arr_key_str_is_decoy(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('is decoy', np.bool)]
        psms = np.array([(s, l, self.is_decoy((s, l))) for s, l in self.psms], dtype=dtype)
        key = psms['score']
        psms = pd.DataFrame(psms)
        self._run_check(psms, key=key, is_decoy='is decoy')

    def test_filter_dataframe_arr_key(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('is decoy', np.bool)]
        psms = np.array([(s, l, self.is_decoy((s, l))) for s, l in self.psms], dtype=dtype)
        key = psms['score']
        psms = pd.DataFrame(psms)
        self._run_check(psms, key=key)

    def test_filter_dataframe_list_key_list_is_decoy(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('is decoy', np.bool)]
        psms = np.array([(s, l, self.is_decoy((s, l))) for s, l in self.psms], dtype=dtype)
        key = list(psms['score'])
        is_decoy = list(psms['is decoy'])
        psms = pd.DataFrame(psms)
        self._run_check(psms, key=key, is_decoy=is_decoy)

    def test_filter_two_lists(self):
        i = np.random.randint(1, len(self.psms)-1)
        psms1 = self.psms[:i]
        psms2 = self.psms[i:]
        self._run_check(psms1, psms2)

    def test_filter_two_arrays(self):
        i = np.random.randint(1, len(self.psms)-1)
        dtype = [('score', np.int8), ('label', np.str_, 1)]
        psms1 = np.array(self.psms[:i], dtype=dtype)
        psms2 = np.array(self.psms[i:], dtype=dtype)
        self._run_check(psms1, psms2)

    def test_filter_two_dataframes(self):
        i = np.random.randint(1, len(self.psms)-1)
        dtype = [('score', np.int8), ('label', np.str_, 1)]
        psms1 = pd.DataFrame(np.array(self.psms[:i], dtype=dtype))
        psms2 = pd.DataFrame(np.array(self.psms[i:], dtype=dtype))
        self._run_check(psms1, psms2)

    def test_filter_two_iters(self):
        i = np.random.randint(1, len(self.psms)-1)
        psms1 = iter(self.psms[:i])
        psms2 = iter(self.psms[i:])
        f11 = aux.filter(psms1, psms2, key=self.key, is_decoy=self.is_decoy, fdr=0.5)
        self.assertEqual(f11.shape[0], 26)
        psms1 = iter(self.psms[:i])
        psms2 = iter(self.psms[i:])
        with aux.filter(psms1, psms2, key=self.key, is_decoy=self.is_decoy, fdr=0.5, full_output=False) as f:
            f11 = list(f)
        self.assertEqual(len(f11), 26)

    def test_filter_iter(self):
        psms = iter(self.psms)
        f11 = aux.filter(psms, key=self.key, is_decoy=self.is_decoy, fdr=0.5)
        self.assertEqual(f11.shape[0], 26)
        psms = iter(self.psms)
        with aux.filter(psms, key=self.key, is_decoy=self.is_decoy, fdr=0.5, full_output=False) as f:
            f11 = list(f)
        self.assertEqual(len(f11), 26)

    def test_filter_two_iters_iter_key_iter_is_decoy(self):
        i = np.random.randint(1, len(self.psms)-1)
        psms1 = iter(self.psms[:i])
        psms2 = iter(self.psms[i:])
        key = iter(self.key(p) for p in self.psms)
        is_decoy = iter(self.is_decoy(p) for p in self.psms)
        f11 = aux.filter(psms1, psms2, key=key, is_decoy=is_decoy, fdr=0.5)
        self.assertEqual(f11.shape[0], 26)
        psms1 = iter(self.psms[:i])
        psms2 = iter(self.psms[i:])
        key = iter(self.key(p) for p in self.psms)
        is_decoy = iter(self.is_decoy(p) for p in self.psms)
        with aux.filter(psms1, psms2, key=key, is_decoy=is_decoy, fdr=0.5, full_output=False) as f:
            f11 = list(f)
        self.assertEqual(len(f11), 26)

    def test_filter_two_arrays_str_key(self):
        i = np.random.randint(1, len(self.psms)-1)
        dtype = [('score', np.int8), ('label', np.str_, 1)]
        psms1 = np.array(self.psms[:i], dtype=dtype)
        psms2 = np.array(self.psms[i:], dtype=dtype)
        self._run_check(psms1, psms2, key='score')

    def test_filter_two_dataframes_str_key(self):
        i = np.random.randint(1, len(self.psms)-1)
        dtype = [('score', np.int8), ('label', np.str_, 1)]
        psms1 = pd.DataFrame(np.array(self.psms[:i], dtype=dtype))
        psms2 = pd.DataFrame(np.array(self.psms[i:], dtype=dtype))
        self._run_check(psms1, psms2, key='score')

    def test_filter_two_arrays_str_key_arr_is_decoy(self):
        i = np.random.randint(1, len(self.psms)-1)
        dtype = [('score', np.int8), ('label', np.str_, 1)]
        psms1 = np.array(self.psms[:i], dtype=dtype)
        psms2 = np.array(self.psms[i:], dtype=dtype)
        is_decoy = np.array([self.is_decoy(p) for p in self.psms])
        self._run_check(psms1, psms2, key='score', is_decoy=is_decoy)

    def test_filter_two_dataframes_str_key_str_is_decoy(self):
        i = np.random.randint(1, len(self.psms)-1)
        dtype = [('score', np.int8), ('label', np.str_, 1), ('is decoy', np.bool)]
        psms = np.array([(s, l, self.is_decoy((s, l))) for s, l in self.psms], dtype=dtype)
        psms1 = pd.DataFrame(psms[:i])
        psms2 = pd.DataFrame(psms[i:])
        self._run_check(psms1, psms2, key='score', is_decoy='is decoy')

    def test_filter_two_dataframes_str_key_arr_is_decoy(self):
        i = np.random.randint(1, len(self.psms)-1)
        dtype = [('score', np.int8), ('label', np.str_, 1), ('is decoy', np.bool)]
        psms = np.array([(s, l, self.is_decoy((s, l))) for s, l in self.psms], dtype=dtype)
        is_decoy = psms['is decoy']
        psms1 = pd.DataFrame(psms[:i])
        psms2 = pd.DataFrame(psms[i:])
        self._run_check(psms1, psms2, key='score', is_decoy=is_decoy)

    def test_filter_two_dataframes_str_key_iter_is_decoy(self):
        i = np.random.randint(1, len(self.psms)-1)
        dtype = [('score', np.int8), ('label', np.str_, 1), ('is decoy', np.bool)]
        psms = np.array([(s, l, self.is_decoy((s, l))) for s, l in self.psms], dtype=dtype)
        is_decoy = iter(psms['is decoy'])
        psms1 = pd.DataFrame(psms[:i])
        psms2 = pd.DataFrame(psms[i:])
        f11 = aux.filter(psms1, psms2, key='score', is_decoy=is_decoy, fdr=0.5)
        self.assertEqual(f11.shape[0], 26)
        is_decoy = iter(psms['is decoy'])
        with aux.filter(psms1, psms2, key='score', is_decoy=is_decoy, fdr=0.5, full_output=False) as f:
            f11 = list(f)
        self.assertEqual(len(f11), 26)


if __name__ == '__main__':
    unittest.main()