import unittest
import string
from itertools import count
import operator as op
import numpy as np
import pandas as pd
import tempfile

from os import path
import pyteomics
pyteomics.__path__ = [path.abspath(path.join(path.dirname(__file__), path.pardir, 'pyteomics'))]
from pyteomics import auxiliary as aux
from pyteomics import tandem
from pyteomics import version

psms = list(zip(count(), string.ascii_uppercase + string.ascii_lowercase,
            np.arange(0.01, 0.062, 0.001)))


class QvalueTest(unittest.TestCase):
    key = staticmethod(op.itemgetter(0))
    is_decoy = staticmethod(lambda x: x[1].islower())
    pep = staticmethod(op.itemgetter(2))

    def setUp(self):
        np.random.shuffle(psms)
        self.psms = iter(psms)

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
        spsms = sorted(self.psms, key=self.key)
        self.assertTrue(np.allclose([self.is_decoy(x) for x in spsms],
            q['is decoy']))
        self.assertTrue(np.allclose([self.key(x) for x in spsms],
            q['score']))
        self.setUp()

    def _run_check_pep(self, q):
        self.assertTrue(np.allclose(q['q'], np.arange(0.01, 0.036, 0.0005)))
        self.setUp()

    def test_qvalues(self):
        q = aux.qvalues(self.psms, key=self.key, is_decoy=self.is_decoy, remove_decoy=True)
        self.assertTrue(np.allclose(q['q'], 0))
        self.assertTrue(np.allclose(q['is decoy'], 0))
        self.assertTrue(np.allclose(q['score'], np.arange(26)))

    def test_qvalues_pep(self):
        q = aux.qvalues(self.psms, pep=self.pep)
        self._run_check_pep(q)
        q = aux.qvalues(self.psms, pep=self.pep, key=self.key)
        self._run_check_pep(q)

    def test_qvalues_with_decoy(self):
        q = aux.qvalues(self.psms, key=self.key, is_decoy=self.is_decoy, remove_decoy=False)
        self._run_check(q, 2)
        q = aux.qvalues(self.psms, key=self.key, is_decoy=self.is_decoy, remove_decoy=False, formula=1)
        self._run_check(q, 1)

    def test_qvalues_full_output(self):
        q = aux.qvalues(self.psms, key=self.key, is_decoy=self.is_decoy, remove_decoy=False, full_output=True)
        self._run_check(q, 2)

    def test_qvalues_pep_full_output(self):
        q = aux.qvalues(self.psms, pep=self.pep, full_output=True)
        self._run_check_pep(q)
        q = aux.qvalues(self.psms, key=self.key, pep=self.pep, full_output=True)
        self._run_check_pep(q)

    def test_qvalues_from_numpy(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64)]
        psms = np.array(list(self.psms), dtype=dtype)
        q = aux.qvalues(psms, key=self.key, is_decoy=self.is_decoy, remove_decoy=False, formula=1)
        self._run_check(q, 1)
        q = aux.qvalues(psms, key=self.key, is_decoy=self.is_decoy, remove_decoy=False, formula=1,
            full_output=True)
        self._run_check(q, 1)
        self.assertTrue(q['psm'].dtype == dtype)

    def test_qvalues_pep_from_numpy(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64)]
        psms = np.array(list(self.psms), dtype=dtype)
        q = aux.qvalues(psms, pep=self.pep)
        self._run_check_pep(q)
        q = aux.qvalues(psms, key=self.key, pep=self.pep, full_output=True)
        self._run_check_pep(q)
        self.assertTrue(q['psm'].dtype == dtype)

    def test_qvalues_from_dataframe(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64)]
        psms = pd.DataFrame(np.array(list(self.psms), dtype=dtype))
        q = aux.qvalues(psms, key=self.key, is_decoy=self.is_decoy, remove_decoy=False, formula=1)
        self._run_check(q, 1)
        q = aux.qvalues(psms, key=self.key, is_decoy=self.is_decoy, remove_decoy=False, formula=1, full_output=True)
        self._run_check(q, 1)

    def test_qvalues_empty_dataframe(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64)]
        psms = pd.DataFrame(np.array([], dtype=dtype))
        q = aux.qvalues(psms, key=self.key, is_decoy=self.is_decoy, remove_decoy=False, formula=1)
        self.assertEqual(q.shape[0], 0)
        q = aux.qvalues(psms, key=self.key, is_decoy=self.is_decoy, remove_decoy=False, formula=1, full_output=True)
        self.assertEqual(q.shape[0], 0)

    def test_qvalues_pep_from_dataframe(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64)]
        psms = pd.DataFrame(np.array(list(self.psms), dtype=dtype))
        q = aux.qvalues(psms, pep=self.pep)
        self._run_check_pep(q)
        q = aux.qvalues(psms, pep=self.pep, full_output=True)
        self._run_check_pep(q)

    def test_qvalues_from_numpy_string_key(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64)]
        psms = np.array(list(self.psms), dtype=dtype)
        q = aux.qvalues(psms, key='score', is_decoy=self.is_decoy, remove_decoy=False, formula=1)
        self._run_check(q, 1)
        q = aux.qvalues(psms, key='score', is_decoy=self.is_decoy, remove_decoy=False, formula=1,
            full_output=True)
        self._run_check(q, 1)
        self.assertTrue(q['psm'].dtype == dtype)

    def test_qvalues_pep_from_numpy_string_pep(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64)]
        psms = np.array(list(self.psms), dtype=dtype)
        q = aux.qvalues(psms, pep='pep')
        self._run_check_pep(q)
        q = aux.qvalues(psms, key='score', pep='pep')
        self._run_check_pep(q)
        q = aux.qvalues(psms, key='score', pep='pep', full_output=True)
        self._run_check_pep(q)

    def test_qvalues_from_dataframe_string_key(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64)]
        psms = pd.DataFrame(np.array(list(self.psms), dtype=dtype))
        q = aux.qvalues(psms, key='score', is_decoy=self.is_decoy, remove_decoy=False, formula=1)
        self._run_check(q, 1)
        q = aux.qvalues(psms, key='score', is_decoy=self.is_decoy, remove_decoy=False, formula=1,
            full_output=True)
        self._run_check(q, 1)

    def test_qvalues_pep_from_dataframe_string_pep(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64)]
        psms = pd.DataFrame(np.array(list(self.psms), dtype=dtype))
        q = aux.qvalues(psms, key=self.key, pep='pep')
        self._run_check_pep(q)
        q = aux.qvalues(psms, pep='pep')
        self._run_check_pep(q)
        q = aux.qvalues(psms, key='score', pep='pep', full_output=True)
        self._run_check_pep(q)

    def test_qvalues_from_dataframe_string_key_and_is_decoy(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64)]
        psms = pd.DataFrame(np.array(list(self.psms), dtype=dtype))
        psms['is decoy'] = [self.is_decoy(row) for _, row in psms.iterrows()]
        q = aux.qvalues(psms, key='score', is_decoy='is decoy', remove_decoy=False, formula=1)
        self._run_check(q, 1)
        q = aux.qvalues(psms, key='score', is_decoy='is decoy', remove_decoy=False, formula=1,
            full_output=True)
        self._run_check(q, 1)

    def test_qvalues_pep_from_dataframe_string_key_and_pep(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64)]
        psms = pd.DataFrame(np.array(list(self.psms), dtype=dtype))
        q = aux.qvalues(psms, key='score', pep='pep')
        self._run_check_pep(q)
        q = aux.qvalues(psms, key='score', pep='pep', full_output=True)
        self._run_check_pep(q)

    def test_qvalues_pep_exceptions(self):
        self.assertRaises(aux.PyteomicsError, aux.qvalues,
            self.psms, pep='pep', is_decoy=self.is_decoy)
        self.assertRaises(aux.PyteomicsError, aux.qvalues,
            self.psms, pep='pep', remove_decoy=False)
        self.assertRaises(aux.PyteomicsError, aux.qvalues,
            self.psms, pep='pep', correction=0)

    def test_qvalues_from_tandem(self):
        psms = tandem.TandemXML('test.t.xml')
        q0 = aux.qvalues(psms, key=op.itemgetter('expect'), is_decoy=tandem.is_decoy)
        with tandem.TandemXML('test.t.xml') as psms:
            q1 = aux.qvalues(psms, key=op.itemgetter('expect'), is_decoy=tandem.is_decoy)
        self.assertTrue(np.allclose(q0['q'], q1['q']))


class FilterTest(unittest.TestCase):
    key = staticmethod(op.itemgetter(0))
    is_decoy = staticmethod(lambda x: x[1].islower())
    pep = staticmethod(op.itemgetter(2))

    def setUp(self):
        self.psms = psms
        np.random.shuffle(self.psms)

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

    def _run_check_pep(self, *args, **kwargs):
        key = kwargs.pop('key', self.key)
        f11 = aux.filter(*args, key=key, fdr=0.02, **kwargs)
        f12 = aux.filter(*args, fdr=0.02, **kwargs)

        self.assertEqual(f11.shape[0], 21)
        self.assertEqual(f12.shape[0], 21)

        with aux.filter(*args, key=key, fdr=0.02, full_output=False, **kwargs) as f:
            f11 = list(f)
        with aux.filter(*args, fdr=0.02, full_output=False, **kwargs) as f:
            f12 = list(f)

        self.assertEqual(len(f11), 21)
        self.assertEqual(len(f12), 21)

    def test_filter(self):
        self._run_check(self.psms)

    def test_filter_pep(self):
        self._run_check_pep(self.psms, pep=self.pep)

    def test_filter_chain(self):
        f = aux.filter.chain(self.psms, self.psms, key=self.key, is_decoy=self.is_decoy, fdr=0.5)
        self.assertEqual(f.shape[0], 52)

    def test_filter_chain_pep(self):
        f = aux.filter.chain(self.psms, self.psms, pep=self.pep, fdr=0.02)
        self.assertEqual(f.shape[0], 42)

    def test_filter_chain_with(self):
        with aux.filter.chain(self.psms, self.psms, key=self.key, is_decoy=self.is_decoy,
            fdr=0.5, full_output=False) as f:
            f1 = list(f)
        self.assertEqual(len(f1), 52)

    def test_filter_pep_chain_with(self):
        with aux.filter.chain(self.psms, self.psms, pep=self.pep,
            fdr=0.02, full_output=False) as f:
            f1 = list(f)
        self.assertEqual(len(f1), 42)

    def test_filter_chain_arr_str_key(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64)]
        psms = np.array(self.psms, dtype=dtype)
        f11 = aux.filter.chain(psms, psms, key='score', is_decoy=self.is_decoy, fdr=0.5)
        self.assertEqual(f11.shape[0], 52)

    def test_filter_pep_chain_arr_str_key(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64)]
        psms = np.array(self.psms, dtype=dtype)
        f = aux.filter.chain(psms, psms, key='score', pep=self.pep, fdr=0.02)
        self.assertEqual(f.shape[0], 42)

    def test_filter_chain_arr_str_key_with(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64)]
        psms = np.array(self.psms, dtype=dtype)
        with aux.filter.chain(psms, psms, key='score', is_decoy=self.is_decoy, fdr=0.5, full_output=False) as f:
            f1 = list(f)
        self.assertEqual(len(f1), 52)

    def test_filter_pep_chain_arr_str_key_with(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64)]
        psms = np.array(self.psms, dtype=dtype)
        with aux.filter.chain(psms, psms, key='score', pep=self.pep, fdr=0.02, full_output=False) as f:
            f1 = list(f)
        self.assertEqual(len(f1), 42)

    def test_filter_chain_from_iterable(self):
        f11 = aux.filter.chain.from_iterable([self.psms, self.psms], key=self.key, is_decoy=self.is_decoy, fdr=0.5)
        self.assertEqual(f11.shape[0], 52)

    def test_filter_pep_chain_from_iterable(self):
        f = aux.filter.chain.from_iterable([self.psms, self.psms], pep=self.pep, fdr=0.02)
        self.assertEqual(f.shape[0], 42)

    def test_filter_chain_from_iterable_with(self):
        with aux.filter.chain.from_iterable([self.psms, self.psms],
            key=self.key, is_decoy=self.is_decoy, fdr=0.5, full_output=False) as f:
            f11 = list(f)
        self.assertEqual(len(f11), 52)

    def test_filter_pep_chain_from_iterable_with(self):
        with aux.filter.chain.from_iterable([self.psms, self.psms],
            key=self.key, pep=self.pep, fdr=0.02, full_output=False) as f:
            f1 = list(f)
        self.assertEqual(len(f1), 42)

    def test_filter_array(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64)]
        psms = np.array(self.psms, dtype=dtype)
        self._run_check(psms)

    def test_filter_pep_array(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64)]
        psms = np.array(self.psms, dtype=dtype)
        self._run_check_pep(psms, pep=self.pep)

    def test_filter_array_str_is_decoy(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64), ('is decoy', np.bool_)]
        psms = np.array([(s, l, p, self.is_decoy((s, l, p))) for s, l, p in self.psms], dtype=dtype)
        self._run_check(psms, is_decoy='is decoy')

    def test_filter_pep_array_str_pep(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64), ('is decoy', np.bool_)]
        psms = np.array([(s, l, p, self.is_decoy((s, l, p))) for s, l, p in self.psms], dtype=dtype)
        self._run_check_pep(psms, pep='pep')

    def test_filter_array_str_is_decoy_str_key(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64), ('is decoy', np.bool_)]
        psms = np.array([(s, l, p, self.is_decoy((s, l, p))) for s, l, p in self.psms], dtype=dtype)
        self._run_check(psms, is_decoy='is decoy', key='score')

    def test_filter_pep_array_str_pep_str_key(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64), ('is decoy', np.bool_)]
        psms = np.array([(s, l, p, self.is_decoy((s, l, p))) for s, l, p in self.psms], dtype=dtype)
        self._run_check_pep(psms, pep='pep', key='score')

    def test_filter_array_list_key(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64)]
        psms = np.array(self.psms, dtype=dtype)
        key = [self.key(psm) for psm in psms]
        self._run_check(psms, key=key)

    def test_filter_pep_array_list_key(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64)]
        psms = np.array(self.psms, dtype=dtype)
        key = [self.key(psm) for psm in psms]
        self._run_check_pep(psms, key=key, pep=self.pep)

    def test_filter_pep_array_list_pep_key(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64)]
        psms = np.array(self.psms, dtype=dtype)
        pep = [self.pep(psm) for psm in psms]
        self._run_check_pep(psms, pep=pep)

    def test_filter_array_gen_key(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64)]
        psms = np.array(self.psms, dtype=dtype)
        key = (self.key(psm) for psm in psms)
        f = aux.filter(psms, key=key, is_decoy=self.is_decoy, fdr=0.5)
        self.assertEqual(f.shape[0], 26)
        key = (self.key(psm) for psm in psms)
        with aux.filter(psms, key=key, is_decoy=self.is_decoy, fdr=0.5, full_output=False) as f:
            f1 = list(f)
        self.assertEqual(len(f1), 26)

    def test_filter_pep_array_gen_key(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64)]
        psms = np.array(self.psms, dtype=dtype)
        key = (self.key(psm) for psm in psms)
        f = aux.filter(psms, key=key, pep=self.pep, fdr=0.02)
        self.assertEqual(f.shape[0], 21)
        key = (self.key(psm) for psm in psms)
        with aux.filter(psms, key=key, pep=self.pep, fdr=0.02, full_output=False) as f:
            f11 = list(f)
        self.assertEqual(len(f11), 21)

    def test_filter_array_iter_key_str_is_decoy(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64), ('is decoy', np.bool_)]
        psms = np.array([(s, l, p, self.is_decoy((s, l, p))) for s, l, p in self.psms], dtype=dtype)
        key = iter([self.key(psm) for psm in psms])
        f11 = aux.filter(psms, key=key, is_decoy='is decoy', fdr=0.5)
        self.assertEqual(f11.shape[0], 26)
        key = iter(self.key(psm) for psm in psms)
        with aux.filter(psms, key=key, is_decoy='is decoy', fdr=0.5, full_output=False) as f:
            f11 = list(f)
        self.assertEqual(len(f11), 26)

    def test_filter_pep_array_iter_key_str_is_decoy(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64), ('is decoy', np.bool_)]
        psms = np.array([(s, l, p, self.is_decoy((s, l, p))) for s, l, p in self.psms], dtype=dtype)
        key = iter([self.key(psm) for psm in psms])
        f = aux.filter(psms, key=key, pep='pep', fdr=0.02)
        self.assertEqual(f.shape[0], 21)
        key = iter(self.key(psm) for psm in psms)
        with aux.filter(psms, key=key, pep='pep', fdr=0.02, full_output=False) as f:
            f1 = list(f)
        self.assertEqual(len(f1), 21)

    def test_filter_array_arr_is_decoy(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64)]
        psms = np.array(self.psms, dtype=dtype)
        is_decoy = np.array([self.is_decoy(psm) for psm in self.psms])
        self._run_check(psms, is_decoy=is_decoy)

    def test_filter_pep_array_arr_is_decoy(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64)]
        psms = np.array(self.psms, dtype=dtype)
        pep = np.array([self.pep(psm) for psm in self.psms])
        self._run_check_pep(psms, pep=pep)

    def test_filter_dataframe(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64)]
        psms = pd.DataFrame(np.array(self.psms, dtype=dtype))
        self._run_check(psms)

    def test_filter_empty_dataframe(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64)]
        psms = pd.DataFrame(np.array([], dtype=dtype))
        f = aux.filter(psms, key=self.key, is_decoy=self.is_decoy, remove_decoy=False, formula=1, fdr=0.1)
        self.assertEqual(f.shape[0], 0)
        f = aux.qvalues(psms, key=self.key, is_decoy=self.is_decoy, remove_decoy=False, formula=1, full_output=True, fdr=0.1)
        self.assertEqual(f.shape[0], 0)

    def test_filter_pep_dataframe(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64)]
        psms = pd.DataFrame(np.array(self.psms, dtype=dtype))
        self._run_check_pep(psms, pep=self.pep)

    def test_filter_dataframe_str_key(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64)]
        psms = pd.DataFrame(np.array(self.psms, dtype=dtype))
        self._run_check(psms, key='score')

    def test_filter_pep_dataframe_str_key(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64)]
        psms = pd.DataFrame(np.array(self.psms, dtype=dtype))
        self._run_check_pep(psms, key='score', pep=self.pep)

    def test_filter_dataframe_str_is_decoy(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64), ('is decoy', np.bool_)]
        psms = np.array([(s, l, p, self.is_decoy((s, l, p))) for s, l, p in self.psms], dtype=dtype)
        psms = pd.DataFrame(psms)
        self._run_check(psms, is_decoy='is decoy')

    def test_filter_pep_dataframe_str_pep(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64), ('is decoy', np.bool_)]
        psms = np.array([(s, l, p, self.is_decoy((s, l, p))) for s, l, p in self.psms], dtype=dtype)
        psms = pd.DataFrame(psms)
        self._run_check(psms, pep='pep', key=self.key)

    def test_filter_dataframe_str_key_str_is_decoy(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64), ('is decoy', np.bool_)]
        psms = np.array([(s, l, p, self.is_decoy((s, l, p))) for s, l, p in self.psms], dtype=dtype)
        psms = pd.DataFrame(psms)
        self._run_check(psms, key='score', is_decoy='is decoy')

    def test_filter_empty_dataframe_str_key_str_is_decoy(self):
        # dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64), ('is decoy', np.bool_)]
        psms = pd.DataFrame({'score': [], 'is decoy': []})
        f = aux.filter(psms, key='score', is_decoy='is decoy', fdr=0.1)
        self.assertEqual(f.shape[0], 0)
        f = aux.qvalues(psms, key='score', is_decoy='is decoy', remove_decoy=False, formula=1, full_output=True, fdr=0.01)
        self.assertEqual(f.shape[0], 0)

    def test_filter_pep_dataframe_str_key_str_pep(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64), ('is decoy', np.bool_)]
        psms = np.array([(s, l, p, self.is_decoy((s, l, p))) for s, l, p in self.psms], dtype=dtype)
        psms = pd.DataFrame(psms)
        self._run_check_pep(psms, key='score', pep='pep')

    def test_filter_dataframe_arr_key_str_is_decoy(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64), ('is decoy', np.bool_)]
        psms = np.array([(s, l, p, self.is_decoy((s, l, p))) for s, l, p in self.psms], dtype=dtype)
        key = psms['score']
        psms = pd.DataFrame(psms)
        self._run_check(psms, key=key, is_decoy='is decoy')

    def test_filter_pep_dataframe_arr_key_str_is_decoy(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64), ('is decoy', np.bool_)]
        psms = np.array([(s, l, p, self.is_decoy((s, l, p))) for s, l, p in self.psms], dtype=dtype)
        key = psms['score']
        psms = pd.DataFrame(psms)
        self._run_check(psms, key=key, pep='pep')

    def test_filter_dataframe_arr_key(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64), ('is decoy', np.bool_)]
        psms = np.array([(s, l, p, self.is_decoy((s, l, p))) for s, l, p in self.psms], dtype=dtype)
        key = psms['score']
        psms = pd.DataFrame(psms)
        self._run_check(psms, key=key)

    def test_filter_pep_dataframe_arr_key(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64), ('is decoy', np.bool_)]
        psms = np.array([(s, l, p, self.is_decoy((s, l, p))) for s, l, p in self.psms], dtype=dtype)
        key = psms['score']
        psms = pd.DataFrame(psms)
        self._run_check_pep(psms, key=key, pep=self.pep)

    def test_filter_dataframe_list_key_list_is_decoy(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64), ('is decoy', np.bool_)]
        psms = np.array([(s, l, p, self.is_decoy((s, l, p))) for s, l, p in self.psms], dtype=dtype)
        key = list(psms['score'])
        is_decoy = list(psms['is decoy'])
        psms = pd.DataFrame(psms)
        self._run_check(psms, key=key, is_decoy=is_decoy)

    def test_filter_pep_dataframe_list_key_list_pep(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64), ('is decoy', np.bool_)]
        psms = np.array([(s, l, p, self.is_decoy((s, l, p))) for s, l, p in self.psms], dtype=dtype)
        key = list(psms['score'])
        pep = list(psms['pep'])
        psms = pd.DataFrame(psms)
        self._run_check(psms, key=key, pep=pep)

    def test_filter_two_lists(self):
        i = np.random.randint(1, len(self.psms)-1)
        psms1 = self.psms[:i]
        psms2 = self.psms[i:]
        self._run_check(psms1, psms2)

    def test_filter_pep_two_lists(self):
        i = np.random.randint(1, len(self.psms)-1)
        psms1 = self.psms[:i]
        psms2 = self.psms[i:]
        self._run_check_pep(psms1, psms2, pep=self.pep)

    def test_filter_two_arrays(self):
        i = np.random.randint(1, len(self.psms)-1)
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64)]
        psms1 = np.array(self.psms[:i], dtype=dtype)
        psms2 = np.array(self.psms[i:], dtype=dtype)
        self._run_check(psms1, psms2)

    def test_filter_pep_two_arrays(self):
        i = np.random.randint(1, len(self.psms)-1)
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64)]
        psms1 = np.array(self.psms[:i], dtype=dtype)
        psms2 = np.array(self.psms[i:], dtype=dtype)
        self._run_check_pep(psms1, psms2, pep=self.pep)

    def test_filter_two_dataframes(self):
        i = np.random.randint(1, len(self.psms)-1)
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64)]
        psms1 = pd.DataFrame(np.array(self.psms[:i], dtype=dtype))
        psms2 = pd.DataFrame(np.array(self.psms[i:], dtype=dtype))
        self._run_check(psms1, psms2)

    def test_filter_pep_two_dataframes(self):
        i = np.random.randint(1, len(self.psms)-1)
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64)]
        psms1 = pd.DataFrame(np.array(self.psms[:i], dtype=dtype))
        psms2 = pd.DataFrame(np.array(self.psms[i:], dtype=dtype))
        self._run_check_pep(psms1, psms2, pep=self.pep)

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

    def test_filter_pep_two_iters(self):
        i = np.random.randint(1, len(self.psms)-1)
        psms1 = iter(self.psms[:i])
        psms2 = iter(self.psms[i:])
        f = aux.filter(psms1, psms2, key=self.key, pep=self.pep, fdr=0.02)
        self.assertEqual(f.shape[0], 21)
        psms1 = iter(self.psms[:i])
        psms2 = iter(self.psms[i:])
        with aux.filter(psms1, psms2, key=self.key, pep=self.pep, fdr=0.02, full_output=False) as f:
            f1 = list(f)
        self.assertEqual(len(f1), 21)

    def test_filter_iter(self):
        psms = iter(self.psms)
        f = aux.filter(psms, key=self.key, is_decoy=self.is_decoy, fdr=0.5)
        self.assertEqual(f.shape[0], 26)
        psms = iter(self.psms)
        with aux.filter(psms, key=self.key, is_decoy=self.is_decoy, fdr=0.5, full_output=False) as f:
            f1 = list(f)
        self.assertEqual(len(f1), 26)

    def test_filter_pep_iter(self):
        psms = iter(self.psms)
        f = aux.filter(psms, key=self.key, pep=self.pep, fdr=0.02)
        self.assertEqual(f.shape[0], 21)
        psms = iter(self.psms)
        with aux.filter(psms, key=self.key, pep=self.pep, fdr=0.02, full_output=False) as f:
            f1 = list(f)
        self.assertEqual(len(f1), 21)

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

    def test_filter_pep_two_iters_iter_key_iter_is_decoy(self):
        i = np.random.randint(1, len(self.psms)-1)
        psms1 = iter(self.psms[:i])
        psms2 = iter(self.psms[i:])
        key = iter(self.key(p) for p in self.psms)
        pep = iter(self.pep(p) for p in self.psms)
        f = aux.filter(psms1, psms2, key=key, pep=pep, fdr=0.02)
        self.assertEqual(f.shape[0], 21)
        psms1 = iter(self.psms[:i])
        psms2 = iter(self.psms[i:])
        key = iter(self.key(p) for p in self.psms)
        pep = iter(self.pep(p) for p in self.psms)
        with aux.filter(psms1, psms2, key=key, pep=pep, fdr=0.02, full_output=False) as f:
            f1 = list(f)
        self.assertEqual(len(f1), 21)

    def test_filter_two_arrays_str_key(self):
        i = np.random.randint(1, len(self.psms)-1)
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64)]
        psms1 = np.array(self.psms[:i], dtype=dtype)
        psms2 = np.array(self.psms[i:], dtype=dtype)
        self._run_check(psms1, psms2, key='score')

    def test_filter_pep_two_arrays_str_key_str_pep(self):
        i = np.random.randint(1, len(self.psms)-1)
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64)]
        psms1 = np.array(self.psms[:i], dtype=dtype)
        psms2 = np.array(self.psms[i:], dtype=dtype)
        self._run_check_pep(psms1, psms2, key='score', pep='pep')

    def test_filter_two_dataframes_str_key(self):
        i = np.random.randint(1, len(self.psms)-1)
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64)]
        psms1 = pd.DataFrame(np.array(self.psms[:i], dtype=dtype))
        psms2 = pd.DataFrame(np.array(self.psms[i:], dtype=dtype))
        self._run_check(psms1, psms2, key='score')

    def test_filter_pep_two_dataframes_str_key(self):
        i = np.random.randint(1, len(self.psms)-1)
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64)]
        psms1 = pd.DataFrame(np.array(self.psms[:i], dtype=dtype))
        psms2 = pd.DataFrame(np.array(self.psms[i:], dtype=dtype))
        self._run_check_pep(psms1, psms2, key='score', pep=self.pep)

    def test_filter_two_arrays_str_key_arr_is_decoy(self):
        i = np.random.randint(1, len(self.psms)-1)
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64)]
        psms1 = np.array(self.psms[:i], dtype=dtype)
        psms2 = np.array(self.psms[i:], dtype=dtype)
        is_decoy = np.array([self.is_decoy(p) for p in self.psms])
        self._run_check(psms1, psms2, key='score', is_decoy=is_decoy)

    def test_filter_pep_two_arrays_str_key_arr_pep(self):
        i = np.random.randint(1, len(self.psms)-1)
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64)]
        psms1 = np.array(self.psms[:i], dtype=dtype)
        psms2 = np.array(self.psms[i:], dtype=dtype)
        pep = np.array([self.pep(p) for p in self.psms])
        self._run_check_pep(psms1, psms2, key='score', pep=pep)

    def test_filter_two_dataframes_str_key_str_is_decoy(self):
        i = np.random.randint(1, len(self.psms)-1)
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64), ('is decoy', np.bool_)]
        psms = np.array([(s, l, p, self.is_decoy((s, l, p))) for s, l, p in self.psms], dtype=dtype)
        psms1 = pd.DataFrame(psms[:i])
        psms2 = pd.DataFrame(psms[i:])
        self._run_check(psms1, psms2, key='score', is_decoy='is decoy')

    def test_filter_pep_two_dataframes_str_key_str_pep(self):
        i = np.random.randint(1, len(self.psms)-1)
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64), ('is decoy', np.bool_)]
        psms = np.array([(s, l, p, self.is_decoy((s, l, p))) for s, l, p in self.psms], dtype=dtype)
        psms1 = pd.DataFrame(psms[:i])
        psms2 = pd.DataFrame(psms[i:])
        self._run_check_pep(psms1, psms2, key='score', pep='pep')

    def test_filter_two_dataframes_str_key_arr_is_decoy(self):
        i = np.random.randint(1, len(self.psms)-1)
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64), ('is decoy', np.bool_)]
        psms = np.array([(s, l, p, self.is_decoy((s, l, p))) for s, l, p in self.psms], dtype=dtype)
        is_decoy = psms['is decoy']
        psms1 = pd.DataFrame(psms[:i])
        psms2 = pd.DataFrame(psms[i:])
        self._run_check(psms1, psms2, key='score', is_decoy=is_decoy)

    def test_filter_pep_two_dataframes_str_key_arr_pep(self):
        i = np.random.randint(1, len(self.psms)-1)
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64), ('is decoy', np.bool_)]
        psms = np.array([(s, l, p, self.is_decoy((s, l, p))) for s, l, p in self.psms], dtype=dtype)
        pep = psms['pep']
        psms1 = pd.DataFrame(psms[:i])
        psms2 = pd.DataFrame(psms[i:])
        self._run_check_pep(psms1, psms2, key='score', pep=pep)

    def test_filter_two_dataframes_str_key_iter_is_decoy(self):
        i = np.random.randint(1, len(self.psms)-1)
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64), ('is decoy', np.bool_)]
        psms = np.array([(s, l, p, self.is_decoy((s, l, p))) for s, l, p in self.psms], dtype=dtype)
        is_decoy = iter(psms['is decoy'])
        psms1 = pd.DataFrame(psms[:i])
        psms2 = pd.DataFrame(psms[i:])
        f11 = aux.filter(psms1, psms2, key='score', is_decoy=is_decoy, fdr=0.5)
        self.assertEqual(f11.shape[0], 26)
        is_decoy = iter(psms['is decoy'])
        with aux.filter(psms1, psms2, key='score', is_decoy=is_decoy, fdr=0.5, full_output=False) as f:
            f11 = list(f)
        self.assertEqual(len(f11), 26)

    def test_filter_pep_two_dataframes_str_key_iter_pep(self):
        i = np.random.randint(1, len(self.psms)-1)
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64), ('is decoy', np.bool_)]
        psms = np.array([(s, l, p, self.is_decoy((s, l, p))) for s, l, p in self.psms], dtype=dtype)
        pep = iter(psms['pep'])
        psms1 = pd.DataFrame(psms[:i])
        psms2 = pd.DataFrame(psms[i:])
        f = aux.filter(psms1, psms2, key='score', pep=pep, fdr=0.02)
        self.assertEqual(f.shape[0], 21)
        pep = iter(psms['pep'])
        with aux.filter(psms1, psms2, key='score', pep=pep, fdr=0.02, full_output=False) as f:
            f1 = list(f)
        self.assertEqual(len(f1), 21)


class FDRTest(unittest.TestCase):

    is_decoy = staticmethod(lambda x: x[1].islower())
    pep = staticmethod(op.itemgetter(2))

    def _run_check(self, psms, **kwargs):
        is_decoy = kwargs.pop('is_decoy', self.is_decoy)
        pep = kwargs.pop('pep', self.pep)
        self.assertAlmostEqual(aux.fdr(psms, is_decoy=is_decoy, formula=1), 1.0)
        self.assertAlmostEqual(aux.fdr(psms, is_decoy=is_decoy, formula=2), 1.0)
        self.assertAlmostEqual(aux.fdr(psms, pep=pep), 0.0355)

    def test_fdr(self):
        self._run_check(psms)

    def test_fdr_iter(self):
        self.assertAlmostEqual(aux.fdr(iter(psms), is_decoy=self.is_decoy), 1.0)
        self.assertAlmostEqual(aux.fdr(iter(psms), pep=self.pep), 0.0355)
        isd = [self.is_decoy((s, l, p)) for s, l, p in psms]
        pep = [self.pep((s, l, p)) for s, l, p in psms]
        self.assertAlmostEqual(aux.fdr(iter(psms), is_decoy=iter(isd)), 1.0)
        self.assertAlmostEqual(aux.fdr(iter(psms), pep=iter(pep)), 0.0355)

    def test_fdr_array_str(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64), ('is decoy', np.bool_)]
        psms_ = np.array([(s, l, p, self.is_decoy((s, l, p))) for s, l, p in psms], dtype=dtype)
        self._run_check(psms_, is_decoy='is decoy', pep='pep')

    def test_fdr_df_str(self):
        dtype = [('score', np.int8), ('label', np.str_, 1), ('pep', np.float64), ('is decoy', np.bool_)]
        psms_ = np.array([(s, l, p, self.is_decoy((s, l, p))) for s, l, p in psms], dtype=dtype)
        psms1 = pd.DataFrame(psms_)
        self._run_check(psms1, is_decoy='is decoy', pep='pep')

    def test_fdr_no_psms(self):
        isd = [self.is_decoy((s, l, p)) for s, l, p in psms]
        pep = [self.pep((s, l, p)) for s, l, p in psms]
        self._run_check(None, is_decoy=isd, pep=pep)

    def test_fdr_list(self):
        isd = [self.is_decoy((s, l, p)) for s, l, p in psms]
        pep = [self.pep((s, l, p)) for s, l, p in psms]
        self._run_check(psms, is_decoy=isd, pep=pep)

    def test_sigma_T(self):
        self.assertAlmostEqual(aux.sigma_T(psms, is_decoy=self.is_decoy), 7.348469228)

    def test_sigma_fdr(self):
        self.assertAlmostEqual(aux.sigma_fdr(psms, is_decoy=self.is_decoy), 0.28263343)


class RegressionTest(unittest.TestCase):
    x = [1, 2, 3]
    y = [3, 5, 7]
    a = 2
    b = 1
    r = 1
    stderr = 0

    def _test_linreg(self, result):
        a, b, r, stderr = result
        self.assertAlmostEqual(a, self.a)
        self.assertAlmostEqual(b, self.b)
        self.assertAlmostEqual(r, self.r)
        self.assertAlmostEqual(stderr, self.stderr)

    def test_linear_regression_simple(self):
        result = aux.linear_regression(self.x, self.y)
        self._test_linreg(result)

    def test_linear_regression_simple_vertical(self):
        result = aux.linear_regression_vertical(self.x, self.y)
        self._test_linreg(result)

    def test_linear_regression_simple_perpendicular(self):
        result = aux.linear_regression_perpendicular(self.x, self.y)
        self._test_linreg(result)

    def test_linear_regression_no_y_list(self):
        x = list(zip(self.x, self.y))
        result = aux.linear_regression(x)
        self._test_linreg(result)

    def test_linear_regression_no_y_list_vertical(self):
        x = list(zip(self.x, self.y))
        result = aux.linear_regression_vertical(x)
        self._test_linreg(result)

    def test_linear_regression_no_y_list_perpendicular(self):
        x = list(zip(self.x, self.y))
        result = aux.linear_regression_perpendicular(x)
        self._test_linreg(result)

    def test_linear_regression_no_y_arr(self):
        x = np.array(list(zip(self.x, self.y)))
        result = aux.linear_regression(x)
        self._test_linreg(result)

    def test_linear_regression_no_y_arr_vertical(self):
        x = np.array(list(zip(self.x, self.y)))
        result = aux.linear_regression_vertical(x)
        self._test_linreg(result)

    def test_linear_regression_no_y_arr_perpendicular(self):
        x = np.array(list(zip(self.x, self.y)))
        result = aux.linear_regression_perpendicular(x)
        self._test_linreg(result)

    def test_linear_regression_shape_exception_vertical(self):
        with self.assertRaises(aux.PyteomicsError):
            aux.linear_regression_vertical(self.x)

    def test_linear_regression_shape_exception_perpendicular(self):
        with self.assertRaises(aux.PyteomicsError):
            aux.linear_regression_perpendicular(self.x)


class OffsetIndexTest(unittest.TestCase):
    def setUp(self):
        self.sequence = [(str(i), i) for i in range(10)]
        self.index = aux.OffsetIndex(self.sequence)

    def test_index_sequence(self):
        self.assertEqual(self.index.index_sequence, tuple(self.sequence))

    def test_find(self):
        self.assertEqual(self.index.find('3'), 3)

    def test_from_index(self):
        self.assertEqual(self.index.from_index(3), '3')
        self.assertEqual(self.index.from_index(4, True), ('4', 4))

    def test_from_slice(self):
        self.assertEqual(self.index.from_slice(slice(1, 3)), ['1', '2'])
        self.assertEqual(self.index.from_slice(slice(1, 3), True), self.sequence[1:3])

    def test_between(self):
        self.assertEqual(self.index.between('1', '3'), ['1', '2', '3'])
        self.assertEqual(self.index.between('1', '3', True), [('1', 1), ('2', 2), ('3', 3)])
        self.assertEqual(self.index.between('3', '1'), ['1', '2', '3'])
        self.assertEqual(self.index.between(None, '3'), ['0', '1', '2', '3'])
        self.assertEqual(self.index.between('8', None), ['8', '9'])


class UseIndexTest(unittest.TestCase):
    class MockFile:
        def __init__(self, seekable, mode):
            if mode is not None:
                self.mode = mode
            if seekable is not None:
                self._seekable = seekable

        def seekable(self):
            return getattr(self, '_seekable', None)

    def _check_file_object(self, fo, value):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            result = aux._check_use_index(fo, None, None)
            self.assertEqual(len(w), 0)
        self.assertEqual(result, value)

    def test_str_name(self):
        for ui in [False, True]:
            for default in [False, True]:
                self.assertEqual(aux._check_use_index('test.mgf', ui, default), ui)

    def test_textfile(self):
        with open('test.fasta') as f:
            self._check_file_object(f, False)

    def test_binfile(self):
        with open('test.mgf', 'rb') as f:
            self._check_file_object(f, True)

    def test_tmpfile_text(self):
        with tempfile.TemporaryFile(mode='r') as f:
            self._check_file_object(f, False)

    def test_tmpfile_bin(self):
        with tempfile.TemporaryFile(mode='wb') as f:
            self._check_file_object(f, True)

    def test_stringio(self):
        try:
            from io import StringIO
        except ImportError:
            from StringIO import StringIO
        with StringIO(u'test') as f:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter('always')
                aux._check_use_index(f, None, None)
                self.assertEqual(len(w), 1)
                self.assertIs(w[0].category, UserWarning)

    def test_error_not_seekable(self):
        source = UseIndexTest.MockFile(False, 'rb')
        self.assertRaises(aux.PyteomicsError, aux._check_use_index, source, None, None)

    def test_warning_not_seekable(self):
        source = UseIndexTest.MockFile(False, 'r')
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            aux._check_use_index(source, True, None)
            self.assertEqual(len(w), 1)
            self.assertIs(w[0].category, UserWarning)
            self.assertIn('is not seekable', str(w[0].message))

    def test_warning_wrong_mode(self):
        for m in ['rb', 'r']:
            source = UseIndexTest.MockFile(True, m)
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter('always')
                aux._check_use_index(source, 'b' not in m, None)
                self.assertEqual(len(w), 1)
                self.assertIs(w[0].category, UserWarning)
                self.assertIn('Setting `use_index` to {}'.format('b' in m), str(w[0].message))

    def test_warning_no_mode(self):
        source = UseIndexTest.MockFile(None, None)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            aux._check_use_index(source, None, None)
            self.assertEqual(len(w), 1)
            self.assertIs(w[0].category, UserWarning)
            self.assertIn('Could not check mode', str(w[0].message))


class VersionTest(unittest.TestCase):

    def test_short_version(self):
        self.assertEqual(version.VersionInfo('1.2'), ('1', '2', None, None, None))

    def test_longer_version(self):
        self.assertEqual(version.VersionInfo('1.2.3'), ('1', '2', '3', None, None))

    def test_short_dev_version(self):
        self.assertEqual(version.VersionInfo('1.2dev3'), ('1', '2', None, 'dev', '3'))

    def test_longer_dev_version(self):
        self.assertEqual(version.VersionInfo('1.2.3dev4'), ('1', '2', '3', 'dev', '4'))


import warnings
if __name__ == '__main__':
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        unittest.main()
