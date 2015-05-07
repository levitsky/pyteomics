import unittest
import string
from itertools import count
import operator as op
import numpy as np
from pyteomics import auxiliary as aux

class QvalueTest(unittest.TestCase):
    def setUp(self):
        self.psms = zip(count(), string.ascii_uppercase + string.ascii_lowercase)
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
        q = aux.qvalues(self.psms, key=self.key, is_decoy=self.is_decoy)
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
        # print(q.dtype)
        self.assertTrue(q['psm'].dtype == dtype)

if __name__ == '__main__':
    unittest.main()