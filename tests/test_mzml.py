from os import path
import pyteomics
pyteomics.__path__ = [path.abspath(path.join(path.dirname(__file__), path.pardir, 'pyteomics'))]
from itertools import product
import unittest
from pyteomics.mzml import *
from data import mzml_spectra
import numpy as np

class MzmlTest(unittest.TestCase):
    maxDiff = None
    path = 'test.mzML'

    def testReadSpectrum(self):
        for rs, it, ui in product([True, False], repeat=3):
            for func in [MzML, read, chain,
                    lambda x, **kw: chain.from_iterable([x], **kw), PreIndexedMzML]:
                with func(self.path, read_schema=rs, iterative=it, use_index=ui) as r:
                    # http://stackoverflow.com/q/14246983/1258041
                    self.assertEqual(mzml_spectra, list(r))

    def test_read_dtype(self):
        dtypes = {'m/z array': np.float32, 'intensity array': np.int32}
        with read(self.path, dtype=dtypes) as f:
            for spec in f:
                for k, v in dtypes.items():
                    self.assertEqual(spec[k].dtype, v)

    def test_has_built_index(self):
        with read(self.path, use_index=True) as f:
            self.assertGreater(len(f._offset_index), 0)
        with read(self.path, use_index=False) as f:
            self.assertEqual(len(f._offset_index), 0)

    def test_unit_extract(self):
        with MzML(self.path) as handle:
            for scan in handle:
                scan_inst = scan['scanList']['scan'][0]
                scan_time = scan_inst['scan start time']
                scan_window_lower_limit = scan_inst['scanWindowList']['scanWindow'][0]['scan window lower limit']
                self.assertEqual(scan_time.unit_info, 'minute')
                self.assertEqual(scan_window_lower_limit.unit_info, 'm/z')

if __name__ == '__main__':
    unittest.main()
