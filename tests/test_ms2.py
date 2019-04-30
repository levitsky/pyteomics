from os import path
import numpy as np
import pyteomics
pyteomics.__path__ = [path.abspath(path.join(path.dirname(__file__), path.pardir, 'pyteomics'))]
import unittest
from pyteomics.ms2 import read, read_header, MS2, IndexedMS2, chain
import data

class MS2Test(unittest.TestCase):
    maxDiff = None
    def setUp(self):
        self.path = 'test.ms2'
        self.header = read_header(self.path)
        self.spectra = list(read(self.path))
        self.ns = len(self.spectra)

    def test_read(self):
        # http://stackoverflow.com/q/14246983/1258041
        self.assertEqual(data.ms2_spectra, list(read(self.path)))
        for reader in [read, MS2, IndexedMS2, chain]:
            with reader(self.path) as reader:
                self.assertEqual(data.ms2_spectra, list(reader))

    def test_read_array_conversion(self):
        with read(self.path, convert_arrays=False) as reader:
            self.assertEqual(data.ms2_spectra_lists, list(reader))
        with read(self.path, convert_arrays=True) as reader:
            s = next(reader)
            self.assertTrue(isinstance(s['m/z array'], np.ndarray))

    def test_header(self):
        self.assertEqual(self.header, data.ms2_header)

    def test_read_dtype(self):
        dtypes = {'m/z array': np.float32, 'intensity array': np.int32}
        with read(self.path, dtype=dtypes) as f:
            for spec in f:
                for k, v in dtypes.items():
                    self.assertEqual(spec[k].dtype, v)

if __name__ == "__main__":
    unittest.main()