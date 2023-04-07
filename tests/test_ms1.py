import os
import numpy as np
import pyteomics
pyteomics.__path__ = [os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir, 'pyteomics'))]
import unittest
import pickle
from pyteomics.ms1 import read, read_header, MS1, IndexedMS1, chain
import data

class MS1Test(unittest.TestCase):
    maxDiff = None
    def setUp(self):
        self.path = 'test.ms1'
        self.header = read_header(self.path)
        self.spectra = list(read(self.path))
        self.ns = len(self.spectra)

    def test_read(self):
        # http://stackoverflow.com/q/14246983/1258041
        self.assertEqual(data.ms1_spectra, list(read(self.path)))
        for reader in [read, MS1, IndexedMS1, chain]:
            with reader(self.path) as reader:
                self.assertEqual(data.ms1_spectra, list(reader))

    def test_read_array_conversion(self):
        with read(self.path, convert_arrays=False) as reader:
            self.assertEqual(data.ms1_spectra_lists, list(reader))
        with read(self.path, convert_arrays=True) as reader:
            s = next(reader)
            self.assertTrue(isinstance(s['m/z array'], np.ndarray))

    def test_header(self):
        self.assertEqual(self.header, data.ms1_header)

    def test_read_dtype(self):
        dtypes = {'m/z array': np.float32, 'intensity array': np.int32}
        with read(self.path, dtype=dtypes) as f:
            for spec in f:
                for k, v in dtypes.items():
                    self.assertEqual(spec[k].dtype, v)

    def test_indexedms1_picklable(self):
        with IndexedMS1(self.path, block_size=12345, dtype=np.float32) as reader:
            spec = pickle.dumps(reader)
        with pickle.loads(spec) as reader:
            self.assertEqual(reader.block_size, 12345)
            self.assertEqual(reader._dtype_dict['m/z array'], np.float32)
            self.assertEqual(data.ms1_spectra, list(reader))

        with IndexedMS1(self.path, use_header=True) as reader:
            spec = pickle.dumps(reader)
        with pickle.loads(spec) as reader:
            self.assertEqual(data.ms1_header, reader.header)

    def test_ms1_picklable(self):
        with MS1(self.path, convert_arrays=0) as reader:
            spec = pickle.dumps(reader)
        with pickle.loads(spec) as reader:
            self.assertEqual(reader._convert_arrays, 0)
            self.assertEqual(data.ms1_spectra_lists, list(reader))


if __name__ == "__main__":
    unittest.main()
