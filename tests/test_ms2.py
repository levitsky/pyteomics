import os
import numpy as np
import pyteomics
pyteomics.__path__ = [os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir, 'pyteomics'))]
import unittest
import copy
import pickle
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

    def test_read_no_charges(self):
        with read(self.path, convert_arrays=False, read_charges=False) as reader:
            lhs = copy.deepcopy(data.ms2_spectra_lists)
            for spec in lhs:
                del spec['charge array']
            self.assertEqual(lhs, list(reader))

        with read(self.path, convert_arrays=1, read_charges=False) as reader:
            lhs = copy.deepcopy(data.ms2_spectra)
            for spec in lhs:
                del spec['charge array']
            self.assertEqual(lhs, list(reader))

    def test_read_no_resolution(self):
        with read(self.path, convert_arrays=False, read_resolutions=False) as reader:
            lhs = copy.deepcopy(data.ms2_spectra_lists)
            for spec in lhs:
                del spec['resolution array']
            self.assertEqual(lhs, list(reader))

        with read(self.path, convert_arrays=1, read_resolutions=False) as reader:
            lhs = copy.deepcopy(data.ms2_spectra)
            for spec in lhs:
                del spec['resolution array']
            self.assertEqual(lhs, list(reader))

    def test_read_array_conversion(self):
        with read(self.path, convert_arrays=0) as reader:
            self.assertEqual(data.ms2_spectra_lists, list(reader))
        with read(self.path, convert_arrays=1) as reader:
            s = next(reader)
            self.assertTrue(isinstance(s['m/z array'], np.ndarray))
        with read(self.path, convert_arrays=2) as reader:
            s = next(reader)
            self.assertTrue(isinstance(s['m/z array'], np.ndarray))
            self.assertTrue(isinstance(s['charge array'], np.ma.core.MaskedArray))

    def test_header(self):
        self.assertEqual(self.header, data.ms2_header)

    def test_read_dtype(self):
        dtypes = {'m/z array': np.float32, 'intensity array': np.int32}
        with read(self.path, dtype=dtypes) as f:
            for spec in f:
                for k, v in dtypes.items():
                    self.assertEqual(spec[k].dtype, v)

    def test_indexedms2_picklable(self):
        with IndexedMS2(self.path, block_size=12345, convert_arrays=1, read_charges=False, read_resolutions=False) as reader:
            spec = pickle.dumps(reader)
        with pickle.loads(spec) as reader:
            self.assertEqual(reader.block_size, 12345)
            self.assertEqual(reader._read_charges, False)
            self.assertEqual(reader._read_resolutions, False)
            lhs = copy.deepcopy(data.ms2_spectra)
            for spec in lhs:
                del spec['resolution array']
                del spec['charge array']
            self.assertEqual(lhs, list(reader))

        with IndexedMS2(self.path, use_header=True) as reader:
            spec = pickle.dumps(reader)
        with pickle.loads(spec) as reader:
            self.assertEqual(data.ms2_header, reader.header)

    def test_ms2_picklable(self):
        with MS2(self.path, convert_arrays=1, read_charges=False, read_resolutions=False) as reader:
            spec = pickle.dumps(reader)
        with pickle.loads(spec) as reader:
            self.assertEqual(reader._read_charges, False)
            self.assertEqual(reader._read_resolutions, False)
            lhs = copy.deepcopy(data.ms2_spectra)
            for spec in lhs:
                del spec['resolution array']
                del spec['charge array']
            self.assertEqual(lhs, list(reader))


if __name__ == "__main__":
    unittest.main()
