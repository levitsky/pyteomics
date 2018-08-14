from os import path
import numpy as np
import pyteomics
pyteomics.__path__ = [path.abspath(path.join(path.dirname(__file__), path.pardir, 'pyteomics'))]
import tempfile
import unittest
import pickle
from pyteomics import mgf
import data

class MGFTest(unittest.TestCase):
    maxDiff = None
    _encoding = 'utf-8'
    def setUp(self):
        self.path = 'test.mgf'
        self.header = mgf.read_header(self.path)
        self.spectra = list(mgf.read(self.path))
        self.tmpfile = tempfile.TemporaryFile(mode='r+')
        mgf.write(header=self.header, spectra=self.spectra, output=self.tmpfile)
        self.tmpfile.seek(0)
        self.header2 = mgf.read_header(self.tmpfile)
        self.tmpfile.seek(0)
        tmpreader = mgf.read(self.tmpfile)
        self.spectra2 = list(tmpreader)
        self.ns = len(self.spectra)
        self.tmpfile.close()

    def test_read(self):
        for func in [mgf.read, mgf.MGF]:
            # http://stackoverflow.com/q/14246983/1258041
            self.assertEqual(data.mgf_spectra_long, list(func(self.path)))
            self.assertEqual(data.mgf_spectra_short, list(func(self.path, False)))
            with func(self.path) as reader:
                self.assertEqual(data.mgf_spectra_long, list(reader))
            with func(self.path, False) as reader:
                self.assertEqual(data.mgf_spectra_short, list(reader))

    def test_read_decoding(self):
        for func in [mgf.read, mgf.MGF, mgf.IndexedMGF]:
            self.assertEqual(data.mgf_spectra_long_decoded,
                list(func(self.path, encoding=self._encoding)))
            self.assertEqual(data.mgf_spectra_short_decoded,
                list(func(self.path, False, encoding=self._encoding)))
            with func(self.path, encoding=self._encoding) as reader:
                self.assertEqual(data.mgf_spectra_long_decoded, list(reader))
            with func(self.path, False, encoding=self._encoding) as reader:
                self.assertEqual(data.mgf_spectra_short_decoded, list(reader))
            self.assertEqual(data.mgf_spectra_long_decoded, list(func(self.path)))

    def test_read_no_charges(self):
        with mgf.read(self.path, read_charges=False) as reader:
            self.assertEqual(data.mgf_spectra_long_no_charges, list(reader))
        with mgf.read(self.path, False, read_charges=False) as reader:
            self.assertEqual(data.mgf_spectra_short_no_charges, list(reader))

    def test_read_array_conversion(self):
        with mgf.read(self.path, convert_arrays=0) as reader:
            self.assertEqual(data.mgf_spectra_lists, list(reader))
        with mgf.read(self.path, convert_arrays=2) as reader:
            s = next(reader)
            self.assertTrue(isinstance(s['charge array'], np.ma.core.MaskedArray))
            self.assertTrue(isinstance(s['m/z array'], np.ndarray))
        with mgf.read(self.path, convert_arrays=1) as reader:
            s = next(reader)
            self.assertTrue(isinstance(s['charge array'], np.ndarray))
            self.assertTrue(isinstance(s['m/z array'], np.ndarray))

    def test_header(self):
        self.assertEqual(self.header, self.header2)

    def test_readwrite_ns(self):
        self.assertEqual(self.ns, len(self.spectra2))

    def test_readwrite_keys(self):
        for s, s2 in zip(self.spectra, self.spectra2):
            self.assertEqual(set(s), set(s2))
            self.assertEqual(set(s),
                    {'intensity array', 'm/z array', 'params', 'charge array'})

    def test_readwrite_params(self):
        for s, s2 in zip(self.spectra, self.spectra2):
            self.assertEqual(s['params'], s2['params'])

    def test_readwrite_msms_len(self):
        for i in range(self.ns):
            al = len(self.spectra[i]['m/z array'])
            self.assertEqual(al, len(self.spectra[i]['intensity array']))
            self.assertEqual(al, len(self.spectra2[i]['m/z array']))
            self.assertEqual(al, len(self.spectra2[i]['intensity array']))
            for j in range(al):
                self.assertEqual(self.spectra[i]['m/z array'][j],
                        self.spectra2[i]['m/z array'][j])
                self.assertEqual(self.spectra[i]['intensity array'][j],
                        self.spectra2[i]['intensity array'][j])

    def test_readwrite_msms(self):
        for i in range(self.ns):
            al = len(self.spectra[i]['m/z array'])
            for j in range(al):
                self.assertEqual(self.spectra[i]['m/z array'][j],
                        self.spectra2[i]['m/z array'][j])
                self.assertEqual(self.spectra[i]['intensity array'][j],
                        self.spectra2[i]['intensity array'][j])

    def test_read_dtype(self):
        dtypes = {'m/z array': np.float32, 'intensity array': np.int32}
        with mgf.read(self.path, dtype=dtypes) as f:
            for spec in f:
                for k, v in dtypes.items():
                    self.assertEqual(spec[k].dtype, v)

    def test_get_spectrum(self):
        key = 'Spectrum 2'
        for klass in [mgf.MGF, mgf.IndexedMGF]:
            f = klass(self.path)
            self.assertEqual(data.mgf_spectra_long[1], f[key])
            self.assertEqual(data.mgf_spectra_long[1], f.get_spectrum(key))
        self.assertEqual(data.mgf_spectra_long[1], mgf.get_spectrum(self.path, key))

    def test_indexedmgf_picklable(self):
        with mgf.IndexedMGF(self.path) as reader:
            spec = pickle.dumps(reader)
        with pickle.loads(spec) as reader:
            self.assertEqual(data.mgf_spectra_long[0], next(reader))

    def test_map(self):
        with mgf.IndexedMGF(self.path) as reader:
            spectra = sorted(list(reader.map()), key=lambda s: s['params']['title'])
        self.assertEqual(data.mgf_spectra_long, spectra)

if __name__ == "__main__":
    unittest.main()
