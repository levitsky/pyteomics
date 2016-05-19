from os import path
import pyteomics
pyteomics.__path__ = [path.abspath(path.join(path.dirname(__file__), path.pardir, 'pyteomics'))]
import tempfile
import unittest
from pyteomics.mgf import *
import data

class MGFTest(unittest.TestCase):
    maxDiff = None
    def setUp(self):
        self.path = 'test.mgf'
        self.header = read_header(self.path)
        self.spectra = list(read(self.path))
        self.tmpfile = tempfile.TemporaryFile(mode='r+')
        write(header=self.header, spectra=self.spectra,
                output=self.tmpfile)
        self.tmpfile.seek(0)
        self.header2 = read_header(self.tmpfile)
        self.tmpfile.seek(0)
        tmpreader = read(self.tmpfile)
        self.spectra2 = list(tmpreader)
        self.ns = len(self.spectra)
        self.tmpfile.close()

    def test_read(self):
        # http://stackoverflow.com/q/14246983/1258041
        self.assertEqual(data.mgf_spectra_long, list(read(self.path)))
        self.assertEqual(data.mgf_spectra_short, list(read(self.path, False)))
        with read(self.path) as reader:
            self.assertEqual(data.mgf_spectra_long, list(reader))
        with read(self.path, False) as reader:
            self.assertEqual(data.mgf_spectra_short, list(reader))

    def test_read_no_charges(self):
        with read(self.path, read_charges=False) as reader:
            self.assertEqual(data.mgf_spectra_long_no_charges, list(reader))
        with read(self.path, False, read_charges=False) as reader:
            self.assertEqual(data.mgf_spectra_short_no_charges, list(reader))

    def test_read_array_conversion(self):
        with read(self.path, convert_arrays=0) as reader:
            self.assertEqual(data.mgf_spectra_lists, list(reader))
        with read(self.path, convert_arrays=2) as reader:
            s = next(reader)
            self.assertTrue(isinstance(s['charge array'], np.ma.core.MaskedArray))
            self.assertTrue(isinstance(s['m/z array'], np.ndarray))
        with read(self.path, convert_arrays=1) as reader:
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

if __name__ == "__main__":
    unittest.main()
