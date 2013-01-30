import tempfile
import unittest
from pyteomics.mgf import *
from data import mgf_spectra_long, mgf_spectra_short

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
        self.tmpfile.close()
        self.ns = len(self.spectra)

    def test_read(self):
        # http://stackoverflow.com/q/14246983/1258041
        self.assertEqual(mgf_spectra_long, list(read(self.path)))
        self.assertEqual(mgf_spectra_short, list(read(self.path, False)))
        with read(self.path) as reader:
            self.assertEqual(mgf_spectra_long, list(reader))
        with read(self.path, False) as reader:
            self.assertEqual(mgf_spectra_short, list(reader))

    def test_header(self):
        self.assertEqual(self.header, self.header2)

    def test_readwrite_ns(self):
        self.assertEqual(self.ns, len(self.spectra2))

    def test_readwrite_keys(self):
        for s, s2 in zip(self.spectra, self.spectra2):
            self.assertEqual(set(s.keys()), set(s2.keys()))
            self.assertEqual(set(s.keys()),
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
