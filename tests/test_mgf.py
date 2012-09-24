import tempfile
import unittest
from pyteomics.mgf import *

class MGFTest(unittest.TestCase):
    def setUp(self):
        path = 'test.mgf'
        self.header = read_header(path)
        self.spectra = [s for s in read(path)]
        self.tmpfile = tempfile.TemporaryFile(mode='r+')
        write(header=self.header, spectra=self.spectra,
                output=self.tmpfile, close=False)
        self.tmpfile.seek(0)
        self.header2 = read_header(self.tmpfile, close=False)
        self.tmpfile.seek(0)
        tmpreader = read(self.tmpfile)
        self.spectra2 = [x for x in tmpreader]
        self.tmpfile.close()
        self.ns = len(self.spectra)

    def test_header(self):
        self.assertEqual(self.header, self.header2)

    def test_readwrite_ns(self):
        self.assertEqual(self.ns, len(self.spectra2))

    def test_readwrite_keys(self):
        for i in range(self.ns):
            self.assertEqual(set(self.spectra[i].keys()),
                    set(self.spectra2[i].keys()))
            self.assertEqual(set(self.spectra[i].keys()),
                    {'intensities', 'masses', 'params', 'charges'})

    def test_readwrite_params(self):
        for i in range(self.ns):
            self.assertEqual(self.spectra[i]['params'], self.spectra2[i]['params'])

    def test_readwrite_msms_len(self):
        for i in range(self.ns):
            al = len(self.spectra[i]['masses'])
            self.assertEqual(al, len(self.spectra[i]['intensities']))
            self.assertEqual(al, len(self.spectra2[i]['masses']))
            self.assertEqual(al, len(self.spectra2[i]['intensities']))
            for j in range(al):
                self.assertEqual(self.spectra[i]['masses'][j],
                        self.spectra2[i]['masses'][j])
                self.assertEqual(self.spectra[i]['intensities'][j],
                        self.spectra2[i]['intensities'][j])

    def test_readwrite_msms(self):
        for i in range(self.ns):
            al = len(self.spectra[i]['masses'])
            for j in range(al):
                self.assertEqual(self.spectra[i]['masses'][j],
                        self.spectra2[i]['masses'][j])
                self.assertEqual(self.spectra[i]['intensities'][j],
                        self.spectra2[i]['intensities'][j])

if __name__ == "__main__":
    unittest.main()
