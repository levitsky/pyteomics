import unittest
from pyteomics.mzid import *
from data import mzid_spectra

class MzidTest(unittest.TestCase):
    def testReadPSM(self):
        with read('test.mzid') as reader:
            psms = list(reader)[0:4]
            self.assertEqual(psms, mzid_spectra)

if __name__ == '__main__':
    unittest.main()
