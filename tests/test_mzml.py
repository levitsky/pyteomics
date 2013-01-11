import unittest
from pyteomics.mzml import *
from data import mzml_spectra
import numpy as np

class MzmlTest(unittest.TestCase):
    maxDiff = None
    def testReadSpectrum(self):
        # http://stackoverflow.com/q/14246983/1258041
        self.assertEqual(mzml_spectra, list(read('test.mzML')))
            
if __name__ == '__main__':
    unittest.main()
