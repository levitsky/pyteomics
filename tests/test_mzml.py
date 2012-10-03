import unittest
from pyteomics.mzml import *

import numpy as np

class MzmlTest(unittest.TestCase):
    def testReadSpectrum(self):
        mz_array = np.load('test_mzml_mz.npy')
        intensity_array = np.load('test_mzml_intensity.npy')
        for spectrum in read('test.mzML'):
            self.assertTrue(
                all(abs(spectrum['m/z array'] - mz_array) < 1.0e-4))
            self.assertTrue(all(np.equal(
                spectrum['intensity array'], intensity_array)))
            self.assertTrue(
                abs(spectrum['highest observed m/z'] - max(mz_array)) < 1e-4)
            self.assertTrue(
                abs(spectrum['lowest observed m/z'] - min(mz_array)) < 1e-4)
            
if __name__ == '__main__':
    unittest.main()
