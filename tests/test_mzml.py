import unittest
from pyteomics.mzml import *

import pylab
import numpy as np

class MzmlTest(unittest.TestCase):
    def setUp(self):
        pass

    def testReadSpectrum(self):
        mz_array = np.load('test_mzml_mz.npy')
        intensity_array = np.load('test_mzml_intensity.npy')
        for spectrum in iter_spectrum('test.mzML'):
            self.assertTrue(
                all(spectrum['m/z array'] - mz_array < 1.0e-4))
            self.assertTrue(all(np.equal(
                spectrum['intensity array'], intensity_array)))

if __name__ == '__main__':
    unittest.main()
