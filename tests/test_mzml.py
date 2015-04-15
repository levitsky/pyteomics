from itertools import product
import unittest
from pyteomics.mzml import *
from data import mzml_spectra
import numpy as np

class MzmlTest(unittest.TestCase):
    maxDiff = None
    def testReadSpectrum(self):
        for rs, it in product([True, False], repeat=2):
            for func in [MzML, read, chain,
                    lambda x, **kw: chain.from_iterable([x], **kw)]:
                with func('test.mzML', read_schema=rs, iterative=it) as r:
                    # http://stackoverflow.com/q/14246983/1258041
                    self.assertEqual(mzml_spectra, list(r))

if __name__ == '__main__':
    unittest.main()
