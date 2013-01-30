import unittest
from pyteomics.mzid import *
from data import mzid_spectra
from itertools import product

class MzidTest(unittest.TestCase):
    def testReadPSM(self):
        for rec, refs in product((True, False), repeat=2):
            with read('test.mzid', recursive=rec, retrieve_refs=refs) as reader:
                psms = list(reader)
                self.assertEqual(psms, mzid_spectra[(rec, refs)])

if __name__ == '__main__':
    unittest.main()
