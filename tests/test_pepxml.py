import unittest
from pyteomics.pepxml import *
from data import pepxml_spectra

class PepxmlTest(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None

    def testReadPSM(self):
        psms = list(read('test.pep.xml'))
        self.assertEqual(psms, pepxml_spectra)

if __name__ == '__main__':
    unittest.main()
