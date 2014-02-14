import unittest
from pyteomics.pepxml import *
from data import pepxml_spectra

class PepxmlTest(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None

    def testReadPSM(self):
        for rs in [True, False]:
            psms = list(read('test.pep.xml', read_schema=rs))
            self.assertEqual(psms, pepxml_spectra)

if __name__ == '__main__':
    unittest.main()
