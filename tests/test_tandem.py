import unittest
from pyteomics.tandem import *
from data import tandem_spectra

class PepxmlTest(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None

    def testReadPSM(self):
        psms = list(read('test.t.xml'))
        self.assertEqual(psms, tandem_spectra)

if __name__ == '__main__':
    unittest.main()
