import unittest
from pyteomics.tandem import *
from data import tandem_spectra

class TandemTest(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None

    def testReadPSM(self):
        for rs in range(2):
            psms = list(read('test.t.xml', read_schema=rs))
            self.assertEqual(psms, tandem_spectra)

if __name__ == '__main__':
    unittest.main()
