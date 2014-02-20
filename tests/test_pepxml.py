import unittest
from pyteomics.pepxml import *
from data import pepxml_spectra

class PepxmlTest(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None

    def testReadPSM(self):
        for rs in [True, False]:
            for func in [read, chain,
                    lambda x, **kw: chain.from_iterable([x], **kw)]:
                with func('test.pep.xml', read_schema=rs) as r:
                    self.assertEqual(list(r), pepxml_spectra)

if __name__ == '__main__':
    unittest.main()
