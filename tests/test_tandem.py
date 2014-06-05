import unittest
from pyteomics.tandem import *
from data import tandem_spectra

class TandemTest(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None

    def testReadPSM(self):
        for func in [read, chain,
                lambda x, **kw: chain.from_iterable([x], **kw)]:
            with func('test.t.xml') as r:
                self.assertEqual(list(r), tandem_spectra)

if __name__ == '__main__':
    unittest.main()
