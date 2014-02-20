import unittest
from pyteomics.tandem import *
from data import tandem_spectra

class TandemTest(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None

    def testReadPSM(self):
        for func in [read, chain,
                lambda x, **kw: chain.from_iterable([x], **kw)]:
            for rs in range(2):
                with func('test.t.xml', read_schema=rs) as r:
                    self.assertEqual(list(r), tandem_spectra)

if __name__ == '__main__':
    unittest.main()
