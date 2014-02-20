import unittest
from pyteomics.mzid import *
from data import mzid_spectra
from itertools import product

class MzidTest(unittest.TestCase):
    maxDiff = None
    def testReadPSM(self):
        for rec, refs, rs in product((True, False), repeat=3):
            for func in [read, chain,
                    lambda x, **kw: chain.from_iterable([x], **kw)]:
                with func('test.mzid', recursive=rec, retrieve_refs=refs,
                        read_schema=rs) as reader:
                    psms = list(reader)
                    self.assertEqual(psms, mzid_spectra[(rec, refs)])

if __name__ == '__main__':
    unittest.main()
