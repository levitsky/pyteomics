import unittest
from pyteomics.mzid import *
from data import mzid_spectra
from itertools import product

class MzidTest(unittest.TestCase):
    maxDiff = None
    def testReadPSM(self):
        for rec, refs, rs, it in product((True, False), repeat=4):
            for func in [MzIdentML, read, chain,
                    lambda x, **kw: chain.from_iterable([x], **kw)]:
                with func('test.mzid', recursive=rec, retrieve_refs=refs,
                        read_schema=rs, iterative=it) as reader:
                    psms = list(reader)
                    self.assertEqual(psms, mzid_spectra[(rec, refs)])

if __name__ == '__main__':
    unittest.main()
