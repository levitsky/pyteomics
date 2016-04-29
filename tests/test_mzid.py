from os import path
import pyteomics
pyteomics.__path__ = [path.abspath(path.join(path.dirname(__file__), path.pardir, 'pyteomics'))]
import unittest
from pyteomics.mzid import *
from data import mzid_spectra
from itertools import product

class MzidTest(unittest.TestCase):
    maxDiff = None
    def testReadPSM(self):
        for rec, refs, rs, it, ui in product((True, False), repeat=5):
            for func in [MzIdentML, read, chain,
                    lambda x, **kw: chain.from_iterable([x], **kw)]:
                with func('test.mzid', recursive=rec, retrieve_refs=refs,
                        read_schema=rs, iterative=it, use_index=ui) as reader:
                    psms = list(reader)
                    self.assertEqual(psms, mzid_spectra[(rec, refs)])

if __name__ == '__main__':
    unittest.main()
