from os import path
import pyteomics
pyteomics.__path__ = [path.abspath(path.join(path.dirname(__file__), path.pardir, 'pyteomics'))]
from itertools import product
import unittest
from pyteomics.mzxml import *
from data import mzxml_spectra
import numpy as np

class MzXMLTest(unittest.TestCase):
    maxDiff = None
    path = 'test.mzXML'

    def testReadSpectrum(self):
        for rs, it, ui in product([True, False], repeat=3):
            for func in [MzXML, read, chain,
                    lambda x, **kw: chain.from_iterable([x], **kw)]:
                with func(self.path, read_schema=rs, iterative=it, use_index=ui) as r:
                    # http://stackoverflow.com/q/14246983/1258041
                    self.assertEqual(mzxml_spectra, list(r))

if __name__ == '__main__':
    unittest.main()
