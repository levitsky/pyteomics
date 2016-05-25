from os import path
import pyteomics
pyteomics.__path__ = [path.abspath(path.join(path.dirname(__file__), path.pardir, 'pyteomics'))]

import unittest
from itertools import product
from data import features
from pyteomics.featurexml import *

class FeatureXMLTest(unittest.TestCase):
    maxDiff = None
    def testRead(self):
        for rs, it, ui in product([True, False], repeat=3):
            for func in [FeatureXML, read, chain,
                    lambda x, **kw: chain.from_iterable([x], **kw)]:
                with func('test.featureXML', read_schema=rs, iterative=it, use_index=ui) as r:
                    self.assertEqual(features, list(r))

if __name__ == '__main__':
    unittest.main()
