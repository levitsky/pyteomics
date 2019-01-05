from os import path
import pyteomics
pyteomics.__path__ = [path.abspath(path.join(path.dirname(__file__), path.pardir, 'pyteomics'))]

import unittest
import operator as op
from itertools import product
from data import features
from pyteomics.openms.featurexml import FeatureXML, read, chain

class FeatureXMLTest(unittest.TestCase):
    maxDiff = None
    path = 'test.featureXML'
    def testRead(self):
        for rs, it, ui in product([True, False], repeat=3):
            for func in [FeatureXML, read, chain,
                    lambda x, **kw: chain.from_iterable([x], **kw)]:
                with func(self.path, read_schema=rs, iterative=it, use_index=ui) as r:
                    self.assertEqual(features, list(r))

    def test_map(self):
        self.assertEqual(sorted(features, key=op.itemgetter('id')),
            sorted(FeatureXML(self.path).map(), key=op.itemgetter('id')))

if __name__ == '__main__':
    unittest.main()
