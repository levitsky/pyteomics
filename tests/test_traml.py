from os import path
import pyteomics
pyteomics.__path__ = [path.abspath(path.join(path.dirname(__file__), path.pardir, 'pyteomics'))]

import unittest
from itertools import product
from data import transitions
from pyteomics.traml import TraML, read, chain

class FeatureXMLTest(unittest.TestCase):
    maxDiff = None
    path = 'ToyExample1.TraML'
    def testRead(self):
        for rs, it, ui in product([True, False], repeat=3):
            for func in [TraML, read, chain,
                    lambda x, **kw: chain.from_iterable([x], **kw)]:
                with func(self.path, read_schema=rs, iterative=it, use_index=ui) as r:
                    self.assertEqual(transitions, list(r))

if __name__ == '__main__':
    unittest.main()
