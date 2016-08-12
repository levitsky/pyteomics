from os import path
import pyteomics
pyteomics.__path__ = [path.abspath(path.join(path.dirname(__file__), path.pardir, 'pyteomics'))]

import unittest
from itertools import product
from data import pairs
from pyteomics.openms.trafoxml import *

class TrafoXMLTest(unittest.TestCase):
    maxDiff = None
    def testRead(self):
        for rs, it in product([True, False], repeat=2):
            for func in [TrafoXML, read, chain,
                    lambda x, **kw: chain.from_iterable([x], **kw)]:
                with func('test.trafoXML', read_schema=rs, iterative=it) as r:
                    self.assertEqual(pairs, list(r))

if __name__ == '__main__':
    unittest.main()
