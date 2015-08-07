from os import path
import pyteomics
pyteomics.__path__ = [path.abspath(path.join(path.dirname(__file__), path.pardir, 'pyteomics'))]
from itertools import product
import unittest
from pyteomics.pepxml import *
from data import pepxml_spectra

class PepxmlTest(unittest.TestCase):
    maxDiff = None
    _kw = {'full_output': False, 'fdr': 1,
        'key': lambda x: min(
            sh['search_score'].get('expect', 1)
            for sh in x['search_hit'])
        }
    def testReadPSM(self):
        for rs, it in product([True, False], repeat=2):
            for func in [PepXML, read, chain,
                    lambda x, **kw: chain.from_iterable([x], **kw),
                    lambda x, **kw: filter(x, **PepxmlTest._kw),
                    lambda x, **kw: filter.chain(x, **PepxmlTest._kw),
                    lambda x, **kw: filter.chain.from_iterable([x], **PepxmlTest._kw)]:
                with func('test.pep.xml', read_schema=rs, iterative=it) as r:
                    self.assertEqual(list(r), pepxml_spectra)

if __name__ == '__main__':
    unittest.main()
