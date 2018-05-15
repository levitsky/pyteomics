from os import path
import pyteomics
pyteomics.__path__ = [path.abspath(path.join(path.dirname(__file__), path.pardir, 'pyteomics'))]
from itertools import product
import unittest
from pyteomics.protxml import ProtXML, read, chain#, filter
from data import protxml_results
import operator as op

class ProtXMLTest(unittest.TestCase):
    maxDiff = None
    _kw = {'full_output': False, 'fdr': 1,
        'key': op.itemgetter('probability'),
        'reverse': True
        }
    path = 'test.prot.xml'

    def test_read(self):
        for rs, it in product([True, False], repeat=2):
            for func in [ProtXML, read, chain,
                    lambda x, **kw: chain.from_iterable([x], **kw),
                    # lambda x, **kw: filter(x, **ProtXMLTest._kw),
                    # lambda x, **kw: filter.chain(x, **ProtXMLTest._kw),
                    # lambda x, **kw: filter.chain.from_iterable([x], **ProtXMLTest._kw)
                    ]:
                with func(self.path, read_schema=rs, iterative=it) as r:
                    self.assertEqual(list(r), protxml_results)

if __name__ == '__main__':
    unittest.main()
