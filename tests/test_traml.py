from os import path
import pyteomics
pyteomics.__path__ = [path.abspath(path.join(path.dirname(__file__), path.pardir, 'pyteomics'))]

import unittest
from itertools import product
import operator as op
from data import transitions
from pyteomics.traml import TraML, read, chain

from psims.controlled_vocabulary.controlled_vocabulary import obo_cache
obo_cache.cache_path = '.'
obo_cache.enabled = True


class TraMLTest(unittest.TestCase):
    maxDiff = None
    path = 'ToyExample1.TraML'

    def test_read(self):
        for rs, it, ui, rr in product([True, False], repeat=4):
            for func in [TraML, read, chain, lambda x, **kw: chain.from_iterable([x], **kw)]:
                with func(self.path, read_schema=rs, iterative=it, use_index=ui, retrieve_refs=rr) as r:
                    self.assertEqual(transitions[rr], list(r))

    def test_map(self):
        with TraML(self.path) as r:
            self.assertEqual(
                sorted(transitions[1], key=op.itemgetter('id')),
                sorted(r.map(), key=op.itemgetter('id')))


if __name__ == '__main__':
    unittest.main()
