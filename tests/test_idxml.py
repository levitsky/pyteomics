from os import path
import pyteomics
pyteomics.__path__ = [path.abspath(path.join(path.dirname(__file__), path.pardir, 'pyteomics'))]
import unittest
from pyteomics.openms.idxml import IDXML, read, chain
from data import idxml_data
from itertools import product

class IdxmlTest(unittest.TestCase):
    maxDiff = None
    path = 'test.idXML'
    def testReadPSM(self):
        for rec, refs, rs, it, ui in product((True, False), repeat=5):
            for func in [IDXML, read, chain,
                    lambda x, **kw: chain.from_iterable([x], **kw)]:
                with func(self.path, recursive=rec, retrieve_refs=refs,
                        read_schema=rs, iterative=it, use_index=ui) as reader:
                    try:
                        psms = list(reader)
                        self.assertEqual(psms, idxml_data[(rec, refs)])
                    except Exception:
                        print('Parameters causing exception: ', rec, refs, rs, it, ui)
                        raise


if __name__ == '__main__':
    unittest.main()
