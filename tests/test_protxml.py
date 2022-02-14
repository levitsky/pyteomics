from os import path
import pyteomics
pyteomics.__path__ = [path.abspath(path.join(path.dirname(__file__), path.pardir, 'pyteomics'))]
from itertools import product
import unittest
from pyteomics import protxml
from data import protxml_results
import operator as op

class ProtXMLTest(unittest.TestCase):
    maxDiff = None
    _kw = {'full_output': False, 'fdr': 1.1,
        'key': op.itemgetter('probability'),
        'reverse': True,
        'remove_decoy': False,
        }
    path = 'test.prot.xml'

    def test_read(self):
        for rs, it in product([True, False], repeat=2):
            for func in [protxml.ProtXML, protxml.read, protxml.chain,
                    lambda x, **kw: protxml.chain.from_iterable([x], **kw),
                    lambda x, **kw: protxml.filter(x, **ProtXMLTest._kw),
                    lambda x, **kw: protxml.filter.chain(x, **ProtXMLTest._kw),
                    lambda x, **kw: protxml.filter.chain.from_iterable([x], **ProtXMLTest._kw)
                    ]:
                with func(self.path, read_schema=rs, iterative=it) as r:
                    self.assertEqual(list(r), protxml_results)

    def test_fdr(self):
        with protxml.ProtXML(self.path) as f:
            self.assertEqual(protxml.fdr(f), 1.0)

    def test_filter(self):
        kw = self._kw.copy()
        kw['remove_decoy'] = True
        x = protxml.filter(self.path, **kw)
        self.assertEqual(list(x), [protxml_results[0]])

    def test_qvalues(self):
        q = protxml.qvalues(self.path, **self._kw)
        self.assertEqual(list(q['q']), [0, 1])

    def test_qvalues_prefix(self):
        q = protxml.qvalues(self.path, decoy_prefix='DECO', **self._kw)
        self.assertEqual(list(q['q']), [0, 1])

    def test_df(self):
        df = protxml.DataFrame(self.path)
        self.assertEqual(df.shape, (2, 15))

    def test_filter_df(self):
        kw = self._kw.copy()
        del kw['full_output']
        del kw['key']
        fdf = protxml.filter_df(self.path, **kw)
        self.assertEqual(fdf.shape, (2, 17))

    def test_filter_df_suffix(self):
        kw = self._kw.copy()
        del kw['full_output']
        del kw['key']
        kw['remove_decoy'] = True
        df = protxml.DataFrame(self.path)
        df['protein_name'] = df.protein_name.str.replace(r'DECOY_(.*)', r'\1_SUF', regex=True)
        fdf = protxml.filter_df(df, decoy_suffix='_SUF', **kw)
        self.assertEqual(fdf.shape, (1, 17))

if __name__ == '__main__':
    unittest.main()
