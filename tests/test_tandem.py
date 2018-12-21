from os import path
import pyteomics
pyteomics.__path__ = [path.abspath(path.join(path.dirname(__file__), path.pardir, 'pyteomics'))]
import unittest
from pyteomics import tandem
from data import tandem_spectra


class TandemTest(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None
        self.path = 'test.t.xml'

    def testReadPSM(self):
        for func in [tandem.TandemXML, tandem.read, tandem.chain,
                lambda x, **kw: tandem.chain.from_iterable([x], **kw),
                lambda x, **kw: tandem.filter(x, fdr=1, full_output=False),
                lambda x, **kw: tandem.filter.chain(x, fdr=1, full_output=False),
                lambda x, **kw: tandem.filter.chain.from_iterable([x], fdr=1, full_output=False)]:
            for it in range(2):
                with func(self.path, iterative=it) as r:
                    self.assertEqual(list(r), tandem_spectra)

    def test_df(self):
        df = tandem.DataFrame(self.path)
        self.assertEqual(df.shape, (1, 29))

if __name__ == '__main__':
    unittest.main()
