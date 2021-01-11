from os import path
import unittest
import pyteomics
pyteomics.__path__ = [path.abspath(path.join(path.dirname(__file__), path.pardir, 'pyteomics'))]
from pyteomics import mztab


class MzTabTest(unittest.TestCase):

    path = 'test.mztab'

    def test_metadata(self):
        reader = mztab.MzTab(self.path)
        self.assertEqual(len(reader.metadata), 208)
        value = reader.metadata['fixed_mod[1]']
        self.assertEqual(value, 'CHEMMOD:57.0214637236')

    def test_iter(self):
        reader = mztab.MzTab(self.path)
        tables = list(reader)
        self.assertEqual(len(tables), 6)
        [self.assertEqual(len(t), 2) for t in tables]

    def test_getitem(self):
        reader = mztab.MzTab(self.path)
        table = reader['psm']
        self.assertIsInstance(table, mztab.pd.DataFrame)

if __name__ == '__main__':
    unittest.main()
