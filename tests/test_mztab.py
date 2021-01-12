from os import path
import unittest
import pyteomics
pyteomics.__path__ = [path.abspath(path.join(path.dirname(__file__), path.pardir, 'pyteomics'))]
from pyteomics import mztab


class MzTabTest(unittest.TestCase):

    path_mztab1 = 'test.mztab'
    path_mztab2 = 'test_mztab2.mztab'

    def test_metadata(self):
        reader_mztab1 = mztab.MzTab(self.path_mztab1)
        self.assertEqual(len(reader_mztab1.metadata), 208)
        value_from_mztab1 = reader_mztab1.metadata['fixed_mod[1]']
        self.assertEqual(value_from_mztab1, 'CHEMMOD:57.0214637236')

        reader_mztab2 = mztab.MzTab(self.path_mztab2)
        self.assertEqual(len(reader_mztab2.metadata), 61)
        value_from_mztab2 = reader_mztab2.metadata['sample_processing[1]']
        self.assertEqual(value_from_mztab2, 'MSIO, MSIO:0000148, high performance liquid chromatography')


    def test_iter(self):
        reader_mztab1 = mztab.MzTab(self.path_mztab1)
        tables = list(reader_mztab1)
        self.assertEqual(len(tables), 4)
        [self.assertEqual(len(t), 2) for t in tables]

        reader_mztab2 = mztab.MzTab(self.path_mztab2)
        tables = list(reader_mztab2)
        self.assertEqual(len(tables), 3)
        [self.assertEqual(len(t), 2) for t in tables]

    def test_getitem(self):
        reader_mztab1 = mztab.MzTab(self.path_mztab1)
        table = reader_mztab1['psm']
        self.assertIsInstance(table, mztab.pd.DataFrame)

        reader_mztab2 = mztab.MzTab(self.path_mztab2)
        table = reader_mztab2['sme']
        self.assertIsInstance(table, mztab.pd.DataFrame)

if __name__ == '__main__':
    unittest.main()
