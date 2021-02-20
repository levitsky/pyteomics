from os import path
import unittest
import warnings
import pyteomics
pyteomics.__path__ = [path.abspath(path.join(path.dirname(__file__), path.pardir, 'pyteomics'))]
from pyteomics import mztab


class MzTabTest(unittest.TestCase):

    path_mztab1 = 'test.mztab'
    path_mztab2 = 'test_mztab2.mztab'

    def test_metadata_mztab1(self):
        reader_mztab1 = mztab.MzTab(self.path_mztab1)
        self.assertEqual(len(reader_mztab1.metadata), 208)
        value_from_mztab1 = reader_mztab1.metadata['fixed_mod[1]']
        self.assertEqual(value_from_mztab1, 'CHEMMOD:57.0214637236')

    def test_metadata_mztab2(self):
        reader_mztab2 = mztab.MzTab(self.path_mztab2)
        self.assertEqual(len(reader_mztab2.metadata), 61)
        value_from_mztab2 = reader_mztab2.metadata['sample_processing[1]']
        self.assertEqual(value_from_mztab2, 'high performance liquid chromatography')

    def test_metadata_variant_P(self):
        reader_mztab1 = mztab.MzTab(self.path_mztab1)
        self.assertEqual(reader_mztab1.variant, 'P')

    def test_metadata_variant_M(self):
        reader_mztab2 = mztab.MzTab(self.path_mztab2)
        self.assertEqual(reader_mztab2.variant, 'M')

    def test_iter_mztab1(self):
        reader_mztab1 = mztab.MzTab(self.path_mztab1)
        tables = list(reader_mztab1)
        self.assertEqual(len(tables), 4)
        [self.assertEqual(len(t), 2) for t in tables]

    def test_iter_mztab2(self):
        reader_mztab2 = mztab.MzTab(self.path_mztab2)
        tables = list(reader_mztab2)
        self.assertEqual(len(tables), 3)
        [self.assertEqual(len(t), 2) for t in tables]

    def test_getitem_mztab1(self):
        reader_mztab1 = mztab.MzTab(self.path_mztab1)
        table = reader_mztab1['psm']
        self.assertIsInstance(table, mztab.pd.DataFrame)

    def test_getitem_mztab2(self):
        reader_mztab2 = mztab.MzTab(self.path_mztab2)
        table = reader_mztab2['sme']
        self.assertIsInstance(table, mztab.pd.DataFrame)

    def test_keys_values_items(self):
        reader_mztab2 = mztab.MzTab(self.path_mztab2, table_format='dict')
        keys = list(reader_mztab2.keys())
        self.assertEqual(keys, [k for k, v in reader_mztab2])
        values = list(reader_mztab2.values())
        self.assertEqual(values, [v for k, v in reader_mztab2])
        items = list(reader_mztab2.items())
        self.assertEqual(items, list(reader_mztab2))

    def test_generated_accessors(self):
        reader = mztab.MzTab(self.path_mztab1)
        self.assertEqual(reader.mode, 'Complete')
        self.assertEqual(reader.version, '1.0.0')
        self.assertEqual(reader.software, {1: ('MaxQuant', '1.6.3.4')})
        ms_runs = reader.ms_runs
        self.assertEqual(len(ms_runs), 63)
        self.assertEqual(
            sorted(ms_runs[1].items()),
            [
                ('format', 'Andromeda:apl file format'),
                ('id_format', 'scan number only nativeID format'),
                ('location', 'file://c:/users/jklein/projects/msv000080527_abelin2017/combined/andromeda/allspectra.hcd.ftms.secpep.sil0_0.apl'),
             ])

    def test_missing_version(self):
        class OverridingMzTab(mztab.MzTab):
            def _parse(self):
                super(OverridingMzTab, self)._parse()
                self.metadata.pop("mzTab-version", None)
        with warnings.catch_warnings(record=True) as w:
            reader = OverridingMzTab(self.path_mztab1)
            assert reader.variant == 'P'
            assert reader.version == '1.0.0'
        assert len(w) > 0

    def test_override(self):
        class OverridingMzTab(mztab.MzTab):
            def mode(self):
                return super(OverridingMzTab, self).mode
        reader = OverridingMzTab(self.path_mztab1)
        self.assertEqual(reader.mode(), 'Complete')


if __name__ == '__main__':
    unittest.main()
