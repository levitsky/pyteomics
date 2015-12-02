from os import path
import pyteomics
pyteomics.__path__ = [path.abspath(path.join(path.dirname(__file__), path.pardir, 'pyteomics'))]
from pyteomics.mass import unimod
import unittest


class UnimodTests(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.handle = unimod.Unimod()

    def test_modifications_have_composition(self):
        for modification in self.handle:
            if modification.composition is None:
                assert modification._composition is None
            else:
                assert modification.composition.mass() != 0

    def test_find_modification_by_name(self):
        for modification in ['Acetyl', 'Deamidated', 'GG', 'GlyGly']:
            mod = self.handle.get(modification)
            names = ({mod.full_name, mod.code_name, mod.ex_code_name} |
                     {alt for alt in mod.alternative_names})
            assert modification in names

    def test_composition_parser_sign(self):
        amid = self.handle['Amidated']
        deamid = self.handle['Deamidated']
        assert amid.composition != deamid.composition

    def test_queries(self):
        mods = self.handle.session.query(
            unimod.Modification).join(
            unimod.Specificity).join(
            unimod.Classification).filter(
            unimod.Classification.classification == 'AA substitution').all()
        for mod in mods:
            assert '->' in mod.full_name

    def test_unimod_get(self):
        mod_by_id = self.handle[2]
        mod_by_name = self.handle['Amidation']
        mod_by_alt_name = self.handle['Top-Down sequencing c-type fragment ion']
        mod_by_name_partial = self.handle.get('Amid', strict=False)
        mod_by_alt_name_partial = self.handle.get('c-type fragment ion', strict=False)
        self.assertEqual(mod_by_id, mod_by_name)
        self.assertEqual(mod_by_id, mod_by_alt_name)
        self.assertEqual(mod_by_id, mod_by_name_partial)
        self.assertEqual(mod_by_id, mod_by_alt_name_partial)

    def test_raise_on_nonmatching_id(self):
        self.assertRaises(KeyError, self.handle.get, -1)
        self.assertRaises(KeyError, self.handle.get, 'NotAModification')

    def test_constructor(self):
        import os
        if os.path.exists('unimod.db'):
            os.remove('unimod.db')
        for name in [None, 'sqlite:///unimod.db']:
            handle = unimod.Unimod(name)
            self.assertEqual(self.handle[1], handle[1])

if __name__ == '__main__':
    unittest.main()
