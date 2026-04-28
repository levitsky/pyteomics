import unittest
import tempfile
from os import path
from sqlalchemy import event
from sqlalchemy.dialects import postgresql
from sqlalchemy.schema import CreateTable
import pyteomics
pyteomics.__path__ = [path.abspath(path.join(path.dirname(__file__), path.pardir, 'pyteomics'))]
from pyteomics.mass import unimod


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


class UnimodCompatibilityTests(unittest.TestCase):

    def _minimal_unimod_xml(self):
        return '''<?xml version="1.0" encoding="UTF-8"?>
<unimod>
  <elements_row record_id="1" avge_mass="1.0079" mono_mass="1.0078" full_name="Hydrogen" element="H"/>
  <bricks_row record_id="1" brick="H" full_name="Hydrogen"/>
  <brick2element_row record_id="1" brick_key="1" num_element="1" element="H"/>
</unimod>
'''

    def test_dependency_ordering_places_parent_tables_first(self):
        ordered_models = list(unimod._iter_models_in_dependency_order())
        model_pos = {model.__name__: i for i, model in enumerate(ordered_models)}
        self.assertLess(model_pos['Element'], model_pos['BrickToElement'])
        self.assertLess(model_pos['Brick'], model_pos['BrickToElement'])

    def test_loader_succeeds_with_strict_foreign_keys(self):
        xml_path = None
        db_path = None
        db_url = None
        old_create_engine = unimod.create_engine
        old_model_registry = unimod.model_registry
        try:
            with tempfile.NamedTemporaryFile('w', suffix='.xml', delete=False) as fh:
                fh.write(self._minimal_unimod_xml())
                xml_path = fh.name
            db_path = xml_path + '.db'
            db_url = 'sqlite:///' + db_path

            def _create_engine_with_fk(*args, **kwargs):
                engine = old_create_engine(*args, **kwargs)
                if engine.dialect.name == 'sqlite':
                    @event.listens_for(engine, 'connect')
                    def _enable_fk_constraints(dbapi_conn, _):
                        cursor = dbapi_conn.cursor()
                        cursor.execute('PRAGMA foreign_keys=ON')
                        cursor.close()
                return engine

            # Deliberately break registration order to emulate backend-sensitive
            # failures that appear with strict FK enforcement.
            unimod.model_registry = [unimod.BrickToElement, unimod.Brick, unimod.Element]
            unimod.create_engine = _create_engine_with_fk

            handle = unimod.load(xml_path, db_url)
            relation = handle.query(unimod.BrickToElement).one()
            self.assertEqual(relation.brick_id, 1)
            self.assertEqual(relation.element, 'H')
            handle.close()
        finally:
            unimod.create_engine = old_create_engine
            unimod.model_registry = old_model_registry
            if xml_path is not None and path.exists(xml_path):
                import os
                os.remove(xml_path)
            if db_path is not None and path.exists(db_path):
                import os
                os.remove(db_path)

    def test_postgresql_ddl_quotes_reserved_column_names(self):
        sql = str(CreateTable(unimod.Specificity.__table__).compile(dialect=postgresql.dialect()))
        self.assertIn('"group"', sql)


if __name__ == '__main__':
    unittest.main()
