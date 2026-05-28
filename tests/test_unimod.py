import unittest
import tempfile
import os
import warnings
from os import path
from sqlalchemy import event
from sqlalchemy import create_engine
from sqlalchemy.dialects import postgresql
from sqlalchemy.schema import CreateTable
import pyteomics
pyteomics.__path__ = [path.abspath(path.join(path.dirname(__file__), path.pardir, 'pyteomics'))]
from pyteomics.mass import unimod

# Shared handle and engine loaded exactly once for the whole module so that
# neither UnimodTests nor the live Postgres test trigger extra downloads.
_UNIMOD_HANDLE = None
_UNIMOD_XML_PATH = None


def setUpModule():
    """Decompress the bundled Unimod fixture and load Unimod from it."""
    global _UNIMOD_HANDLE, _UNIMOD_XML_PATH
    import gzip
    fixture = path.join(path.dirname(path.abspath(__file__)), 'unimod_tables.xml.gz')
    with tempfile.NamedTemporaryFile('wb', suffix='.xml', delete=False) as fh:
        with gzip.open(fixture, 'rb') as gz:
            fh.write(gz.read())
        _UNIMOD_XML_PATH = fh.name

    old_url = unimod._unimod_xml_download_url
    unimod._unimod_xml_download_url = _UNIMOD_XML_PATH
    try:
        _UNIMOD_HANDLE = unimod.Unimod()
    finally:
        unimod._unimod_xml_download_url = old_url


def tearDownModule():
    global _UNIMOD_XML_PATH
    if _UNIMOD_XML_PATH is not None and path.exists(_UNIMOD_XML_PATH):
        os.remove(_UNIMOD_XML_PATH)
        _UNIMOD_XML_PATH = None


class UnimodTests(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.handle = _UNIMOD_HANDLE

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
        if os.path.exists('unimod.db'):
            os.remove('unimod.db')
        # None → in-memory; 'sqlite:///unimod.db' → file-backed.
        # Both should produce data identical to the module-level shared handle.
        for name in [None, 'sqlite:///unimod.db']:
            handle = unimod.Unimod(name)
            self.assertEqual(_UNIMOD_HANDLE[1], handle[1])


class UnimodCompatibilityTests(unittest.TestCase):

    OXIDATION_MONOISOTOPIC_MASS = 15.994915
    OXIDATION_AVERAGE_MASS = 15.9994

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

    def test_loader_skips_row_errors_and_warns(self):
        xml_path = None
        db_path = None
        db_url = None
        old_create_engine = unimod.create_engine
        try:
            with tempfile.NamedTemporaryFile('w', suffix='.xml', delete=False) as fh:
                fh.write('''<?xml version="1.0" encoding="UTF-8"?>
<unimod>
  <elements_row record_id="1" avge_mass="1.0079" mono_mass="1.0078" full_name="Hydrogen" element="H"/>
  <bricks_row record_id="1" brick="H" full_name="Hydrogen"/>
  <brick2element_row record_id="1" brick_key="1" num_element="1" element="H"/>
  <brick2element_row record_id="2" brick_key="999" num_element="1" element="H"/>
</unimod>
''')
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

            unimod.create_engine = _create_engine_with_fk

            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter('always')
                handle = unimod.load(xml_path, db_url)
                self.assertEqual(handle.query(unimod.BrickToElement).count(), 1)
                handle.close()

            warning_texts = [str(w.message) for w in caught]
            self.assertTrue(any('Skipped' in text for text in warning_texts))
        finally:
            unimod.create_engine = old_create_engine
            if xml_path is not None and path.exists(xml_path):
                import os
                os.remove(xml_path)
            if db_path is not None and path.exists(db_path):
                import os
                os.remove(db_path)

    @unittest.skipUnless(
        os.environ.get('PYTEOMICS_TEST_POSTGRES_URL'),
        'Set PYTEOMICS_TEST_POSTGRES_URL to run live PostgreSQL integration tests.',
    )
    def test_live_postgresql_loader(self):
        """Load full Unimod into PostgreSQL and verify data integrity.

        Data is loaded from the shared XML file downloaded once in setUpModule,
        so no additional network access is needed.
        """
        pg_url = os.environ['PYTEOMICS_TEST_POSTGRES_URL']
        pg_engine = create_engine(pg_url)
        pg_session = None
        try:
            unimod.Base.metadata.drop_all(pg_engine)

            with warnings.catch_warnings(record=True):
                warnings.simplefilter('always')
                pg_session = unimod.load(_UNIMOD_XML_PATH, pg_url)

            self.assertEqual(pg_engine.dialect.name, 'postgresql')
            self.assertGreater(pg_session.query(unimod.Modification).count(), 1000)

            oxidation = pg_session.query(unimod.Modification).filter(
                unimod.Modification.ex_code_name == 'Oxidation').one()
            self.assertAlmostEqual(
                float(oxidation.monoisotopic_mass),
                self.OXIDATION_MONOISOTOPIC_MASS, places=6)
            self.assertAlmostEqual(
                float(oxidation.average_mass),
                self.OXIDATION_AVERAGE_MASS, places=4)
            self.assertEqual(dict(sorted(oxidation.composition.items())), {'O': 1})
        finally:
            if pg_session is not None:
                pg_session.close()
                if pg_session.bind is not None:
                    pg_session.bind.dispose()
            try:
                unimod.Base.metadata.drop_all(pg_engine)
            finally:
                pg_engine.dispose()


if __name__ == '__main__':
    unittest.main()
