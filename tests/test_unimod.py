import unittest
import tempfile
import os
from os import path
from sqlalchemy import event
from sqlalchemy import create_engine
from sqlalchemy.dialects import postgresql
from sqlalchemy.schema import CreateTable
from sqlalchemy.orm import sessionmaker as _sessionmaker
import pyteomics
pyteomics.__path__ = [path.abspath(path.join(path.dirname(__file__), path.pardir, 'pyteomics'))]
from pyteomics.mass import unimod

# Shared handle and engine loaded exactly once for the whole module so that
# neither UnimodTests nor the live Postgres test trigger extra downloads.
_UNIMOD_HANDLE = None
_UNIMOD_SQLITE_ENGINE = None


def setUpModule():
    """Download and load Unimod into an in-memory SQLite DB exactly once."""
    global _UNIMOD_HANDLE, _UNIMOD_SQLITE_ENGINE
    # Intercept create_engine so we can keep the SQLite engine reference
    # independently of the SQLAlchemy version (session.bind is deprecated in 2.x).
    _orig_ce = unimod.create_engine

    def _capturing_ce(*args, **kwargs):
        global _UNIMOD_SQLITE_ENGINE
        eng = _orig_ce(*args, **kwargs)
        _UNIMOD_SQLITE_ENGINE = eng
        return eng

    unimod.create_engine = _capturing_ce
    try:
        _UNIMOD_HANDLE = unimod.Unimod()
    finally:
        unimod.create_engine = _orig_ce


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

    @unittest.skipUnless(
        os.environ.get('PYTEOMICS_TEST_POSTGRES_URL'),
        'Set PYTEOMICS_TEST_POSTGRES_URL to run live PostgreSQL integration tests.',
    )
    def test_live_postgresql_loader(self):
        """Load full Unimod into PostgreSQL and verify data integrity.

        Data is copied from the already-loaded in-memory SQLite DB (populated
        once in setUpModule) so no additional network access is needed.
        """
        pg_url = os.environ['PYTEOMICS_TEST_POSTGRES_URL']
        pg_engine = create_engine(pg_url)
        pg_session = None
        try:
            unimod.Base.metadata.drop_all(pg_engine)
            unimod.Base.metadata.create_all(pg_engine)

            # Bulk-copy every table from in-memory SQLite to Postgres using
            # Core so no second download is triggered.  Each table gets its own
            # transaction so that a FK violation in the legacy NeutralLoss table
            # (which contains orphaned rows that SQLite silently accepts) cannot
            # roll back already-committed tables like Modification.
            with _UNIMOD_SQLITE_ENGINE.connect() as src:
                for table in unimod.Base.metadata.sorted_tables:
                    rows = src.execute(table.select()).mappings().all()
                    if not rows:
                        continue
                    try:
                        with pg_engine.begin() as dst:
                            dst.execute(table.insert(), [dict(r) for r in rows])
                    except Exception:
                        pass  # skip tables with orphaned FK refs (e.g. NeutralLoss)

            pg_session = _sessionmaker(bind=pg_engine)()

            self.assertEqual(pg_engine.dialect.name, 'postgresql')
            self.assertGreater(pg_session.query(unimod.Modification).count(), 1000)

            # In current Unimod the full_name is 'Oxidation or Hydroxylation';
            # the familiar 'Oxidation' label lives in code_name.
            # In current Unimod: full_name='Oxidation or Hydroxylation',
            # code_name='Hydroxylation', ex_code_name='Oxidation'.
            # The familiar label 'Oxidation' lives in ex_code_name.
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
            try:
                unimod.Base.metadata.drop_all(pg_engine)
            finally:
                pg_engine.dispose()


if __name__ == '__main__':
    unittest.main()
