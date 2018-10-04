from os import path
import unittest
import pyteomics
pyteomics.__path__ = [path.abspath(path.join(path.dirname(__file__), path.pardir, 'pyteomics'))]
from pyteomics import peff


class PEFFTest(unittest.TestCase):
    maxDiff = None

    def setUp(self):
        self.peff_file = 'test.peff'

    def test_parse(self, reader=None):
        if reader is None:
            reader = peff.IndexedPEFF(self.peff_file)
        self.assertEqual(reader.number_of_entries, 5)
        self.assertEqual(len(reader.header_blocks), 1)

        protein = next(reader)
        self.assertEqual(protein.description.Tag, "NX_P07585-1")
        self.assertEqual(protein, reader.get_entry("NX_P07585-1"))

        protein2 = reader.get_entry("NX_P07585-3")
        self.assertEqual(protein, protein)
        self.assertNotEqual(protein, protein2)

        self.assertEqual(protein.description.TaxName, "Homo Sapiens")
        self.assertEqual(protein.description["NcbiTaxId"], 9606)
        self.assertEqual(len(protein.description.ModResPsi), 2)

if __name__ == '__main__':
    unittest.main()
