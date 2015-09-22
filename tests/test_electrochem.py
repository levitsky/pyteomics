from os import path
import pyteomics
pyteomics.__path__ = [path.abspath(path.join(path.dirname(__file__), path.pardir, 'pyteomics'))]
import unittest
from pyteomics.electrochem import charge, pI
from pyteomics.auxiliary import PyteomicsError

class ElectrochemTest(unittest.TestCase):
    def setUp(self):
        pass

    def test_charge_calculations_str(self):
        self.assertTrue(
            abs(charge('AAA', 5.0,
                       pK={'H-': [(9., 1)], '-OH': [(8., -1)]},
                       pK_nterm={'H-': {'A': [(3., 1)]}})) < 0.01)
        self.assertTrue(
            abs(charge('H-AAA-OH', 0.0) - 1.0) < 0.01)
        self.assertTrue(
            abs(charge('H-AAA-OH', 14.0) + 1.0) < 0.01)
        self.assertTrue(
            abs(charge('H-AAA-OH', (2.34 + 9.69) / 2.0)) < 0.01)


    def test_charge_calculations_list(self):
        self.assertRaises(PyteomicsError,
            charge, ['A','A','A'], 5.0,
            pK={'H-': [(9., 1)], '-OH': [(8., -1)]},
            pK_nterm={'H-': {'A': [(3., 1)]}})
        self.assertTrue(
            abs(charge(['H-','A','A','A','-OH'], 0.0) - 1.0) < 0.01)
        self.assertTrue(
            abs(charge(['H-','A','A','A','-OH'], 14.0) + 1.0) < 0.01)
        self.assertTrue(
            abs(charge(['H-','A','A','A','-OH'], (2.34 + 9.69) / 2.0)) < 0.01)

    def test_charge_calculations_dict(self):
        self.assertRaises(PyteomicsError, charge, {'H-': 1, '-OH': 1, 'E': 1},
                          7, pK_nterm={'H-': {'A': [(9., 1)]}})
        self.assertTrue(
                abs(charge({'A': 3, 'H-': 1, '-OH': 1}, 14.0) + 1.0) < 0.01)
        self.assertTrue(
                abs(charge({'A': 1, 'H-': 1, '-OH': 1, 'ntermB': 1, 'ctermA': 1},
                    14.0, pK={'H-': [(9., 1)], '-OH': [(8., -1)]},
                    pK_nterm={'H-': {'A': [(3., 1)], 'B': [(3., 1)]}}) + 1.0)
                < 0.01)
        self.assertRaises(PyteomicsError, charge,
                {'A': 1, 'H-': 1, '-OH': 1, 'ctermA': 1}, 14.0,
                    pK={'H-': [(9., 1)], '-OH': [(8., -1)]},
                    pK_nterm={'H-': {'A': [(3., 1)]}})
        self.assertRaises(PyteomicsError, charge,
                {'A': 1, 'H-': 1, '-OH': 1, 'ntermA': 1}, 14.0,
                    pK={'H-': [(9., 1)], '-OH': [(8., -1)]},
                    pK_nterm={'H-': {'A': [(3., 1)]}})
        self.assertRaises(PyteomicsError, charge,
                {'A': 1, 'H-': 1, '-OH': 1, 'ntermA': 2, 'ctermA': 1}, 14.0,
                    pK={'H-': [(9., 1)], '-OH': [(8., -1)]},
                    pK_nterm={'H-': {'A': [(3., 1)]}})
        self.assertRaises(PyteomicsError, charge,
                {'A': 1, 'H-': 1, 'ntermA': 1, 'ctermA': 1}, 14.0,
                    pK={'H-': [(9., 1)], '-OH': [(8., -1)]},
                    pK_nterm={'H-': {'A': [(3., 1)]}})

    def test_pI_calculations(self):
        self.assertTrue(
            abs(pI('H-AAA-OH') - (2.34 + 9.69) / 2.0) < 0.01)

    def test_pI_precision(self):
        pI_best = pI('PEPTIDE', precision_pI=1e-15)
        for i in range(16):
            precision = 10 ** (-i)
            self.assertTrue(
                abs(pI('PEPTIDE', precision_pI=precision) - pI_best) < precision)

    def test_charge_input(self):
        for i in range(0, 14):
            self.assertAlmostEqual(
                charge('H-ACDEFGH-OH', i),
                charge(['H-', 'A', 'C', 'D', 'E', 'F', 'G', 'H', '-OH'], i))
        for i in range(0, 14):
            self.assertAlmostEqual(
                charge('H-ACDEFGH-OH', i),
                charge({'H-': 1, 'A': 1, 'C': 1, 'D': 1,
                        'E': 1, 'F': 1, 'G': 1, 'H': 1, '-OH': 1}, i))




if __name__ == '__main__':
    unittest.main()
