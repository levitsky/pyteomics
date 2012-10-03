import unittest
from pyteomics.electrochem import charge, pI
from pyteomics.auxiliary import PyteomicsError

class ElectrochemTest(unittest.TestCase):
    def setUp(self):
        pass

    def test_charge_calculations(self):
        self.assertTrue(
            abs(charge('H-AAA-OH', 0.0) - 1.0) < 0.01)
        self.assertTrue(
            abs(charge('H-AAA-OH', 14.0) + 1.0) < 0.01)
        self.assertTrue(
            abs(charge('H-AAA-OH', (2.34 + 9.69) / 2.0)) < 0.01)
        self.assertRaises(PyteomicsError, charge, 'O', 7)
        
    def test_charge_input(self):
        for i in range(0, 14):
            self.assertAlmostEqual(
                charge('H-ACDEFGH-OH', i),
                charge(['H-', 'A', 'C', 'D', 'E', 'F', 'G', 'H', '-OH'], i))
        for i in range(0, 14):
            self.assertAlmostEqual(
                charge('H-ACDEFGH-OH', i),
                charge({'H-':1, 'A':1, 'C':1, 'D':1,
                        'E':1, 'F':1, 'G':1, 'H':1, '-OH':1}, i))

    def test_pI_calculations(self):
        self.assertTrue(
            abs(pI('H-AAA-OH') - (2.34 + 9.69) / 2.0) < 0.01)

    def test_pI_precision(self):
        pI_best = pI('PEPTIDE', precision_pI = 1e-15)
        for i in range(16):
            precision = 10**(-i)
            self.assertTrue(
                abs(pI('PEPTIDE', precision_pI = precision) - pI_best) < precision)

if __name__ == '__main__':
    unittest.main()
                         
