from os import path
import pyteomics
pyteomics.__path__ = [path.abspath(path.join(path.dirname(__file__), path.pardir, 'pyteomics'))]
import unittest
import random
from pyteomics import mass, auxiliary, parser
import gzip

class MassTest(unittest.TestCase):
    def setUp(self):
        self.mass_data = {
            'A' : {0: (1.0, 1.0),
                   1: (1.0, 0.5),
                   2: (2.0, 0.5)},
            'B' : {0: (2.0, 1.0),
                   2: (2.0, 0.5),
                   3: (3.0, 0.5)},
            'C' : {0: (3.0, 1.0),
                   3: (3.0, 0.5),
                   4: (4.0, 0.5)},
            'D' : {0: (4.0, 1.0),
                   4: (4.0, 0.5),
                   5: (5.0, 0.5)},
            'E' : {0: (5.0, 1.0),
                   5: (5.0, 0.5),
                   6: (6.0, 0.5)},
            'F' : {0: (6.0, 1.0),
                   6: (6.0, 0.7),
                   7: (7.0, 0.3)},
            'H+': {0: (5.0, 1.0),
                   5: (5.0, 1.0)},
            }

        self.mass_H = mass.nist_mass['H'][0][0]
        self.mass_O = mass.nist_mass['O'][0][0]
        self.test_aa_mass = {'X': 1.0, 'Y': 2.0, 'Z': 3.0}
        self.random_peptides = [
            ''.join([random.choice('XYZ') for i in range(20)])
            for i in range(10)]

        self.aa_comp = {'X':   mass.Composition({'A': 1},
                                                mass_data=self.mass_data),
                        'Y':   mass.Composition({'B': 1},
                                                mass_data=self.mass_data),
                        'Z':   mass.Composition({'C': 1},
                                                mass_data=self.mass_data),
                        'H-':  mass.Composition({'D': 1},
                                                mass_data=self.mass_data),
                        '-OH': mass.Composition({'E': 1},
                                                mass_data=self.mass_data),
                        }

        self.ion_comp = {'M': mass.Composition({},
                                               mass_data=self.mass_data),
                         'a': mass.Composition({'A': -1},
                                               mass_data=self.mass_data)}
        self.mods = {'a': mass.Composition(A=1),
                'b': mass.Composition(B=1)}
        self.d = {atom: 1 for atom in 'ABCDE'}

    def test_fast_mass(self):
        for pep in self.random_peptides:
            self.assertAlmostEqual(
                mass.fast_mass(pep, aa_mass=self.test_aa_mass),
                sum(pep.count(aa) * m
                     for aa, m in self.test_aa_mass.items())
                + self.mass_H * 2.0 + self.mass_O)

    def test_fast_mass2(self):
        for pep in self.random_peptides:
            self.assertAlmostEqual(
                mass.fast_mass2(pep, aa_mass=self.test_aa_mass),
                sum(pep.count(aa) * m
                     for aa, m in self.test_aa_mass.items())
                + self.mass_H * 2.0 + self.mass_O)


    def test_Composition_dict(self):
        # Test Composition from a dict.
        self.assertEqual(
            mass.Composition(self.d, mass_data=self.mass_data), self.d)

    def test_Composition_formula(self):
        # Test Composition from a formula.
        self.assertEqual(self.d,
            mass.Composition(formula='ABCDE',
                         mass_data={atom: {0: (1.0, 1.0)} for atom in 'ABCDE'}))

    def test_Composition_seq(self):
        # Test Composition from a sequence.
        self.assertEqual(self.d,
            mass.Composition(sequence='XYZ', aa_comp=self.aa_comp))

    def test_Composition_pseq(self):
        # Test Composition from a parsed sequence.
        self.assertEqual(
            mass.Composition(parsed_sequence=['X', 'Y', 'Z'],
                             aa_comp=self.aa_comp),
            {atom: 1 for atom in 'ABC'})

    def test_Composition_sseq(self):
        # Test Composition from a split sequence.
        self.assertEqual(
            mass.Composition(split_sequence=[('X',), ('Y',), ('Z',)],
                             aa_comp=self.aa_comp),
            {atom: 1 for atom in 'ABC'})

    def test_Composition_sum(self):
        # Test sum of Composition objects.
        self.assertEqual(
            mass.Composition(sequence='XXY', aa_comp=self.aa_comp)
            + mass.Composition(sequence='YZZ', aa_comp=self.aa_comp),
            {atom: 2 for atom in 'ABCDE'})

    def test_Composition_sub(self):
        # Test subtraction of Composition objects
        self.assertEqual({}
            - mass.Composition(sequence='XYZ', aa_comp=self.aa_comp),
            {atom: -1 for atom in 'ABCDE'})

    def test_Composition_mul(self):
        # Test multiplication of Composition by integers
        self.assertEqual(
                2 * mass.Composition(sequence='XYZ', aa_comp=self.aa_comp),
                {atom: 2 for atom in 'ABCDE'})
        self.assertEqual(
                mass.Composition(sequence='XYZ', aa_comp=self.aa_comp) * 2,
                {atom: 2 for atom in 'ABCDE'})

    def test_Composition_positional(self):
        # Test creation from positional args
        ac = self.aa_comp.copy()
        ac.update(self.mods)
        self.assertEqual(mass.Composition('aXbYZ', aa_comp=ac),
                {'A': 2, 'B': 2, 'C': 1, 'D': 1, 'E': 1})
        self.assertEqual(mass.Composition('AB2C3', mass_data=self.mass_data),
                {'A': 1, 'B': 2, 'C': 3})

    def test_calculate_mass(self):
        # Calculate mass by a formula.
        self.assertEqual(
            mass.calculate_mass(formula='ABCDE', mass_data=self.mass_data),
            sum([self.mass_data[atom][0][0] for atom in 'ABCDE']))

        # Calculate mass by a sequence.
        self.assertEqual(
            mass.calculate_mass(sequence='XYZ',
                                aa_comp=self.aa_comp,
                                mass_data=self.mass_data),
            sum([self.mass_data[atom][0][0] for atom in 'ABCDE']))

        # Calculate mass by a parsed sequence.
        self.assertEqual(
            mass.calculate_mass(parsed_sequence=['H-','X','Y','Z','-OH'],
                                aa_comp=self.aa_comp,
                                mass_data=self.mass_data),
            sum([self.mass_data[atom][0][0] for atom in 'ABCDE']))

        # Calculate average mass by a formula.
        self.assertEqual(
            mass.calculate_mass(formula='ABCDE',
                                average=True,
                                mass_data=self.mass_data),
            sum([self.mass_data[atom][isotope][0]
                 * self.mass_data[atom][isotope][1]
                 for atom in 'ABCDE'
                 for isotope in self.mass_data[atom] if isotope != 0]))

        # Calculate m/z of an ion.
        for charge in [1,2,3]:
            self.assertEqual(
                mass.calculate_mass(formula='ABCDE',
                                    ion_type='M',
                                    charge=charge,
                                    mass_data=self.mass_data),
                mass.calculate_mass(formula='ABCDE'+'H+%d' % (charge,),
                                    mass_data=self.mass_data))

            self.assertEqual(
                mass.calculate_mass(formula='ABCDE',
                                    ion_type='M',
                                    charge=charge,
                                    mass_data=self.mass_data),
                (mass.calculate_mass(formula='ABCDE',
                                    mass_data=self.mass_data)
                 + self.mass_data['H+'][0][0] * charge
                 ) / charge)

            self.assertRaises(
                auxiliary.PyteomicsError,
                mass.calculate_mass,
                **{'formula': 'ABCDEH+%d' % charge,
                   'ion_type': 'M',
                   'charge': charge,
                   'mass_data': self.mass_data})

        # Sanity check.
        for pep in self.random_peptides:
            self.assertEqual(mass.calculate_mass(
                sequence=pep, aa_comp=self.aa_comp, mass_data=self.mass_data,
                ion_comp=self.ion_comp),
                mass.calculate_mass(
                    parsed_sequence=parser.parse(
                        pep, labels=['X', 'Y', 'Z'], show_unmodified_termini=True),
                    aa_comp=self.aa_comp, mass_data=self.mass_data,
                    ion_comp=self.ion_comp))

    def test_most_probable_isotopic_composition(self):
        self.assertEqual(
            mass.most_probable_isotopic_composition(
                formula='F',
                mass_data=self.mass_data),
            (mass.Composition({'F[6]': 1, 'F[7]': 0},
                             mass_data=self.mass_data), 0.7))

        self.assertEqual(
            mass.most_probable_isotopic_composition(
                formula='F10',
                mass_data=self.mass_data),
            (mass.Composition({'F[6]': 7, 'F[7]': 3},
                             mass_data=self.mass_data),
            (0.3)**3 * (0.7)**7 * 120))

        self.assertEqual(
            mass.most_probable_isotopic_composition(
                formula='A20F10',
                elements_with_isotopes = ['F'],
                mass_data=self.mass_data),
            (mass.Composition({'A': 20, 'F[6]': 7, 'F[7]': 3},
                             mass_data=self.mass_data),
            (0.3)**3 * (0.7)**7 * 120))

    def test_isotopic_composition_abundance(self):
        for peplen in range(1,10):
            self.assertAlmostEqual(
                mass.isotopic_composition_abundance(formula='F[6]' * peplen,
                                                    mass_data=self.mass_data),
                self.mass_data['F'][6][1] ** peplen)

            self.assertAlmostEqual(
                mass.isotopic_composition_abundance(formula='AF[6]' * peplen,
                                                    mass_data=self.mass_data),
                self.mass_data['F'][6][1] ** peplen)

            self.assertAlmostEqual(
                mass.isotopic_composition_abundance(
                    formula='A[1]F[6]' * peplen,
                    mass_data=self.mass_data),
                (self.mass_data['A'][1][1]
                 * self.mass_data['F'][6][1] ) ** peplen)

    def test_Unimod(self):
        db = mass.Unimod(gzip.open('unimod.xml.gz'))
        for x in db.mods:
            self.assertGreater(0.00001,
                abs(x['mono_mass'] - mass.calculate_mass(x['composition'],
                    mass_data=db.mass_data)))

    def test_nist_mass(self):
        self.assertTrue(all(abs(g[0][1]-1) < 1e-6 for g in mass.nist_mass.values()))
        for g in mass.nist_mass.values():
            s = sum(p[1] for num, p in g.items() if num)
            self.assertTrue(abs(s-1) < 1e-6 or abs(s) < 1e-6)

    def test_composition_objects_are_pickleable(self):
        dict_ = mass.Composition(self.d, mass_data=self.mass_data)
        formula = mass.Composition(formula='ABCDE',
                         mass_data={atom: {0: (1.0, 1.0)} for atom in 'ABCDE'})
        sequence = mass.Composition(sequence='XYZ', aa_comp=self.aa_comp)
        parsed_sequence = mass.Composition(parsed_sequence=['X', 'Y', 'Z'],
                             aa_comp=self.aa_comp)
        split_sequence = mass.Composition(split_sequence=[('X',), ('Y',), ('Z',)],
                             aa_comp=self.aa_comp)

        import pickle
        self.assertEqual(dict_, pickle.loads(pickle.dumps(dict_)))
        self.assertEqual(formula, pickle.loads(pickle.dumps(formula)))
        self.assertEqual(sequence, pickle.loads(pickle.dumps(sequence)))
        self.assertEqual(parsed_sequence, pickle.loads(pickle.dumps(parsed_sequence)))
        self.assertEqual(split_sequence, pickle.loads(pickle.dumps(split_sequence)))


if __name__ == '__main__':
    unittest.main()
