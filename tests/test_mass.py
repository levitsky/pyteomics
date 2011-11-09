import unittest
import random
from pyteomics import mass, auxiliary, parser

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
        self.test_aa_mass = {'X':1.0, 'Y':2.0, 'Z':3.0}
        self.random_peptides = [
            ''.join([random.choice('XYZ') for i in range(20)])
            for i in range(10)]

        self.aa_comp = {'X':   mass.Composition({'A':1},
                                                mass_data=self.mass_data),
                        'Y':   mass.Composition({'B':1},
                                                mass_data=self.mass_data),
                        'Z':   mass.Composition({'C':1},
                                                mass_data=self.mass_data),
                        'H-':  mass.Composition({'D':1},
                                                mass_data=self.mass_data),
                        '-OH': mass.Composition({'E':1},
                                                mass_data=self.mass_data),
                        }

        self.ion_comp = {'M': mass.Composition({},
                                               mass_data=self.mass_data),
                         'a': mass.Composition({'A':-1},
                                               mass_data=self.mass_data)}
        

    def test_fast_mass(self):
        for pep in self.random_peptides:
            self.assertAlmostEqual(
                mass.fast_mass(pep, aa_mass=self.test_aa_mass),
                sum([pep.count(aa) * self.test_aa_mass[aa]
                     for aa in self.test_aa_mass])
                + self.mass_H * 2.0 + self.mass_O )

    def test_Composition(self):
        # Test Composition from a dict.
        self.assertEqual(
            mass.Composition({atom:1 for atom in 'ABCDE'},
                             mass_data=self.mass_data),
            {atom:1 for atom in 'ABCDE'})

        # Test Composition from a formula.
        self.assertEqual(
            mass.Composition(formula='ABCDE',
                             mass_data={atom:{0:(1.0,1.0)}
                                        for atom in 'ABCDE'}),
            {atom:1 for atom in 'ABCDE'})

        # Test Composition from a sequence.        
        self.assertEqual(
            mass.Composition(sequence='XYZ', aa_comp=self.aa_comp),
            {atom:1 for atom in 'ABCDE'})

        # Test Composition from a parsed sequence.
        self.assertEqual(
            mass.Composition(parsed_sequence=['X','Y','Z'],
                             aa_comp=self.aa_comp),
            {atom:1 for atom in 'ABC'})

        # Test sum of Composition objects.
        self.assertEqual(
            mass.Composition(sequence='XXY', aa_comp=self.aa_comp)
            + mass.Composition(sequence='YZZ', aa_comp=self.aa_comp),
            {atom:2 for atom in 'ABCDE'})

        # Test subtraction of Composition objects
        self.assertEqual(
            mass.Composition(formula='', aa_comp=self.aa_comp)
            - mass.Composition(sequence='XYZ', aa_comp=self.aa_comp),
            {atom:-1 for atom in 'ABCDE'})

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
                    parsed_sequence=parser.parse_sequence(
                        pep, labels=['X', 'Y', 'Z'], show_unmodified_termini=True),
                    aa_comp=self.aa_comp, mass_data=self.mass_data,
                    ion_comp=self.ion_comp))

    def test_most_probable_isotopic_composition(self):
        self.assertEqual(
            mass.most_probable_isotopic_composition(
                formula='F',
                mass_data=self.mass_data),
            mass.Composition({'F[6]': 1.0, 'F[7]': 0.0},
                             mass_data=self.mass_data))
        
        self.assertEqual(
            mass.most_probable_isotopic_composition(
                formula='F10',
                mass_data=self.mass_data),
            mass.Composition({'F[6]': 7.0, 'F[7]': 3.0},
                             mass_data=self.mass_data))

        self.assertEqual(
            mass.most_probable_isotopic_composition(
                formula='A20F10',
                elements_with_isotopes = ['F'],
                mass_data=self.mass_data),
            mass.Composition({'A': 20.0, 'F[6]': 7.0, 'F[7]': 3.0},
                             mass_data=self.mass_data))
                        
    def test_isotopic_composition_abundance(self):
        for peplen in range(1,10):
            self.assertEqual(
                mass.isotopic_composition_abundance(formula='F[6]' * peplen,
                                                    mass_data=self.mass_data),
                self.mass_data['F'][6][1] ** peplen)
            
            self.assertEqual(
                mass.isotopic_composition_abundance(formula='AF[6]' * peplen,
                                                    mass_data=self.mass_data),
                self.mass_data['F'][6][1] ** peplen)
            
            self.assertEqual(
                mass.isotopic_composition_abundance(
                    formula='A[1]F[6]' * peplen,
                    mass_data=self.mass_data),
                (self.mass_data['A'][1][1]
                 * self.mass_data['F'][6][1] ) ** peplen)

                                            
if __name__ == '__main__':
    unittest.main()
