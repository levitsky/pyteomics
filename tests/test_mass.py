from os import path
import pyteomics
pyteomics.__path__ = [path.abspath(path.join(path.dirname(__file__), path.pardir, 'pyteomics'))]
import unittest
import random
import pickle
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

        self.aa_comp = {
            'X':   mass.Composition({'A': 1}),
            'Y':   mass.Composition({'B': 1}),
            'Z':   mass.Composition({'C': 1}),
            'F':   mass.Composition({'F': 1}),
            'H-':  mass.Composition({'D': 1}),
            '-OH': mass.Composition({'E': 1}),
        }

        self.ion_comp = {
            'M': mass.Composition({}),
            'a': mass.Composition({'A': -1})
        }

        self.mods = {'a': mass.Composition(A=1), 'b': mass.Composition(B=1)}
        self.d = {atom: 1 for atom in 'ABCDE'}

    def test_fast_mass(self):
        for pep in self.random_peptides:
            self.assertAlmostEqual(
                mass.fast_mass(pep, aa_mass=self.test_aa_mass),
                sum(pep.count(aa) * m for aa, m in self.test_aa_mass.items()) + self.mass_H * 2.0 + self.mass_O)

    def test_fast_mass2_no_term(self):
        for pep in self.random_peptides:
            self.assertAlmostEqual(
                mass.fast_mass2(pep, aa_mass=self.test_aa_mass),
                sum(pep.count(aa) * m for aa, m in self.test_aa_mass.items()) + self.mass_H * 2.0 + self.mass_O)

    def test_fast_mass2_sanity(self):
        self.assertAlmostEqual(mass.fast_mass2('PEPTIDE'), mass.fast_mass('PEPTIDE'))
        self.assertAlmostEqual(mass.fast_mass2('PEPTIDE'), 799.359964)

    def test_fast_mass2_term(self):
        for pep in self.random_peptides:
            nterm = 'AB2C3-'
            cterm = '-DE2F3'
            self.assertAlmostEqual(
                mass.fast_mass2(nterm + pep + cterm, aa_mass=self.test_aa_mass, mass_data=self.mass_data),
                sum(pep.count(aa) * m for aa, m in self.test_aa_mass.items()) + (
                    self.mass_data['A'][0][0] + self.mass_data['B'][0][0] * 2 + self.mass_data['C'][0][0] * 3 +
                    self.mass_data['D'][0][0] + self.mass_data['E'][0][0] * 2 + self.mass_data['F'][0][0] * 3))

    def test_fast_mass2_term_label(self):
        mass_data = dict(self.mass_data)
        mass_data['H'] = {0: (self.mass_H, 1.0)}
        mass_data['O'] = {0: (self.mass_O, 1.0)}
        aa_mass = self.test_aa_mass.copy()
        aa_mass.update({k: mass.calculate_mass(composition=v, mass_data=mass_data) for k, v in self.mods.items()})
        for pep in self.random_peptides:
            for mlabel, mcomp in self.mods.items():
                mpep = mlabel + '-' + pep + '-' + mlabel
                self.assertRaises(auxiliary.PyteomicsError,
                    mass.fast_mass2, mpep, mass_data=mass_data, aa_mass=aa_mass)

    def test_composition_term(self):
        aa_comp = self.aa_comp.copy()
        aa_comp.update(self.mods)
        for pep in self.random_peptides:
            for mlabel, mcomp in self.mods.items():
                mpep = mlabel + '-' + pep + '-' + mlabel
                self.assertRaises(auxiliary.PyteomicsError, mass.Composition, sequence=mpep, aa_comp=aa_comp)

    def test_composition_term_sseq(self):
        aa_comp = self.aa_comp.copy()
        aa_comp.update(self.mods)
        for pep in self.random_peptides:
            for mlabel, mcomp in self.mods.items():
                split_sequence = parser.parse(pep, split=True)
                self.assertRaises(auxiliary.PyteomicsError, mass.Composition, split_sequence=[
                    (mlabel + '-',) + split_sequence[0]] + split_sequence[1:-1] + [split_sequence[-1] + ('-' + mlabel,)], aa_comp=aa_comp)

    def test_Composition_dict(self):
        # Test Composition from a dict.
        self.assertEqual(mass.Composition(self.d, mass_data=self.mass_data), self.d)

    def test_Composition_formula(self):
        # Test Composition from a formula.
        self.assertEqual(self.d, mass.Composition(formula='ABCDE', mass_data={atom: {0: (1.0, 1.0)} for atom in 'ABCDE'}))

    def test_Composition_seq(self):
        # Test Composition from a sequence.
        self.assertEqual(self.d, mass.Composition(sequence='XYZ', aa_comp=self.aa_comp))

    def test_Composition_pseq(self):
        # Test Composition from a parsed sequence.
        self.assertEqual(
            mass.Composition(parsed_sequence=['X', 'Y', 'Z'], aa_comp=self.aa_comp),
            {atom: 1 for atom in 'ABC'})

    def test_Composition_sseq(self):
        # Test Composition from a split sequence.
        self.assertEqual(
            mass.Composition(split_sequence=[('X',), ('Y',), ('Z',)], aa_comp=self.aa_comp),
            {atom: 1 for atom in 'ABC'})

    def test_Composition_term_formula(self):
        self.assertEqual(mass.Composition(sequence='A2B-XYZ-DE2F3', aa_comp=self.aa_comp),
            {'A': 3, 'B': 2, 'C': 1, 'D': 1, 'E': 2, 'F': 3})

    def test_Composition_nterm_formula(self):
        self.assertEqual(mass.Composition(sequence='AB-XYZ', aa_comp=self.aa_comp),
            {'A': 2, 'B': 2, 'C': 1, 'E': 1})

    def test_Composition_cterm_formula(self):
        self.assertEqual(mass.Composition(sequence='XYZ-AB', aa_comp=self.aa_comp),
            {'A': 2, 'B': 2, 'C': 1, 'D': 1})

    def test_Composition_sum(self):
        # Test sum of Composition objects.
        self.assertEqual(
            mass.Composition(sequence='XXY', aa_comp=self.aa_comp) + mass.Composition(sequence='YZZ', aa_comp=self.aa_comp),
            {atom: 2 for atom in 'ABCDE'})

    def test_Composition_sub(self):
        # Test subtraction of Composition objects
        self.assertEqual({} - mass.Composition(sequence='XYZ', aa_comp=self.aa_comp),
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
        self.assertEqual(mass.Composition('aXbYZ', aa_comp=ac), {'A': 2, 'B': 2, 'C': 1, 'D': 1, 'E': 1})
        self.assertEqual(mass.Composition('AB2C3', mass_data=self.mass_data), {'A': 1, 'B': 2, 'C': 3})

    def test_calculate_mass(self):
        # Calculate mass by a formula.
        self.assertEqual(
            mass.calculate_mass(formula='ABCDE', mass_data=self.mass_data),
            sum(self.mass_data[atom][0][0] for atom in 'ABCDE'))

        # Calculate mass by a sequence.
        self.assertEqual(
            mass.calculate_mass(sequence='XYZ',
                                aa_comp=self.aa_comp,
                                mass_data=self.mass_data),
            sum(self.mass_data[atom][0][0] for atom in 'ABCDE'))

        # Calculate mass by a parsed sequence.
        self.assertEqual(
            mass.calculate_mass(parsed_sequence=['H-', 'X', 'Y', 'Z', '-OH'], aa_comp=self.aa_comp, mass_data=self.mass_data),
            sum(self.mass_data[atom][0][0] for atom in 'ABCDE'))

        # Calculate mass by composition
        self.assertEqual(
            mass.calculate_mass(composition={'A': 1, 'B': 1, 'C': 1, 'D': 1, 'E': 1}, mass_data=self.mass_data),
            sum(self.mass_data[atom][0][0] for atom in 'ABCDE'))

        # Calculate average mass by a formula.
        self.assertEqual(
            mass.calculate_mass(formula='ABCDE', average=True, mass_data=self.mass_data),
            sum(self.mass_data[atom][isotope][0] * self.mass_data[atom][isotope][1]
                for atom in 'ABCDE'for isotope in self.mass_data[atom] if isotope != 0))

        # Calculate m/z of an ion.
        for charge in [1, 2, 3]:
            self.assertEqual(
                mass.calculate_mass(formula='ABCDE', ion_type='M', charge=charge, mass_data=self.mass_data),
                mass.calculate_mass(formula='ABCDE' + 'H+%d' % (charge,), mass_data=self.mass_data))

            self.assertEqual(
                mass.calculate_mass(formula='ABCDE', ion_type='M', charge=charge * 2, charge_carrier='AB+2', mass_data=self.mass_data),
                (mass.calculate_mass(formula='ABCDE', mass_data=self.mass_data) + charge * (
                    self.mass_data['A'][0][0] + self.mass_data['B'][0][0])) / charge / 2)

            self.assertEqual(
                mass.calculate_mass(formula='ABCDE', ion_type='M', charge=charge * 2, charge_carrier={'A': 1, 'B': 1},
                    carrier_charge=2, mass_data=self.mass_data),
                (mass.calculate_mass(formula='ABCDE', mass_data=self.mass_data) + charge * (
                    self.mass_data['A'][0][0] + self.mass_data['B'][0][0])) / charge / 2)

            self.assertEqual(
                mass.calculate_mass(formula='ABCDE', ion_type='M', charge=charge, mass_data=self.mass_data),
                (mass.calculate_mass(formula='ABCDE', mass_data=self.mass_data) + self.mass_data['H+'][0][0] * charge) / charge)

            self.assertAlmostEqual(
                mass.calculate_mass(composition={'A': 1, 'B': 1, 'C': 1, 'D': 1, 'E': 1}, charge=charge, charge_carrier='BC+', mass_data=self.mass_data),
                mass.calculate_mass(composition={'A': 1, 'B': 1 + charge, 'C': 1 + charge, 'D': 1, 'E': 1}, mass_data=self.mass_data) / charge)

            self.assertRaises(auxiliary.PyteomicsError, mass.calculate_mass, **{'formula': 'ABCDEH+%d' % charge,
                   'ion_type': 'M', 'charge': charge, 'mass_data': self.mass_data})

        self.assertRaises(auxiliary.PyteomicsError, mass.calculate_mass, **{'formula': 'ABCDE',
                'ion_type': 'M', 'charge': 3, 'carrier_charge': 2, 'mass_data': self.mass_data})

        # Sanity check.
        for pep in self.random_peptides:
            self.assertEqual(
                mass.calculate_mass(sequence=pep, aa_comp=self.aa_comp, mass_data=self.mass_data, ion_comp=self.ion_comp),
                mass.calculate_mass(parsed_sequence=parser.parse(pep, labels=['X', 'Y', 'Z'], show_unmodified_termini=True),
                    aa_comp=self.aa_comp, mass_data=self.mass_data, ion_comp=self.ion_comp))

    def test_calculate_proforma_mass(self):
        seq_modX = 'PEPTIcamCIDE'
        aa_comp = mass.std_aa_comp.copy()
        aa_comp['cam'] = mass.Composition(formula='H3C2NO')
        seq_proforma = 'PEPTIC[+57.021464]IDE'
        for charge in [None, 0, 1, 2]:
            self.assertAlmostEqual(
                mass.calculate_mass(sequence=seq_modX, charge=charge, aa_comp=aa_comp),
                mass.calculate_mass(proforma=seq_proforma, charge=charge),
                places=6)

    def test_most_probable_isotopic_composition(self):
        self.assertEqual(
            mass.most_probable_isotopic_composition(formula='F', mass_data=self.mass_data),
            (mass.Composition({'F[6]': 1, 'F[7]': 0}, mass_data=self.mass_data), 0.7))

        self.assertEqual(
            mass.most_probable_isotopic_composition(formula='F10', mass_data=self.mass_data),
            (mass.Composition({'F[6]': 7, 'F[7]': 3}, mass_data=self.mass_data), (0.3)**3 * (0.7)**7 * 120))

        self.assertEqual(
            mass.most_probable_isotopic_composition(formula='A20F10', elements_with_isotopes=['F'], mass_data=self.mass_data),
            (mass.Composition({'A': 20, 'F[6]': 7, 'F[7]': 3}, mass_data=self.mass_data), (0.3)**3 * (0.7)**7 * 120))

    def test_isotopic_composition_abundance(self):
        for peplen in range(1, 10):
            self.assertAlmostEqual(
                mass.isotopic_composition_abundance(formula='F[6]' * peplen, mass_data=self.mass_data),
                self.mass_data['F'][6][1] ** peplen)

            self.assertAlmostEqual(
                mass.isotopic_composition_abundance(formula='AF[6]' * peplen, mass_data=self.mass_data),
                self.mass_data['F'][6][1] ** peplen)

            self.assertAlmostEqual(
                mass.isotopic_composition_abundance(formula='A[1]F[6]' * peplen, mass_data=self.mass_data),
                (self.mass_data['A'][1][1] * self.mass_data['F'][6][1]) ** peplen)

    def test_Unimod_mass(self):
        db = mass.Unimod(gzip.open('unimod.xml.gz'))
        for x in db.mods:
            self.assertGreater(0.00001,
                abs(x['mono_mass'] - mass.calculate_mass(x['composition'], mass_data=db.mass_data)))

    def test_Unimod_methods(self):
        db = mass.Unimod(gzip.open('unimod.xml.gz'))
        rec_id = 1
        rec_name = 'Acetylation'
        rec_title = 'Acetyl'
        record = db.by_id(rec_id)
        self.assertEqual(record['title'], rec_title)
        self.assertEqual(record['full_name'], rec_name)
        self.assertEqual(record, db[rec_id])
        self.assertEqual(record, db.by_title(rec_title))
        self.assertEqual(record, db.by_name(rec_name))

    def test_nist_mass(self):
        self.assertTrue(all(abs(g[0][1] - 1) < 1e-6 for g in mass.nist_mass.values()))
        for g in mass.nist_mass.values():
            s = sum(p[1] for num, p in g.items() if num)
            self.assertTrue(abs(s - 1) < 1e-6 or abs(s) < 1e-6)

    def test_composition_objects_are_pickleable(self):
        dict_ = mass.Composition(self.d, mass_data=self.mass_data)
        formula = mass.Composition(formula='ABCDE',
                         mass_data={atom: {0: (1.0, 1.0)} for atom in 'ABCDE'})
        sequence = mass.Composition(sequence='XYZ', aa_comp=self.aa_comp)
        parsed_sequence = mass.Composition(parsed_sequence=['X', 'Y', 'Z'],
                             aa_comp=self.aa_comp)
        split_sequence = mass.Composition(split_sequence=[('X',), ('Y',), ('Z',)],
                             aa_comp=self.aa_comp)

        self.assertEqual(dict_, pickle.loads(pickle.dumps(dict_)))
        self.assertEqual(formula, pickle.loads(pickle.dumps(formula)))
        self.assertEqual(sequence, pickle.loads(pickle.dumps(sequence)))
        self.assertEqual(parsed_sequence, pickle.loads(pickle.dumps(parsed_sequence)))
        self.assertEqual(split_sequence, pickle.loads(pickle.dumps(split_sequence)))

    def test_aa_mass(self):
        h2o = mass.calculate_mass(formula='H2O')
        for aa, m in mass.std_aa_mass.items():
            self.assertEqual(m + h2o, mass.fast_mass(aa))

    def test_isotopologues(self):
        peptide = 'XYF'
        states = [{'F[6]': 1, 'A': 1, 'B': 1, 'D': 1, 'E': 1}, {'F[7]': 1, 'A': 1, 'B': 1, 'D': 1, 'E': 1}]
        abundances = [0.7, 0.3]
        kw_common = dict(elements_with_isotopes='F', aa_comp=self.aa_comp, mass_data=self.mass_data)
        kwlist = [
            {},
            {'sequence': 'XYF'},
            {'parsed_sequence': parser.parse('XYF', show_unmodified_termini=True)},
            {'split_sequence': parser.parse('XYF', show_unmodified_termini=True, split=True)},
            {'formula': 'ABDEF'},
            {'composition': mass.Composition(sequence='XYF', aa_comp=self.aa_comp)}]
        arglist = [(peptide,), (), (), (), (), ()]
        for args, kw in zip(arglist, kwlist):
            kwargs = kw_common.copy()
            kwargs.update(kw)
            isotopologues = mass.isotopologues(*args, **kwargs)
            for state in isotopologues:
                i = states.index(state)
                self.assertNotEqual(i, -1)
                self.assertAlmostEqual(abundances[i], mass.isotopic_composition_abundance(state,
                    aa_comp=self.aa_comp, mass_data=self.mass_data))

    def test_isotopologues_with_abundances(self):
        peptide = 'XYF'
        states = [{'F[6]': 1, 'A': 1, 'B': 1, 'D': 1, 'E': 1}, {'F[7]': 1, 'A': 1, 'B': 1, 'D': 1, 'E': 1}]
        abundances = [0.7, 0.3]
        for state, abundance in mass.isotopologues(peptide, elements_with_isotopes='F',
                aa_comp=self.aa_comp, mass_data=self.mass_data, report_abundance=True):
            i = states.index(state)
            self.assertNotEqual(i, -1)
            self.assertAlmostEqual(abundances[i], abundance)

    def test_std_aa_mass(self):
        for key, value in mass.std_aa_mass.items():
            self.assertAlmostEqual(value, mass.calculate_mass(parsed_sequence=[key]), places=4)

    def test_fragment_series_modx_vs_proforma(self):
        """Test that fragment_series produces identical results for modX and ProForma formats."""
        # Define sequences with phosphorylation on threonine
        # modX format: phosphorylated threonine
        seq_modx = 'PEPpTIDE'
        # ProForma format: threonine with +79.966331 (phosphorylation)
        seq_proforma = 'PEPT[+79.966331]IDE'

        # Define amino acid compositions including phosphorylated threonine
        aa_comp_modx = mass.std_aa_comp.copy()
        aa_comp_modx['pT'] = mass.Composition(formula='H8C4O5NP')  # T + phosphate group

        # Test parameters
        ion_types = ('b', 'y')
        maxcharge = 2

        # Generate fragments for modX sequence
        fragments_modx = mass.fragment_series(
            seq_modx,
            ion_types=ion_types,
            maxcharge=maxcharge,
            aa_mass={k: mass.calculate_mass(composition=v) for k, v in aa_comp_modx.items()}
        )

        # Generate fragments for ProForma sequence
        fragments_proforma = mass.fragment_series(
            seq_proforma,
            ion_types=ion_types,
            maxcharge=maxcharge
        )

        # Verify that both formats produce the same ion types
        self.assertEqual(set(fragments_modx.keys()), set(fragments_proforma.keys()),
                        "Ion types should be identical for modX and ProForma")

        # Verify that m/z values are identical (within tolerance)
        for ion_type in ion_types:
            modx_fragments = fragments_modx[ion_type]
            proforma_fragments = fragments_proforma[ion_type]

            self.assertEqual(len(modx_fragments), len(proforma_fragments),
                           f"Number of {ion_type} ions should be identical")

            # Check that fragment names match
            self.assertEqual(set(modx_fragments), set(proforma_fragments),
                           f"{ion_type} ion names should be identical")

            # Check m/z values are close (allowing for small numerical differences)
            for fragment_name in modx_fragments:
                self.assertAlmostEqual(modx_fragments[fragment_name], proforma_fragments[fragment_name], places=5,
                                     msg=f"{ion_type} ion {fragment_name} m/z values should be identical")

        # Additional checks to ensure we got reasonable results
        self.assertGreater(len(fragments_modx['b']), 0, "Should generate b ions")
        self.assertGreater(len(fragments_modx['y']), 0, "Should generate y ions")

        # Check that we have the expected number of fragments for a 7-residue peptide
        # For each charge state, we expect 6 b ions and 6 y ions
        expected_fragments_per_charge = 6
        expected_total_fragments = expected_fragments_per_charge * maxcharge

        self.assertEqual(len(fragments_modx['b']), expected_total_fragments,
                        f"Expected {expected_total_fragments} b ions")
        self.assertEqual(len(fragments_modx['y']), expected_total_fragments,
                        f"Expected {expected_total_fragments} y ions")

    def test_fragment_series_basic_functionality(self):
        """Test basic functionality of fragment_series with a simple sequence."""
        seq = 'PEPTIDE'

        # Test with default parameters
        fragments = mass.fragment_series(seq)

        # Should have b and y ions by default
        self.assertIn('b', fragments)
        self.assertIn('y', fragments)

        # For PEPTIDE (7 residues), with maxcharge=1, expect 6 b ions and 6 y ions
        self.assertEqual(len(fragments['b']), 6)
        self.assertEqual(len(fragments['y']), 6)

        # Check that names are correct
        expected_b_names = {'b1+', 'b2+', 'b3+', 'b4+', 'b5+', 'b6+'}
        expected_y_names = {'y1+', 'y2+', 'y3+', 'y4+', 'y5+', 'y6+'}
        self.assertEqual(set(fragments['b'].keys()), expected_b_names)
        self.assertEqual(set(fragments['y'].keys()), expected_y_names)

        # Check that m/z values are reasonable (all positive)
        for ion_type in ['b', 'y']:
            for fragment_name, mz_val in fragments[ion_type].items():
                self.assertGreater(mz_val, 0, f"m/z value for {fragment_name} should be positive")

        # Test with different ion types
        fragments_abc = mass.fragment_series(seq, ion_types=('a', 'b', 'c'))
        self.assertIn('a', fragments_abc)
        self.assertIn('b', fragments_abc)
        self.assertIn('c', fragments_abc)

        # Test with higher charge states
        fragments_z2 = mass.fragment_series(seq, maxcharge=2)
        # With maxcharge=2, should have twice as many fragments
        self.assertEqual(len(fragments_z2['b']), 12)  # 6 fragments × 2 charge states
        self.assertEqual(len(fragments_z2['y']), 12)


if __name__ == '__main__':
    unittest.main()
