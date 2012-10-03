import unittest
from pyteomics.parser import *
from string import ascii_uppercase as uppercase
import random
class ParserTest(unittest.TestCase):
    def setUp(self):
        self.simple_sequences = [''.join(random.choice(uppercase) for i in range(
            int(random.uniform(1, 20)))) for j in range(10)]
        self.labels = ['A', 'B', 'C', 'N', 'X']
        self.extlabels = self.labels[:]
        self.potential = {'pot': ['X', 'A', 'B'], 'otherpot': ['A', 'C'], 'N-': ['N'], '-C': ['C']}
        self.constant = {'const': ['B']}
        self.extlabels.extend(('pot', 'otherpot', 'const', '-C', 'N-'))

    def test_parse_simple(self):
        for seq in self.simple_sequences:
            self.assertEqual(seq, ''.join(parse(seq, labels=uppercase)))

    def test_parse(self):
        self.assertEqual(
                [('P',), ('E',), ('P',), ('T',), ('I',), ('D',), ('E',)],
                parse('PEPTIDE', split=True))
        self.assertEqual(['P', 'E', 'P', 'T', 'I', 'D', 'E'],
                parse('H-PEPTIDE'))
        self.assertEqual(['H-', 'P', 'E', 'P', 'T', 'I', 'D', 'E', '-OH'],
                parse('PEPTIDE', show_unmodified_termini=True))
        self.assertEqual(['T', 'E', 'pS', 'T', 'oxM'],
                parse('TEpSToxM', labels=std_labels + ['pS', 'oxM']))
        self.assertEqual(
                [('H-', 'z', 'P'), ('E',), ('P',), ('z', 'T'), ('I',), ('D',),
                    ('z', 'E', '-OH')],
                parse(
                    'zPEPzTIDzE', True, True, labels=std_labels+['z']))

    def test_tostring(self):
        for seq in self.simple_sequences:
            self.assertEqual(seq, tostring(parse(seq, labels=uppercase)))
            self.assertEqual(seq, tostring(parse(
                seq, True, True, labels=uppercase), False))

    def test_amino_acid_composition_simple(self):
        for seq in self.simple_sequences:
            comp = amino_acid_composition(seq, labels=uppercase)
            for aa in set(seq):
                self.assertEqual(seq.count(aa), comp[aa])

    def test_amino_acid_composition(self):
        for seq in self.simple_sequences:
            comp = amino_acid_composition(seq, term_aa=True, labels=uppercase)
            comp_default = amino_acid_composition(seq, labels=uppercase)
            self.assertEqual(1, comp['nterm'+seq[0]])
            if len(seq) > 1:
                self.assertEqual(1, comp['cterm'+seq[-1]])
            self.assertEqual(sum(comp_default.values()),
                             sum(comp.values()))

    def test_cleave(self):
        for seq in self.simple_sequences:
            for elem in cleave(
                    seq, expasy_rules['trypsin'], int(random.uniform(1, 10))):
                self.assertIn(elem, seq)
            self.assertTrue(any(elem == seq
                for elem in cleave(seq, expasy_rules['trypsin'], len(seq))))

    def test_isoforms_simple(self):
        self.assertEqual(
                set(isoforms('PEPTIDE', variable_mods={'xx': ['A', 'B', 'P', 'E']})),
                set(['PEPTIDE', 'PEPTIDxxE', 'PExxPTIDE', 'PExxPTIDxxE', 'PxxEPTIDE',
                     'PxxEPTIDxxE', 'PxxExxPTIDE', 'PxxExxPTIDxxE', 'xxPEPTIDE', 'xxPEPTIDxxE',
                     'xxPExxPTIDE', 'xxPExxPTIDxxE', 'xxPxxEPTIDE', 'xxPxxEPTIDxxE',
                     'xxPxxExxPTIDE', 'xxPxxExxPTIDxxE']))

    def test_isoforms_len(self):
        for j in range(50):
            L = random.randint(1, 10)
            peptide = ''.join([random.choice(self.labels) for _ in range(L)])
#           print peptide
            modseqs = isoforms(peptide, variable_mods=self.potential,
                    fixed_mods=self.constant, labels=self.labels)
            forms = sum(1 for x in modseqs)
            pp = parse(peptide, labels=self.extlabels)
            N = 0
            if pp[0] =='N': N += 1
            if pp[-1] == 'C': N += 1
            for p in modseqs:
                self.assertEqual(len(pp),
                        peptide_length(p, labels=self.extlabels))
            self.assertEqual(forms, (3**pp.count('A')) *
                    (2**(pp.count('X')+pp.count('C'))) * 2**N)

if __name__ == '__main__':
    unittest.main()
