from os import path
import pyteomics
pyteomics.__path__ = [path.abspath(path.join(path.dirname(__file__), path.pardir, 'pyteomics'))]
import unittest
from pyteomics import parser
from string import ascii_uppercase as uppercase
import random


class ParserTest(unittest.TestCase):
    def setUp(self):
        self.simple_sequences = [''.join(random.choice(uppercase) for i in range(
            int(random.uniform(1, 20)))) for j in range(10)]
        self.labels = ['A', 'B', 'C', 'N', 'X']
        self.extlabels = self.labels[:]
        self.potential = {'pot': ['X', 'A', 'B'], 'otherpot': ['A', 'C'],
                'N-': ['N'], '-C': ['C']}
        self.constant = {'const': ['B']}
        self.extlabels.extend(('pot', 'otherpot', 'const', '-C', 'N-'))

    def test_parse_simple(self):
        for seq in self.simple_sequences:
            self.assertEqual(seq, ''.join(parser.parse(seq, labels=uppercase)))

    def test_parse(self):
        self.assertEqual(
            [('P',), ('E',), ('P',), ('T',), ('I',), ('D',), ('E',)],
            parser.parse('PEPTIDE', split=True))
        self.assertEqual(['P', 'E', 'P', 'T', 'I', 'D', 'E'],
            parser.parse('H-PEPTIDE'))
        for seq in ['PEPTIDE', 'H-PEPTIDE', 'PEPTIDE-OH', 'H-PEPTIDE-OH']:
            self.assertEqual(['H-', 'P', 'E', 'P', 'T', 'I', 'D', 'E', '-OH'],
                parser.parse(seq, show_unmodified_termini=True))
        self.assertEqual(['T', 'E', 'pS', 'T', 'oxM'],
                parser.parse('TEpSToxM', labels=parser.std_labels + ['pS', 'oxM']))
        self.assertEqual(
            [('H-', 'z', 'P'), ('E',), ('P',), ('z', 'T'), ('I',), ('D',), ('z', 'E', '-OH')],
            parser.parse('zPEPzTIDzE', True, True, labels=parser.std_labels + ['z']))

    def test_tostring(self):
        for seq in self.simple_sequences:
            self.assertEqual(seq, parser.tostring(parser.parse(seq, labels=uppercase)))
            self.assertEqual(seq, parser.tostring(parser.parse(
                seq, True, True, labels=uppercase), False))

    def test_amino_acid_composition_simple(self):
        for seq in self.simple_sequences:
            comp = parser.amino_acid_composition(seq, labels=uppercase)
            for aa in set(seq):
                self.assertEqual(seq.count(aa), comp[aa])

    def test_amino_acid_composition(self):
        for seq in self.simple_sequences:
            comp = parser.amino_acid_composition(seq, term_aa=True, labels=uppercase)
            comp_default = parser.amino_acid_composition(seq, labels=uppercase)
            self.assertEqual(1, comp['nterm' + seq[0]])
            if len(seq) > 1:
                self.assertEqual(1, comp['cterm' + seq[-1]])
            self.assertEqual(sum(comp_default.values()), sum(comp.values()))

    def test_cleave(self):
        self.assertEqual(parser.xcleave('PEPTIDEKS', parser.expasy_rules['trypsin']), [(0, 'PEPTIDEK'), (8, 'S')])
        self.assertEqual(parser.xcleave('PEPTIDEKS', 'trypsin'), [(0, 'PEPTIDEK'), (8, 'S')])
        self.assertEqual(parser.xcleave('PEPTIDEKS', 'Trypsin'), [(0, 'PEPTIDEK'), (8, 'S')])
        for seq in self.simple_sequences:
            for elem in parser.cleave(
                    seq, 'trypsin', int(random.uniform(1, 10))):
                self.assertIn(elem, seq)
            self.assertTrue(any(elem == seq
                for elem in parser.cleave(seq, parser.expasy_rules['trypsin'], len(seq))))

    def test_cleave_semi(self):
        self.assertEqual(parser.xcleave('PEPTIDEKS', 'trypsin', semi=True),
            [(0, 'PEPTIDEK'), (0, 'P'), (0, 'PE'), (0, 'PEP'), (0, 'PEPT'), (0, 'PEPTI'), (0, 'PEPTID'), (0, 'PEPTIDE'),
             (1, 'EPTIDEK'), (2, 'PTIDEK'), (3, 'TIDEK'), (4, 'IDEK'), (5, 'DEK'), (6, 'EK'), (7, 'K'), (8, 'S')])
        self.assertEqual(parser.cleave('PEPTIDEKS', parser.expasy_rules['trypsin'], semi=True),
            {'PEPTIDEK', 'P', 'PE', 'PEP', 'PEPT', 'PEPTI', 'PEPTID', 'PEPTIDE', 'EPTIDEK', 'PTIDEK', 'TIDEK', 'IDEK', 'DEK', 'EK', 'K', 'S'})

    def test_cleave_min_length(self):
        for seq in self.simple_sequences:
            ml = random.uniform(1, 5)
            for elem in parser.cleave(
                    seq, parser.expasy_rules['trypsin'], int(random.uniform(1, 10)), ml):
                self.assertTrue(len(elem) >= ml)

    def test_num_sites(self):
        self.assertEqual(parser.num_sites('RKCDE', 'K'), 1)
        self.assertEqual(parser.num_sites('RKCDE', 'E'), 0)
        self.assertEqual(parser.num_sites('RKCDE', 'R'), 1)
        self.assertEqual(parser.num_sites('RKCDE', 'Z'), 0)

    def test_isoforms_simple(self):
        self.assertEqual(
            list(parser.isoforms('PEPTIDE', variable_mods={'xx': ['A', 'B', 'P', 'E']})),
            ['PEPTIDE', 'PEPTIDxxE', 'PExxPTIDE', 'PExxPTIDxxE', 'PxxEPTIDE', 'PxxEPTIDxxE', 'PxxExxPTIDE',
            'PxxExxPTIDxxE', 'xxPEPTIDE', 'xxPEPTIDxxE', 'xxPExxPTIDE', 'xxPExxPTIDxxE', 'xxPxxEPTIDE',
            'xxPxxEPTIDxxE', 'xxPxxExxPTIDE', 'xxPxxExxPTIDxxE'])

    def test_isoforms_fixed_simple(self):
        self.assertEqual(
            list(parser.isoforms('PEPTIDE', fixed_mods={'n-': True, '-c': True, 'x': ['P', 'T']})),
            ['n-xPExPxTIDE-c'])

    def test_isoforms_simple_2(self):
        self.assertEqual(list(parser.isoforms('PEPTIDE', variable_mods={'x': 'T', 'y': 'T'})),
            ['PEPTIDE', 'PEPxTIDE', 'PEPyTIDE'])

    def test_isoforms_universal(self):
        self.assertEqual(set(parser.isoforms('PEPTIDE', variable_mods={'xx-': True})), {'PEPTIDE', 'xx-PEPTIDE'})
        self.assertEqual(set(parser.isoforms('PEPTIDE', variable_mods={'-xx': True})), {'PEPTIDE', 'PEPTIDE-xx'})
        for seq in self.simple_sequences:
            self.assertEqual(sum(1 for _ in parser.isoforms(seq, variable_mods={'x': True})), 2**len(seq))

    def test_isoforms_terminal(self):
        self.assertEqual(set(parser.isoforms('PEPTIDE', variable_mods={'xx': ['ntermP'], 'yy-': 'P'})),
            {'PEPTIDE', 'xxPEPTIDE', 'yy-PEPTIDE', 'yy-xxPEPTIDE'})

    def test_isoforms_len(self):
        for j in range(50):
            L = random.randint(1, 10)
            peptide = ''.join(random.choice(self.labels) for _ in range(L))
            modseqs = list(parser.isoforms(peptide, variable_mods=self.potential,
                    fixed_mods=self.constant, labels=self.labels))
            pp = parser.parse(peptide, labels=self.extlabels)
            N = (pp[0] == 'N') + (pp[-1] == 'C')
            for p in modseqs:
                self.assertEqual(len(pp), parser.length(p, labels=self.extlabels))
            self.assertEqual(len(modseqs), (3 ** pp.count('A')) * (2 ** (pp.count('X') + pp.count('C') + N)))

    def test_isoforms_maxmods(self):
        for j in range(50):
            L = random.randint(1, 10)
            M = random.randint(1, 10)
            peptide = ''.join([random.choice(self.labels) for _ in range(L)])
            modseqs = parser.isoforms(peptide, variable_mods=self.potential,
                    labels=self.labels, max_mods=M, format='split')
            pp = parser.parse(peptide, labels=self.extlabels, split=True)
            for ms in modseqs:
                self.assertEqual(len(pp), len(ms))
                self.assertLessEqual(sum(i != j for i, j in zip(pp, ms)), M)

    def test_fast_valid(self):
        for j in range(50):
            L = random.randint(1, 10)
            peptide = ''.join([random.choice(self.labels) for _ in range(L)])
            self.assertTrue(parser.fast_valid(peptide, labels=self.labels))
            self.assertTrue(parser.valid(peptide, labels=self.labels))
            self.assertTrue(parser.valid(peptide))
            for aa in set(peptide):
                bad = peptide.replace(aa, 'Z')
                self.assertFalse(parser.fast_valid(bad, labels=self.labels))
                self.assertFalse(parser.valid(bad, labels=self.labels))

    def test_valid(self):
        for j in range(50):
            L = random.randint(1, 10)
            peptide = ''.join([random.choice(self.labels) for _ in range(L)])
            modseqs = parser.isoforms(peptide, variable_mods=self.potential,
                    fixed_mods=self.constant, labels=self.labels)
            self.assertFalse(parser.valid('H-' + peptide, labels=self.labels))
            for s in modseqs:
                self.assertTrue(parser.valid(s, labels=self.extlabels))
                for aa in set(peptide):
                    bad = s.replace(aa, 'Z')
                    self.assertFalse(parser.fast_valid(bad, labels=self.labels))
                    self.assertFalse(parser.valid(bad, labels=self.labels))


if __name__ == '__main__':
    import doctest
    doctest.testmod(parser)
    unittest.main()
