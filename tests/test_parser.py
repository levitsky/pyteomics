import unittest
from pyteomics.parser import *
from string import ascii_uppercase as uppercase
import random
class ParserTest(unittest.TestCase):
    def setUp(self):
        self.simple_sequences = [''.join([random.choice(uppercase) for i in range(
            int(random.uniform(1, 20)))]) for j in range(10)]
        self.labels = ['A', 'B', 'C', 'N', 'X']
        self.extlabels = self.labels[:]
        self.potential = {'pot': ['X', 'A', 'B'], 'otherpot': ['A', 'C']}
        self.constant = {'const': ['B']}
        self.extlabels.extend(['potX', 'potA', 'potB', 'otherpotA', 'otherpotC'])
        self.extlabels.extend(['constB'])

    def test_parse_sequence_simple(self):
        for seq in self.simple_sequences:
            self.assertEqual(seq, ''.join(parse_sequence(seq, labels=uppercase)))

    def test_amino_acid_composition_simple(self):
        for seq in self.simple_sequences:
            comp = amino_acid_composition(seq, labels=uppercase)
            for aa in set(seq):
                self.assertEqual(seq.count(aa), comp[aa])

    def test_isoforms_simple(self):
        self.assertEqual(
                isoforms('PEPTIDE', variable_mods={'xx': ['A', 'B', 'P', 'E']}),
                set(['PEPTIDE', 'PEPTIDxxE', 'PExxPTIDE', 'PExxPTIDxxE', 'PxxEPTIDE',
                     'PxxEPTIDxxE', 'PxxExxPTIDE', 'PxxExxPTIDxxE', 'xxPEPTIDE', 'xxPEPTIDxxE',
                     'xxPExxPTIDE', 'xxPExxPTIDxxE', 'xxPxxEPTIDE', 'xxPxxEPTIDxxE',
                     'xxPxxExxPTIDE', 'xxPxxExxPTIDxxE']))

    def test_isoforms_len(self):
        for j in range(10):
            L = random.randint(5, 15)
            peptide = ''.join([random.choice(self.labels) for _ in range(L)])
            modseqs = isoforms(peptide, variable_mods=self.potential,
                    fixed_mods=self.constant, labels=self.labels)
            pp = parse_sequence(peptide, labels=self.labels)
            for p in modseqs:
                self.assertEqual(len(pp),
                        len(parse_sequence(p, labels=self.extlabels)))
            self.assertEqual(len(modseqs), (3**pp.count('A'))*
                    (2**(pp.count('X')+pp.count('C'))))

if __name__ == '__main__':
    unittest.main()
