import unittest
from pyteomics.parser import *
from string import uppercase
import random
class ParserTest(unittest.TestCase):
    def setUp(self):
        self.simple_sequences = [''.join([random.choice(uppercase) for i in range(
            int(random.uniform(1, 20)))]) for j in range(10)]
        self.labels = ['A', 'B', 'C', 'N', 'X']
        self.extlabels = self.labels[:]
        self.potential = {'pot': ['X', 'A', 'B'], 'otherpot': ['A', 'C']}
        self.constant = {'const': ['B']}
        add_modifications(self.extlabels, self.potential)
        add_modifications(self.extlabels, self.constant)

    def test_parse_sequence_simple(self):
        for seq in self.simple_sequences:
            self.assertEqual(seq, ''.join(parse_sequence(seq, labels=uppercase)))

    def test_amino_acid_composition_simple(self):
        for seq in self.simple_sequences:
            comp = amino_acid_composition(seq, labels=uppercase)
            for aa in set(seq):
                self.assertEqual(seq.count(aa), comp[aa])

    def test_modify_peptide_simple(self):
        self.assertEqual(
                modify_peptide('PEPTIDE', potential={'xx': ['A', 'B', 'P', 'E']}),
                set(['PEPTIDE', 'PEPTIDxxE', 'PExxPTIDE', 'PExxPTIDxxE', 'PxxEPTIDE',
                     'PxxEPTIDxxE', 'PxxExxPTIDE', 'PxxExxPTIDxxE', 'xxPEPTIDE', 'xxPEPTIDxxE',
                     'xxPExxPTIDE', 'xxPExxPTIDxxE', 'xxPxxEPTIDE', 'xxPxxEPTIDxxE',
                     'xxPxxExxPTIDE', 'xxPxxExxPTIDxxE']))

    def test_add_modifications(self):
        self.assertEqual(len(self.extlabels)-len(self.labels),
            sum([len(x) for x in self.potential.values()]) +
            sum([len(x) for x in self.constant.values()]))

    def test_modify_peptide_len(self):
        for j in xrange(10):
            L = random.randint(5, 15)
            peptide = ''.join([random.choice(self.labels) for _ in xrange(L)])
            modseqs = modify_peptide(peptide, potential=self.potential,
                    constant=self.constant, labels=self.labels)
            pp = parse_sequence(peptide, labels=self.labels)
            for p in modseqs:
                self.assertEqual(len(pp),
                        len(parse_sequence(p, labels=self.extlabels)))
            self.assertEqual(len(modseqs), (3**pp.count('A'))*
                    (2**(pp.count('X')+pp.count('C'))))

if __name__ == '__main__':
    unittest.main()
