import unittest
from pyteomics.parser import *
from string import uppercase
import random
class ParserTest(unittest.TestCase):
    def setUp(self):
        self.simple_sequences = [''.join([random.choice(uppercase) for i in range(
            int(random.uniform(1, 20)))]) for j in range(10)]

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

#    def test_modify_peptide_len(self):
#        L = random.randint([5, 15])
#        labels = ['A', 'B', 'C', 'N']
#        peptide = ''.join([random.choice(labels) for _ in xrange(L)])
#        potential = {'pot': ['X', 'A', 'B'], 'otherpot': ['C', 'D', 'E']}
#        constant = {'const': ['B']}
#        modseqs = modify_peptide(peptide, potential=potential,
#                constant=constant, labels=labels)
#        for p in modseqs:
#            self.assertEqual(len(

if __name__ == '__main__':
    unittest.main()
