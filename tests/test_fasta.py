import tempfile
import unittest
from pyteomics.fasta import *
import random
import string

class FastaTest(unittest.TestCase):
    def setUp(self):
        self.fasta_file = tempfile.TemporaryFile(mode='r+')
        self.fasta_file.write('''
            ;test sequence
            ;test sequence 2
            TEST
            
            >test sequence 3
            TE
            ST*
            >test sequence 4
            TEST
            ''')
        self.fasta_file.seek(0)
        self.fasta_entries_short = list(read(self.fasta_file, ignore_comments=True))
        self.fasta_file.seek(0)
        self.fasta_entries_long = list(read(self.fasta_file))

    def test_simple_read_long_comments(self):
        self.assertEqual(self.fasta_entries_long,
                         [('test sequence test sequence 2', 'TEST'),
                          ('test sequence 3', 'TEST'),
                          ('test sequence 4', 'TEST')])

    def test_simple_read_short_comments(self):
        self.assertEqual(self.fasta_entries_short,
                         [('test sequence', 'TEST'),
                          ('test sequence 3', 'TEST'),
                          ('test sequence 4', 'TEST')])

    def test_decoy_sequence_reverse(self):
        sequence = ''.join(random.choice(string.ascii_uppercase)
                             for i in range(random.randint(1,100)))
        self.assertEqual(decoy_sequence(sequence, 'reverse'),
                sequence[::-1])

    def test_decoy_sequence_shuffle(self):
        sequences = (''.join(random.choice(string.ascii_uppercase)
                             for i in range(random.randint(1,100)))
                                for j in range(10))
        test = True
        for s in sequences:
            ss = decoy_sequence(s, 'shuffle')
            self.assertEqual(sorted(list(s)), sorted(list(ss)))
            if not all(a == b for a, b in zip(s, ss)):
                test = False
        self.assertFalse(test)

    def test_read_and_write_fasta_short(self):
        self.fasta_file.seek(0)
        new_fasta_file = tempfile.TemporaryFile(mode='r+')
        write(read(self.fasta_file, ignore_comments=True),
                new_fasta_file)
        new_fasta_file.seek(0)
        new_entries = list(read(new_fasta_file, ignore_comments=True))
        self.fasta_file.seek(0)
        self.assertEqual(new_entries, self.fasta_entries_short)
        new_fasta_file.close()

    def test_read_and_write_long(self):
        self.fasta_file.seek(0)
        new_fasta_file = tempfile.TemporaryFile(mode='r+')
        write(read(self.fasta_file), new_fasta_file)
        new_fasta_file.seek(0)
        new_entries = list(read(new_fasta_file))
        self.fasta_file.seek(0)
        self.assertEqual(new_entries, self.fasta_entries_long)
        new_fasta_file.close()
        
    def test_write_decoy_db(self):
        self.fasta_file.seek(0)
        decdb = tempfile.TemporaryFile(mode='r+')
        write_decoy_db(self.fasta_file, decdb, decoy_only=False, prefix='PREFIX_')
        decdb.seek(0)
        all_entries = list(read(decdb, False))
        decdb.close()
        self.assertEqual(all_entries, self.fasta_entries_long + 
                [('PREFIX_' + a, b[::-1]) for a, b in self.fasta_entries_long])
        self.fasta_file.seek(0)

if __name__ == '__main__':
    unittest.main()
                         
