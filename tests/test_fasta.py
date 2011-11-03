import tempfile
import unittest
from pyteomics.fasta import *

class FastaTest(unittest.TestCase):
    def setUp(self):
        self.fasta_file = tempfile.TemporaryFile()
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
        self.fasta_entries = [i for i in read_fasta(self.fasta_file)]
        self.fasta_file.seek(0)
        self.fasta_entries_long = [i for i in read_fasta(self.fasta_file, False)]

    def test_simple_read_long_comments(self):
        self.assertEqual(self.fasta_entries_long,
                         [('test sequence test sequence 2', 'TEST'),
                          ('test sequence 3', 'TEST'),
                          ('test sequence 4', 'TEST')])
    def test_simple_read_short_comments(self):
        self.assertEqual(self.fasta_entries,
                         [('test sequence', 'TEST'),
                          ('test sequence 3', 'TEST'),
                          ('test sequence 4', 'TEST')])
    def test_decoy_sequence_reverse(self):
        sequence = 'ABCDEF123'
        self.assertEqual(decoy_sequence(sequence, 'reverse'),
                sequence[::-1])
    def test_read_and_write_fasta(self):
        self.fasta_file.seek(0)
        new_fasta_file = tempfile.TemporaryFile()
        write_fasta(read_fasta(self.fasta_file), new_fasta_file, False)
        new_fasta_file.seek(0)
        self.assertEqual(self.fasta_file.read(),
                new_fasta_file.read())
        self.fasta_file.seek(0)




if __name__ == '__main__':
    unittest.main()
                         
