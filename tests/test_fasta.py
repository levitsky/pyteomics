import tempfile
import unittest
from pyteomics.fasta import read_fasta

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

    def test_simple_read_long_comments(self):
        fasta_entries = [i for i in read_fasta(self.fasta_file, False)]
        self.assertEqual(fasta_entries,
                         [('test sequence test sequence 2', 'TEST'),
                          ('test sequence 3', 'TEST'),
                          ('test sequence 4', 'TEST')])
        self.fasta_file.seek(0)
    def test_simple_read_short_comments(self):
        fasta_entries = [i for i in read_fasta(self.fasta_file)]
        self.assertEqual(fasta_entries,
                         [('test sequence', 'TEST'),
                          ('test sequence 3', 'TEST'),
                          ('test sequence 4', 'TEST')])
        self.fasta_file.seek(0)
        
if __name__ == '__main__':
    unittest.main()
                         
