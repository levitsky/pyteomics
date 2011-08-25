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

    def test_simple_read(self):
        fasta_entries = [i for i in read_fasta(self.fasta_file)]
        self.assertEqual(fasta_entries,
                         [('a test sequencea test sequence 2', 'TEST'),
                          ('a test sequence 3', 'TEST'),
                          ('a test sequence 4', 'TEST')])
        
if __name__ == '__main__':
    unittest.main()
                         
