from os import path
import pyteomics
pyteomics.__path__ = [path.abspath(path.join(path.dirname(__file__), path.pardir, 'pyteomics'))]
import tempfile
import unittest
from pyteomics.fasta import *
import random
import string

class FastaTest(unittest.TestCase):
    def setUp(self):
        self.fasta_file = 'test.fasta'
        self.fasta_entries_short = list(read(self.fasta_file,
            ignore_comments=True))
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
                             for i in range(random.randint(1, 50)))
        self.assertEqual(decoy_sequence(sequence, 'reverse'),
                sequence[::-1])

    def test_decoy_sequence_shuffle(self):
        sequences = [''.join(random.choice(string.ascii_uppercase)
                             for i in range(random.randint(1, 50)))
                                for j in range(10)]
        test = True
        for s in sequences:
            ss = decoy_sequence(s, 'shuffle')
            self.assertEqual(sorted(list(s)), sorted(list(ss)))
            if not all(a == b for a, b in zip(s, ss)):
                test = False
        self.assertFalse(test)

    def test_decoy_sequence_fused(self):
        sequences = [''.join(random.choice(string.ascii_uppercase)
                             for i in range(random.randint(1, 50)))
                                for j in range(10)]
        for s in sequences:
            ss = decoy_sequence(s, 'fused')
            self.assertEqual(ss, s[::-1] + 'R' + s)

    def test_decoy_keep_nterm(self):
        sequences = [''.join(random.choice(string.ascii_uppercase)
                             for i in range(random.randint(1, 50)))
                                for j in range(10)]
        for mode in ('shuffle', 'reverse'):
            for seq in sequences:
                self.assertEqual(seq[0], decoy_sequence(seq, mode, keep_nterm=True)[0])

    def test_read_and_write_fasta_short(self):
        with tempfile.TemporaryFile(mode='r+') as new_fasta_file:
            write(read(self.fasta_file, ignore_comments=True),
                new_fasta_file)
            new_fasta_file.seek(0)
            new_entries = list(read(new_fasta_file, ignore_comments=True))
            self.assertEqual(new_entries, self.fasta_entries_short)

    def test_read_and_write_long(self):
        with tempfile.TemporaryFile(mode='r+') as new_fasta_file:
            write(read(self.fasta_file), new_fasta_file)
            new_fasta_file.seek(0)
            new_entries = list(read(new_fasta_file))
            self.assertEqual(new_entries, self.fasta_entries_long)

    def test_write_decoy_db(self):
        with tempfile.TemporaryFile(mode='r+') as decdb:
            write_decoy_db(self.fasta_file, decdb,
                    decoy_only=False, prefix='PREFIX_')
            decdb.seek(0)
            all_entries = list(read(decdb, False))
        self.assertEqual(all_entries, self.fasta_entries_long +
                [('PREFIX_' + a, b[::-1]) for a, b in self.fasta_entries_long])

    def test_parser_uniprotkb_decoydb(self):
        header = ('sp|P27748|ACOX_RALEH Acetoin catabolism protein X OS=Ralstonia'
            ' eutropha (strain ATCC 17699 / H16 / DSM 428 / Stanier 337)'
            ' GN=acoX PE=4 SV=2')
        sequence = 'SEQUENCE'
        with tempfile.TemporaryFile(mode='r+') as db:
            write([(header, sequence)], db)
            db.seek(0)
            entries = list(decoy_db(db, prefix='PREFIX_', parser=parse, decoy_only=True))

        parsed = {'GN': 'acoX',
                 'OS': 'Ralstonia eutropha '
                    '(strain ATCC 17699 / H16 / DSM 428 / Stanier 337)',
                 'PE': 4,
                 'SV': 2,
                 'db': 'PREFIX_sp',
                 'entry': 'ACOX_RALEH',
                 'id': 'P27748',
                 'gene_id': 'ACOX',
                 'name': 'Acetoin catabolism protein X',
                 'taxon': 'RALEH'}
        self.assertEqual(entries[0][0], parsed)
        self.assertEqual(entries[0][1], 'SEQUENCE'[::-1])
        self.assertEqual(len(entries), 1)

    def test_parser_uniprotkb(self):
        header = ('sp|P27748|ACOX_RALEH Acetoin catabolism protein X OS=Ralstonia'
            ' eutropha (strain ATCC 17699 / H16 / DSM 428 / Stanier 337)'
            ' GN=acoX PE=4 SV=2')
        parsed = {'GN': 'acoX',
                 'OS': 'Ralstonia eutropha '
                    '(strain ATCC 17699 / H16 / DSM 428 / Stanier 337)',
                 'PE': 4,
                 'SV': 2,
                 'db': 'sp',
                 'entry': 'ACOX_RALEH',
                 'id': 'P27748',
                 'gene_id': 'ACOX',
                 'name': 'Acetoin catabolism protein X',
                 'taxon': 'RALEH'}
        self.assertEqual(parse(header), parsed)

    def test_parser_uniptokb_isoform(self):
        header = ('sp|Q4R572-2|1433B_MACFA Isoform Short of 14-3-3 protein beta'
                '/alpha OS=Macaca fascicularis GN=YWHAB')
        parsed = {'GN': 'YWHAB',
                 'OS': 'Macaca fascicularis',
                 'db': 'sp',
                 'entry': '1433B_MACFA',
                 'gene_id': '1433B',
                 'id': 'Q4R572-2',
                 'name': 'Isoform Short of 14-3-3 protein beta/alpha',
                 'taxon': 'MACFA'}
        self.assertEqual(parse(header), parsed)

    def test_parser_uniref(self):
        header = ('>UniRef100_A5DI11 Elongation factor 2 n=1 '
                'Tax=Pichia guilliermondii RepID=EF2_PICGU')
        parsed = {'RepID': 'EF2_PICGU',
                 'taxon': 'PICGU',
                 'gene_id': 'EF2',
                 'Tax': 'Pichia guilliermondii',
                 'cluster': 'Elongation factor 2',
                 'id': 'UniRef100_A5DI11',
                 'type': 'UniRef100',
                 'accession': 'A5DI11',
                 'n': 1}
        self.assertEqual(parse(header), parsed)

    def test_parser_uniparc(self):
        header = '>UPI0000000005 status=active'
        parsed = {'id': 'UPI0000000005', 'status': 'active'}
        self.assertEqual(parse(header), parsed)

    def test_parser_unimes(self):
        header = ('MES00000000005 Putative uncharacterized protein GOS_3018412 '
            '(Fragment) OS=marine metagenome Pep=JCVI_PEP_1096688850003 SV=1')
        parsed = {'OS': 'marine metagenome',
                 'Pep': 'JCVI_PEP_1096688850003',
                 'SV': 1,
                 'id': 'MES00000000005',
                 'name': 'Putative uncharacterized protein GOS_3018412 (Fragment)'}
        self.assertEqual(parse(header), parsed)

    def test_parser_spd(self):
        header = ('>P31947|1433S_HUMAN| 14-3-3 protein sigma (Stratifin) '
            '(Epithelial cell marker protein 1).')
        parsed = {'description': '14-3-3 protein sigma (Stratifin) '
                '(Epithelial cell marker protein 1).',
                 'gene': '1433S_HUMAN',
                 'gene_id': '1433S',
                 'id': 'P31947',
                 'taxon': 'HUMAN'}
        self.assertEqual(parse(header), parsed)

    def test_parser_spd_mult_ids(self):
        header = ('>P02763 Q8TC16|A1AG1_HUMAN| Alpha-1-acid glycoprotein 1 '
            'precursor (AGP 1) (Orosomucoid-1) (OMD 1)')
        parsed = {'description': 'Alpha-1-acid glycoprotein 1 precursor (AGP 1)'
                ' (Orosomucoid-1) (OMD 1)',
                 'gene': 'A1AG1_HUMAN',
                 'gene_id': 'A1AG1',
                 'id': 'P02763 Q8TC16',
                 'taxon': 'HUMAN'}
        self.assertEqual(parse(header), parsed)

if __name__ == '__main__':
    unittest.main()
