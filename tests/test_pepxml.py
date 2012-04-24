import unittest
from pyteomics.pepxml import *

import numpy as np

class PepxmlTest(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None

    def testReadPSM(self):
        psms = [psm for psm in iter_psm('test.pep.xml')]
        self.assertEqual(psms, [
            {'spectrum': 'pps_sl20060731_18mix_25ul_r1_1154456409.0100.0100.1',
             'end_scan': 100.0, 
             'start_scan': 100.0, 
             'index': 1.0,
             'assumed_charge': 1.0,
             'precursor_neutral_mass': 860.392,
             'search_hits': [{
             'num_missed_cleavages': 0.0, 
             'tot_num_ions': 12.0, 
             'is_rejected': '0',
             'deltacn': 0.081,
             'hit_rank': 1.0, 
             'num_matched_ions': 11.0, 
             'sprank': 1.0, 
             'num_tot_proteins': 1.0, 
             'peptide': 'SLNGEWR', 
             'massdiff': -0.5, 
             'peptideprophet': 0.96, 
             'deltacnstar': 0.0, 
             'spscore': 894.0, 
             'modifications': [],
             'modified_peptide': 'SLNGEWR', 
             'proteins': [{'num_tol_term': 2.0,
                           'protein': 'sp|P00722|BGAL_ECOLI',
                           'peptide_prev_aa': 'R',
                           'protein_descr': 'BETA-GALACTOSIDASE (EC 3.2.1.23) '
                                            '(LACTASE) - Escherichia coli.',
                           'peptide_next_aa': 'F'}], 
             'calc_neutral_pep_mass': 860.892, 
             'xcorr': 1.553}]},
            {'precursor_neutral_mass': 677.392,
             'spectrum': 'pps_sl20060731_18mix_25ul_r1_1154456409.0040.0040.1',
             'start_scan': 40.0, 
             'assumed_charge': 1.0, 
             'index': 2.0,
             'end_scan': 40.0,
             'search_hits': [{'tot_num_ions': 10.0,
             'num_missed_cleavages': 1.0,
             'is_rejected': '0',
             'hit_rank': 1.0,
             'num_matched_ions': 8.0,
             'sprank': 1.0,
             'num_tot_proteins': 1.0,
             'peptide': 'GKKFAK', 
             'massdiff': -0.5,
             'peptideprophet': 0.548,
             'deltacn': 0.165,
             'deltacnstar': 0.0,
             'spscore': 427.0,
             'modifications': [],
             'modified_peptide': 'GKKFAK', 
             'proteins': [{'num_tol_term': 1.0,
                           'protein': 'gi|3212198|gb|AAC22319.1|', 
                           'peptide_prev_aa': 'N', 
                           'protein_descr': 'hemoglobin-binding protein '
                                            '[Haemophilus influenzae Rd]',
                           'peptide_next_aa': 'I'}],
             'calc_neutral_pep_mass': 677.892,
             'xcorr': 1.644}]}])

if __name__ == '__main__':
    unittest.main()
