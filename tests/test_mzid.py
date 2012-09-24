import unittest
from pyteomics.mzid import *

class MzidTest(unittest.TestCase):
    def testReadPSM(self):
        psms = [psm for psm in read('test.mzid')][0:4]
        self.assertEqual(psms, [
           {'SpectrumIdentificationItem': [
               {'ProteinScape:IntensityCoverage': 0.3919545603809718, 
                'PeptideEvidenceRef': [{'peptideEvidence_ref': 'PE1_SEQ_spec1_pep1'}],
                'passThreshold': True, 
                'rank': 1, 'chargeState': 1, 
                'calculatedMassToCharge': 1507.695, 
                'peptide_ref': 'prot1_pep1', 
                'experimentalMassToCharge': 1507.696, 
                'id': 'SEQ_spec1_pep1', 
                'ProteinScape:SequestMetaScore': 7.59488518903425}], 
            'spectrumID': 'databasekey=1', 
            'id': 'SEQ_spec1', 
            'spectraData_ref': 'LCMALDI_spectra'},
 
            {'SpectrumIdentificationItem': [
                {'ProteinScape:IntensityCoverage': 0.5070386909133888, 
                'PeptideEvidenceRef': [{'peptideEvidence_ref': 'PE1_SEQ_spec2a_pep1'}], 
                'passThreshold': True, 
                'rank': 1, 
                'chargeState': 1, 
                'calculatedMassToCharge': 1920.9224, 
                'peptide_ref': 'prot1_pep2', 
                'experimentalMassToCharge': 1920.923, 
                'id': 'SEQ_spec2a_pep1', 'ProteinScape:SequestMetaScore': 10.8810331335713}], 
            'spectrumID': 'databasekey=2', 
            'id': 'SEQ_spec2a', 
            'spectraData_ref': 'LCMALDI_spectra'}, 

            {'SpectrumIdentificationItem': [
                {'ProteinScape:IntensityCoverage': 0.43376827663349576, 
                'PeptideEvidenceRef': [{'peptideEvidence_ref': 'PE1_SEQ_spec3a_pep1'}], 
                'passThreshold': True, 
                'rank': 1, 
                'chargeState': 1, 
                'calculatedMassToCharge': 864.4752, 
                'peptide_ref': 'prot1_pep3', 
                'experimentalMassToCharge': 864.474, 
                'id': 'SEQ_spec3a_pep1', 
                'ProteinScape:SequestMetaScore': 6.1021771936508955}], 
            'spectrumID': 'databasekey=3', 
            'id': 'SEQ_spec3a', 
            'spectraData_ref': 'LCMALDI_spectra'}, 

            {'SpectrumIdentificationItem': 
                [{'ProteinScape:IntensityCoverage': 0.16164593872706742, 
                'PeptideEvidenceRef': [{'peptideEvidence_ref': 'PE1_SEQ_spec10_pep1'}], 
                'passThreshold': True, 
                'rank': 1, 
                'chargeState': 1, 
                'calculatedMassToCharge': 1832.862115, 
                'peptide_ref': 'prot1_pep4', 
                'experimentalMassToCharge': 1832.863, 
                'id': 'SEQ_spec10_pep1', 'ProteinScape:SequestMetaScore': 5.635013787097159}], 
                'spectrumID': 'databasekey=10', 
                'id': 'SEQ_spec10', 
                'spectraData_ref': 'LCMALDI_spectra'}])

if __name__ == '__main__':
    unittest.main()
