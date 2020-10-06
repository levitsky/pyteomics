from data import usi_proxi_data
from os import path
import pyteomics
pyteomics.__path__ = [path.abspath(path.join(path.dirname(__file__), path.pardir, 'pyteomics'))]

import unittest
from itertools import product
import operator as op

import numpy as np

from pyteomics.usi import USI, proxi


class USITest(unittest.TestCase):
    def test_parse(self):
        usi_str = "mzspec:MSV000085202:210320_SARS_CoV_2_T:scan:131256"
        inst = USI.parse(usi_str)
        assert str(inst) == usi_str
        assert inst.protocol == 'mzspec'
        assert inst.dataset == "MSV000085202"
        assert inst.datafile == "210320_SARS_CoV_2_T"
        assert inst.scan_identifier_type == "scan"
        assert inst.scan_identifier == "131256"
        assert inst.interpretation == None


class PROXITest(unittest.TestCase):
    def test_request(self):
        usi_str = "mzspec:MSV000085202:210320_SARS_CoV_2_T:scan:131256"
        response = proxi(usi_str, backend='peptide_atlas')

        assert usi_proxi_data.keys() <= response.keys()
        assert np.allclose(response['m/z array'] - usi_proxi_data['m/z array'], 0)
        assert np.allclose(response['intensity array'] - usi_proxi_data['intensity array'], 0)




if __name__ == "__main__":
    unittest.main()
