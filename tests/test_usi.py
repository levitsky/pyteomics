from data import usi_proxi_data
from os import path
import pyteomics
pyteomics.__path__ = [path.abspath(path.join(path.dirname(__file__), path.pardir, 'pyteomics'))]

import unittest

from pyteomics.usi import USI, proxi, AGGREGATOR_KEY
from pyteomics.auxiliary import PyteomicsError


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
        response = proxi(usi_str, backend='massive')

        assert set(usi_proxi_data.keys()) <= set(response.keys())

        for a, b in zip(response['m/z array'], usi_proxi_data['m/z array']):
            self.assertAlmostEqual(a, b, 3)

        for a, b in zip(response['intensity array'], usi_proxi_data['intensity array']):
            self.assertAlmostEqual(a, b, 3)

    def test_errors(self):
        usi_str = "mzspec:MSV000085202:210320_SARS_CoV_2_T:scan:131256"
        with self.assertRaises(TypeError, msg='Unrecognized backend type: NoneType'):
            proxi(usi_str, backend=None)
        with self.assertRaises(PyteomicsError, msg='Unknown PROXI backend name: BackendName'):
            proxi(usi_str, backend='BackendName')


class PROXIAggregatorTest(unittest.TestCase):
    def test_request(self):
        usi_str = "mzspec:MSV000085202:210320_SARS_CoV_2_T:scan:131256"
        response = proxi(usi_str, backend=AGGREGATOR_KEY)

        assert set(usi_proxi_data.keys()) <= set(response.keys())

        for a, b in zip(response['m/z array'], usi_proxi_data['m/z array']):
            self.assertAlmostEqual(a, b, 3)

        for a, b in zip(response['intensity array'], usi_proxi_data['intensity array']):
            self.assertAlmostEqual(a, b, 3)


if __name__ == "__main__":
    unittest.main()
