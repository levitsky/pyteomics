from os import path
import pyteomics
pyteomics.__path__ = [path.abspath(path.join(path.dirname(__file__), path.pardir, 'pyteomics'))]
import unittest
from pyteomics.mzid import MzIdentML, read, chain
from pyteomics import auxiliary as aux
from data import mzid_spectra
from itertools import product

class MzidTest(unittest.TestCase):
    maxDiff = None
    path = 'test.mzid'
    def testReadPSM(self):
        for rec, refs, rs, it, ui in product((True, False), repeat=5):
            for func in [MzIdentML, read, chain,
                    lambda x, **kw: chain.from_iterable([x], **kw)]:
                with func(self.path, recursive=rec, retrieve_refs=refs,
                        read_schema=rs, iterative=it, use_index=ui) as reader:
                    psms = list(reader)
                    self.assertEqual(psms, mzid_spectra[(rec, refs)])

    def test_skip_empty_values(self):
        with MzIdentML(self.path, skip_empty_cvparam_values=True, recursive=True, retrieve_refs=True) as f:
            self.assertEqual(list(f), mzid_spectra[(True, True)])

    def test_unit_info(self):
        with MzIdentML(self.path) as handle:
            for protocol in handle.iterfind("SpectrumIdentificationProtocol"):
                fragment_tolerance = protocol['FragmentTolerance']
                self.assertEqual(fragment_tolerance['search tolerance minus value'].unit_info, 'dalton')
                parent_tolerance = protocol['ParentTolerance']
                self.assertEqual(parent_tolerance['search tolerance plus value'].unit_info, 'parts per million')

    def test_structure_normalization(self):
        gen = read('mzid_snippet.xml').iterfind("SpectraData")
        datum = next(gen)
        index = aux.cvquery(datum)
        assert index['MS:1000768'] == 'Thermo nativeID format'
        datum = next(gen)
        index = aux.cvquery(datum)
        assert index['MS:1000774'] == 'multiple peak list nativeID format'


if __name__ == '__main__':
    unittest.main()
