from os import path
import pyteomics
pyteomics.__path__ = [path.abspath(path.join(path.dirname(__file__), path.pardir, 'pyteomics'))]
import unittest
from pyteomics.mzid import MzIdentML, read, chain
from pyteomics import auxiliary as aux
from data import mzid_spectra
from itertools import product
import operator as op
from psims.controlled_vocabulary.controlled_vocabulary import obo_cache
obo_cache.cache_path = '.'
obo_cache.enabled = True


class MzidTest(unittest.TestCase):
    maxDiff = None
    path = 'test.mzid'

    def test_read(self):
        for rec, refs, rs, it, ui in product((True, False), repeat=5):
            for func in [MzIdentML, read, chain, lambda x, **kw: chain.from_iterable([x], **kw)]:
                with self.subTest(rec=rec, refs=refs, rs=rs, it=it, ui=ui, func=func):
                    with func(self.path, recursive=rec, retrieve_refs=refs, read_schema=rs, iterative=it, use_index=ui) as reader:
                        psms = list(reader)
                        self.assertEqual(psms, mzid_spectra[(rec, refs)])

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

    def test_map(self):
        def key(x):
            return (x['spectrumID'], x['SpectrumIdentificationItem'][0]['PeptideSequence'], len(x['SpectrumIdentificationItem']))

        with MzIdentML(self.path) as reader:
            reader.chunksize = 1
            for method in ['t', 'p']:
                with self.subTest(method=method):
                    self.assertEqual(sorted(mzid_spectra[(1, 1)], key=key),
                                    sorted(reader.map(method=method), key=key))

    def test_iterfind_map(self):
        with MzIdentML(self.path) as r:
            self.assertEqual(
                len(mzid_spectra[(1, 1)]),
                sum(1 for _ in r.iterfind("SpectrumIdentificationResult").map()))


if __name__ == '__main__':
    unittest.main()
