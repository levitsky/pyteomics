from os import path
import pyteomics
pyteomics.__path__ = [path.abspath(path.join(path.dirname(__file__), path.pardir, 'pyteomics'))]
import unittest
from pyteomics.mzid import MzIdentML, read, chain
from pyteomics import auxiliary as aux, mzid
from data import mzid_spectra
from itertools import product
import os
from psims.controlled_vocabulary.controlled_vocabulary import obo_cache


CACHE_PATH = os.environ.get('OBO_CACHE_PATH')
if CACHE_PATH:
    obo_cache.cache_path = CACHE_PATH

obo_cache.enabled = bool(os.environ.get("OBO_CACHE_ENABLED", True))


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

    def test_streaming_filter_preserves_psm_identity(self):
        score_name = 'ProteinScape:SequestMetaScore'

        def score(psm):
            return psm['SpectrumIdentificationItem'][0][score_name]

        with MzIdentML(self.path) as reader:
            psms = [p for p in reader if score_name in p['SpectrumIdentificationItem'][0]]
        key = 'spectrumID'
        expected = ['databasekey=1', 'databasekey=2', 'databasekey=11', 'databasekey=12']
        options = dict(key=score, is_decoy=mzid.is_decoy, reverse=True, fdr=0.01)
        self.assertCountEqual([psm[key] for psm in aux.filter(psms, **options)], expected)
        with aux.filter(psms, full_output=False, **options) as filtered:
            self.assertEqual([psm[key] for psm in filtered], expected)

    def test_unit_info(self):
        with MzIdentML(self.path) as handle:
            for protocol in handle.iterfind("SpectrumIdentificationProtocol"):
                fragment_tolerance = protocol['FragmentTolerance']
                self.assertEqual(fragment_tolerance['search tolerance minus value'].unit_info, 'dalton')
                parent_tolerance = protocol['ParentTolerance']
                self.assertEqual(parent_tolerance['search tolerance plus value'].unit_info, 'parts per million')

    def test_map(self):
        def key(x):
            return (x['spectrumID'], x['SpectrumIdentificationItem'][0]['PeptideSequence'], len(x['SpectrumIdentificationItem']))

        with MzIdentML(self.path) as reader:
            for method in ['t', 'p']:
                with self.subTest(method=method):
                    self.assertEqual(sorted(mzid_spectra[(1, 1)], key=key),
                                    sorted(reader.map(method=method), key=key))

    def test_iterfind_map(self):
        with MzIdentML(self.path) as r:
            for method in ['t', 'p']:
                with self.subTest(method=method):
                    self.assertEqual(
                        len(mzid_spectra[(1, 1)]),
                        sum(1 for _ in r.iterfind("SpectrumIdentificationResult").map(method=method)))

    def test_index(self):
        with MzIdentML(self.path) as reader:
            self.assertEqual(len(mzid_spectra[(1, 1)]), len(reader.index[MzIdentML._default_iter_tag]))


class ExplicitNSMzidTest(MzidTest):
    path = 'test.explicit-ns.mzid'


class OtherTest(unittest.TestCase):
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
