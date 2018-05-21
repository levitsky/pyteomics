import os
import shutil
from os import path
import tempfile
import pyteomics
pyteomics.__path__ = [path.abspath(path.join(path.dirname(__file__), path.pardir, 'pyteomics'))]
from itertools import product
import unittest
from pyteomics.mzml import MzML, PreIndexedMzML, read, chain
from pyteomics import auxiliary as aux, xml
from data import mzml_spectra, mzml_spectra_skip_empty_values
import numpy as np

class MzmlTest(unittest.TestCase):
    maxDiff = None
    path = 'test.mzML'

    def test_read(self):
        for rs, it, ui in product([True, False], repeat=3):
            for func in [MzML, read, chain,
                    lambda x, **kw: chain.from_iterable([x], **kw), PreIndexedMzML]:
                with func(self.path, read_schema=rs, iterative=it, use_index=ui) as r:
                    # http://stackoverflow.com/q/14246983/1258041
                    self.assertEqual(mzml_spectra, list(r))

    def test_read_skip_empty_values(self):
        with MzML(self.path, skip_empty_cvparam_values=True) as r:
            self.assertEqual(mzml_spectra_skip_empty_values, list(r))

    def test_decoding(self):
        with MzML(self.path, decode_binary=True) as reader:
            spectrum = next(reader)
            self.assertIsNotNone(spectrum['m/z array'])
        validation = spectrum['m/z array']
        with MzML(self.path) as reader:
            spectrum = next(reader)
            self.assertIsNotNone(spectrum['m/z array'])
            self.assertTrue(np.allclose(spectrum['m/z array'], validation))
        with MzML(self.path, decode_binary=False) as reader:
            spectrum = next(reader)
            record = spectrum['m/z array']
            self.assertEqual(record.compression, "no compression")
            self.assertEqual(record.dtype, "d")
            array = record.decode()
            self.assertTrue(np.allclose(validation, array))
            record = spectrum['intensity array']
            self.assertEqual(record.dtype, "f")
            self.assertEqual(record.compression, "no compression")

            spectrum = next(reader)
            record = spectrum['intensity array']
            self.assertEqual(record.compression, "zlib compression")
            self.assertEqual(mzml_spectra[1]['intensity array'], record.decode())

    def test_read_dtype(self):
        dtypes = {'m/z array': np.float32, 'intensity array': np.int32}
        with read(self.path, dtype=dtypes) as f:
            for spec in f:
                for k, v in dtypes.items():
                    self.assertEqual(spec[k].dtype, v)

    def test_has_built_index(self):
        with read(self.path, use_index=True) as f:
            self.assertGreater(len(f._offset_index), 0)
        with read(self.path, use_index=False) as f:
            self.assertEqual(len(f._offset_index), 0)

    def test_prebuild_index(self):
        test_dir = tempfile.mkdtemp()
        work_path = path.join(test_dir, self.path)
        with open(work_path, 'w') as dest, open(self.path) as source:
            dest.write(source.read())
        assert dest.closed
        with MzML(work_path, use_index=True) as inst:
            offsets_exist = path.exists(inst._byte_offset_filename)
            self.assertEqual(offsets_exist, inst._check_has_byte_offset_file())
            self.assertTrue(isinstance(inst._offset_index, xml.FlatTagSpecificXMLByteIndex))
            self.assertTrue(not isinstance(inst._offset_index, xml.PrebuiltOffsetIndex))
        self.assertTrue(inst._source.closed)
        MzML.prebuild_byte_offset_file(work_path)
        with MzML(work_path, use_index=True) as inst:
            offsets_exist = path.exists(inst._byte_offset_filename)
            self.assertTrue(offsets_exist)
            self.assertEqual(offsets_exist, inst._check_has_byte_offset_file())
            self.assertTrue(isinstance(inst._offset_index, xml.PrebuiltOffsetIndex))
        self.assertTrue(inst._source.closed)
        os.remove(inst._byte_offset_filename)
        with MzML(work_path, use_index=True) as inst:
            offsets_exist = path.exists(inst._byte_offset_filename)
            self.assertEqual(offsets_exist, inst._check_has_byte_offset_file())
            self.assertTrue(isinstance(inst._offset_index, xml.FlatTagSpecificXMLByteIndex))
            self.assertTrue(not isinstance(inst._offset_index, xml.PrebuiltOffsetIndex))
        self.assertTrue(inst._source.closed)
        shutil.rmtree(test_dir, True)

    def test_unit_extract(self):
        with MzML(self.path) as handle:
            for scan in handle:
                scan_inst = scan['scanList']['scan'][0]
                scan_time = scan_inst['scan start time']
                scan_window_lower_limit = scan_inst['scanWindowList']['scanWindow'][0]['scan window lower limit']
                self.assertEqual(scan_time.unit_info, 'minute')
                self.assertEqual(scan_window_lower_limit.unit_info, 'm/z')

    def test_cv_query(self):
        with MzML(self.path) as handle:
            scan = next(handle)
            index = aux.cvquery(scan)
            self.assertEqual(index['MS:1000511'], 1)
            self.assertEqual(aux.cvquery(scan, "MS:1000511"), 1)

            # test deep traversal
            self.assertEqual(index['MS:1000016'], 0.004935)
            self.assertEqual(aux.cvquery(scan, 'MS:1000016'), 0.004935)

    def test_retrieve_refs(self):
        with MzML(self.path) as reader:
            derefed = list(reader.iterfind("instrumentConfiguration", retrieve_refs=True))
            reader.reset()
            raw = list(reader.iterfind("instrumentConfiguration", retrieve_refs=False))
            self.assertEqual(raw[0].get("ref"), 'CommonInstrumentParams')
            self.assertNotIn("ref", derefed[0])
            self.assertEqual(derefed[0].get('instrument serial number'), 'SN06061F')


if __name__ == '__main__':
    unittest.main()
