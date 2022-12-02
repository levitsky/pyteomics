import os
import shutil
import tempfile
import pyteomics
from io import BytesIO
pyteomics.__path__ = [os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir, 'pyteomics'))]
from itertools import product
import unittest
from pyteomics.mzml import MzML, PreIndexedMzML, read, chain
from pyteomics import auxiliary as aux, xml
from data import mzml_spectra
import numpy as np
import pickle
import operator as op
import pynumpress
import base64
import zlib

class MzmlTest(unittest.TestCase):
    maxDiff = None
    path = 'test.mzML'

    def test_read(self):
        for rs, it, ui in product([True, False], repeat=3):
            if rs: continue # temporarily disable retrieval of schema
            for func in [MzML, read, chain,
                    lambda x, **kw: chain.from_iterable([x], **kw), PreIndexedMzML]:
                with func(self.path, read_schema=rs, iterative=it, use_index=ui) as r:
                    # http://stackoverflow.com/q/14246983/1258041
                    self.assertEqual(mzml_spectra, list(r))

    def test_mp_read(self):
        key = op.itemgetter('index')
        with MzML(self.path) as f:
            self.assertEqual(sorted(mzml_spectra, key=key), sorted(list(f.map()), key=key))

    def test_mp_requires_index(self):
        with MzML(self.path, use_index=False) as r:
            self.assertRaises(aux.PyteomicsError, r.map)

    def test_map_qsize(self):
        key = op.itemgetter('index')
        with MzML(self.path, queue_size=1000) as f:
            self.assertEqual(f._queue_size, 1000)
            self.assertEqual(sorted(mzml_spectra, key=key), sorted(list(f.map()), key=key))

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
            self.assertEqual(record.dtype, np.float64)
            array = record.decode()
            self.assertTrue(np.allclose(validation, array))
            record = spectrum['intensity array']
            self.assertEqual(record.dtype, np.float32)
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
            self.assertEqual(f._offset_index, None)

    def test_prebuild_index(self):
        test_dir = tempfile.mkdtemp()
        work_path = os.path.join(test_dir, self.path)
        with open(work_path, 'w') as dest, open(self.path) as source:
            dest.write(source.read())
        assert dest.closed
        with MzML(work_path, use_index=False) as inst:
            self.assertRaises(IOError, inst._read_byte_offsets)
            with open(inst._byte_offset_filename, 'wt') as fh:
                fh.write("{}")
            self.assertRaises(TypeError, inst._read_byte_offsets)
            os.remove(inst._byte_offset_filename)
        with MzML(work_path, use_index=True) as inst:
            offsets_exist = os.path.exists(inst._byte_offset_filename)
            self.assertEqual(offsets_exist, inst._check_has_byte_offset_file())
            self.assertTrue(isinstance(inst._offset_index, xml.HierarchicalOffsetIndex))
        self.assertTrue(inst._source.closed)
        MzML.prebuild_byte_offset_file(work_path)
        with open(inst._byte_offset_filename, 'rt') as fh:
            index = MzML._index_class.load(fh)
            assert inst._offset_index['spectrum'] == index['spectrum']
        with MzML(work_path, use_index=True) as inst:
            offsets_exist = os.path.exists(inst._byte_offset_filename)
            self.assertTrue(offsets_exist)
            self.assertEqual(offsets_exist, inst._check_has_byte_offset_file())
            self.assertTrue(isinstance(inst._offset_index, xml.HierarchicalOffsetIndex))
        self.assertTrue(inst._source.closed)
        os.remove(inst._byte_offset_filename)
        with MzML(work_path, use_index=True) as inst:
            offsets_exist = os.path.exists(inst._byte_offset_filename)
            self.assertEqual(offsets_exist, inst._check_has_byte_offset_file())
            self.assertTrue(isinstance(inst._offset_index, xml.HierarchicalOffsetIndex))
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
            self.assertEqual(raw[0].get("softwareRef"),
                             {'ref': 'Xcalibur'})
            self.assertNotIn("ref", derefed[0]['softwareRef'])
            self.assertEqual(derefed[0].get('softwareRef'), {
                             'version': '1.1 Beta 7', 'Xcalibur': ''})

    def test_in_memory_buffer(self):
        with open(self.path, 'rb') as fh:
            data_buffer = BytesIO(fh.read())
        with MzML(data_buffer) as reader:
            spectrum = next(reader)
            self.assertEqual(spectrum['id'], 'controllerType=0 controllerNumber=1 scan=1')
        data_buffer.seek(0)
        with MzML(data_buffer, use_index=True) as reader:
            spectrum = next(reader)
            self.assertEqual(spectrum['id'], 'controllerType=0 controllerNumber=1 scan=1')

    def test_picklable(self):
        with MzML(self.path) as reader:
            expected_data = next(reader)
            spec = pickle.dumps(reader)
        with pickle.loads(spec) as reader:
            self.assertEqual(next(reader)['id'], expected_data['id'])

    def test_indexing(self):
        with MzML(self.path) as reader:
            self.assertEqual(mzml_spectra[0], reader[0])
            self.assertEqual(mzml_spectra[0], reader['controllerType=0 controllerNumber=1 scan=1'])
            self.assertEqual(mzml_spectra, reader[0:2])
            self.assertEqual(mzml_spectra,
                [reader['controllerType=0 controllerNumber=1 scan=1'],
                 reader['controllerType=0 controllerNumber=1 scan=2']])
            self.assertEqual(mzml_spectra, reader[[0, 1]])
            self.assertEqual(mzml_spectra, reader[
                ['controllerType=0 controllerNumber=1 scan=1', 'controllerType=0 controllerNumber=1 scan=2']])
            self.assertEqual(mzml_spectra, reader[
                'controllerType=0 controllerNumber=1 scan=2':'controllerType=0 controllerNumber=1 scan=1'])

    def test_time_locator(self):
        with MzML(self.path) as reader:
            self.assertEqual(mzml_spectra[0], reader.time[0])
            self.assertEqual(mzml_spectra[1], reader.time[0.1])
            self.assertEqual(mzml_spectra, reader.time[0:0.1])

    def test_numpress_slof(self):
        data = mzml_spectra[0]['intensity array']
        encoded = base64.b64encode(pynumpress.encode_slof(data, pynumpress.optimal_slof_fixed_point(data)).tobytes()).decode('ascii')
        record = aux.BinaryDataArrayTransformer()._make_record(encoded, 'MS-Numpress short logged float compression', data.dtype)
        self.assertTrue(np.allclose(data, record.decode(), rtol=0.001))

    def test_numpress_slof_zlib(self):
        data = mzml_spectra[0]['intensity array']
        encoded = base64.b64encode(zlib.compress(pynumpress.encode_slof(data, pynumpress.optimal_slof_fixed_point(data)).tobytes())).decode('ascii')
        record = aux.BinaryDataArrayTransformer()._make_record(encoded, 'MS-Numpress short logged float compression followed by zlib compression', data.dtype)
        self.assertTrue(np.allclose(data, record.decode(), rtol=0.001))

    def test_numpress_linear(self):
        data = mzml_spectra[0]['intensity array']
        encoded = base64.b64encode(pynumpress.encode_linear(data, pynumpress.optimal_linear_fixed_point(data)).tobytes()).decode('ascii')
        record = aux.BinaryDataArrayTransformer()._make_record(encoded, 'MS-Numpress linear prediction compression', data.dtype)
        self.assertTrue(np.allclose(data, record.decode(), rtol=0.001))

    def test_numpress_linear_zlib(self):
        data = mzml_spectra[0]['intensity array']
        encoded = base64.b64encode(zlib.compress(pynumpress.encode_linear(data, pynumpress.optimal_linear_fixed_point(data)).tobytes())).decode('ascii')
        record = aux.BinaryDataArrayTransformer()._make_record(encoded, 'MS-Numpress linear prediction compression followed by zlib compression', data.dtype)
        self.assertTrue(np.allclose(data, record.decode(), rtol=0.001))

    def test_numpress_pic(self):
        data = mzml_spectra[0]['intensity array']
        encoded = base64.b64encode(pynumpress.encode_pic(data).tobytes()).decode('ascii')
        record = aux.BinaryDataArrayTransformer()._make_record(encoded, 'MS-Numpress positive integer compression', data.dtype)
        self.assertTrue(np.allclose(data, record.decode(), atol=0.6))

    def test_numpress_pic_zlib(self):
        data = mzml_spectra[0]['intensity array']
        encoded = base64.b64encode(zlib.compress(pynumpress.encode_pic(data).tobytes())).decode('ascii')
        record = aux.BinaryDataArrayTransformer()._make_record(encoded, 'MS-Numpress positive integer compression followed by zlib compression', data.dtype)
        self.assertTrue(np.allclose(data, record.decode(), atol=0.6))


if __name__ == '__main__':
    unittest.main()
