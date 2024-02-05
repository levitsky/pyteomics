import os
import numpy as np
import pyteomics

pyteomics.__path__ = [os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir, 'pyteomics'))]

import tempfile
import unittest
import pickle
import shutil
import json
from collections import OrderedDict
import warnings
from pyteomics import mgf, auxiliary as aux
import data


class MGFTest(unittest.TestCase):
    maxDiff = None
    _encoding = 'utf-8'

    def setUp(self):
        self.path = 'test.mgf'
        self.header = mgf.read_header(self.path)
        with mgf.read(self.path) as f:
            self.spectra = list(f)
        self.tmpfile = tempfile.TemporaryFile(mode='r+')
        mgf.write(header=self.header, spectra=self.spectra, output=self.tmpfile)
        self.tmpfile.seek(0)
        self.header2 = mgf.read_header(self.tmpfile)
        self.tmpfile.seek(0)
        tmpreader = mgf.read(self.tmpfile)
        self.spectra2 = list(tmpreader)
        self.ns = len(self.spectra)
        self.tmpfile.close()
        self.path_annotated = 'test_annotated.mgf'
        self.header_annotated = mgf.read_header(self.path_annotated)
        with mgf.read(self.path_annotated, read_ions=True) as f:
            self.spectra_annotated = list(f)

    def test_read(self):
        for func in [mgf.read, mgf.MGF, mgf.IndexedMGF]:
            # http://stackoverflow.com/q/14246983/1258041
            self.assertEqual(data.mgf_spectra_long, list(func(self.path)))
            self.assertEqual(data.mgf_spectra_short, list(func(self.path, False)))
            with func(self.path) as reader:
                self.assertEqual(data.mgf_spectra_long, list(reader))
            with func(self.path, False) as reader:
                self.assertEqual(data.mgf_spectra_short, list(reader))

    def test_read_source_kw(self):
        for func in [mgf.read, mgf.MGF, mgf.IndexedMGF]:
            self.assertEqual(data.mgf_spectra_long, list(func(source=self.path)))

    def test_read_decoding(self):
        for func in [mgf.read, mgf.MGF, mgf.IndexedMGF]:
            self.assertEqual(data.mgf_spectra_long_decoded,
                             list(func(self.path, encoding=self._encoding)))
            self.assertEqual(data.mgf_spectra_short_decoded,
                             list(func(self.path, False, encoding=self._encoding)))
            with func(self.path, encoding=self._encoding) as reader:
                self.assertEqual(data.mgf_spectra_long_decoded, list(reader))
            with func(self.path, False, encoding=self._encoding) as reader:
                self.assertEqual(data.mgf_spectra_short_decoded, list(reader))
            self.assertEqual(data.mgf_spectra_long_decoded, list(func(self.path)))

    def test_read_no_charges(self):
        with mgf.read(self.path, read_charges=False) as reader:
            self.assertEqual(data.mgf_spectra_long_no_charges, list(reader))
        with mgf.read(self.path, False, read_charges=False) as reader:
            self.assertEqual(data.mgf_spectra_short_no_charges, list(reader))

    def test_read_with_ions(self):
        for spec_data, spec_read in zip(data.mgf_spectra_annotated_long, list(self.spectra_annotated)):
            # Check that the spectra have the same dict keys
            self.assertEqual(spec_data.keys(), spec_read.keys())
            for key in spec_data.keys():
                if type(spec_data[key]) == dict:
                    self.assertDictEqual(spec_data[key], spec_read[key])
                else:
                    np.testing.assert_array_equal(spec_data[key], spec_read[key])

    def test_read_write_with_ions(self):
        formats = ['{:.6f} {:.6f} {}', '%.6f %.6f %s']
        for use_numpy in range(2):
            with tempfile.TemporaryFile(mode='r+') as f:
                mgf.write(self.spectra_annotated, f, write_ions=True, use_numpy=use_numpy,
                    fragment_format=formats[use_numpy])
                f.seek(0)
                spectra = list(mgf.read(f, read_ions=True))
            for spec_data, spec_read in zip(data.mgf_spectra_annotated_long, spectra):
                # Check that the spectra have the same dict keys
                self.assertEqual(spec_data.keys(), spec_read.keys())
                for key in spec_data.keys():
                    if type(spec_data[key]) == dict:
                        self.assertDictEqual(spec_data[key], spec_read[key])
                    else:
                        np.testing.assert_array_equal(spec_data[key], spec_read[key])

    def test_read_array_conversion(self):
        with mgf.read(self.path, convert_arrays=0) as reader:
            self.assertEqual(data.mgf_spectra_lists, list(reader))
        with mgf.read(self.path, convert_arrays=2) as reader:
            s = next(reader)
            self.assertTrue(isinstance(s['charge array'], np.ma.core.MaskedArray))
            self.assertTrue(isinstance(s['m/z array'], np.ndarray))
        with mgf.read(self.path, convert_arrays=1) as reader:
            s = next(reader)
            self.assertTrue(isinstance(s['charge array'], np.ndarray))
            self.assertTrue(isinstance(s['m/z array'], np.ndarray))

    def test_header(self):
        self.assertEqual(self.header, self.header2)

    def test_readwrite_ns(self):
        self.assertEqual(self.ns, len(self.spectra2))

    def test_readwrite_keys(self):
        for s, s2 in zip(self.spectra, self.spectra2):
            self.assertEqual(set(s), set(s2))
            self.assertEqual(set(s), {'intensity array', 'm/z array', 'params', 'charge array'})

    def test_readwrite_params(self):
        for s, s2 in zip(self.spectra, self.spectra2):
            self.assertEqual(s['params'], s2['params'])

    def test_readwrite_msms_len(self):
        for i in range(self.ns):
            al = len(self.spectra[i]['m/z array'])
            self.assertEqual(al, len(self.spectra[i]['intensity array']))
            self.assertEqual(al, len(self.spectra2[i]['m/z array']))
            self.assertEqual(al, len(self.spectra2[i]['intensity array']))
            for j in range(al):
                self.assertEqual(self.spectra[i]['m/z array'][j],
                                 self.spectra2[i]['m/z array'][j])
                self.assertEqual(self.spectra[i]['intensity array'][j],
                                 self.spectra2[i]['intensity array'][j])

    def test_readwrite_msms(self):
        for i in range(self.ns):
            al = len(self.spectra[i]['m/z array'])
            for j in range(al):
                self.assertEqual(self.spectra[i]['m/z array'][j],
                                 self.spectra2[i]['m/z array'][j])
                self.assertEqual(self.spectra[i]['intensity array'][j],
                                 self.spectra2[i]['intensity array'][j])

    def test_write_single(self):
        tmpfile = tempfile.TemporaryFile(mode='r+')
        with warnings.catch_warnings(record=True) as ws:
            warnings.simplefilter("always")
            for spectrum in self.spectra:
                mgf.write(spectra=spectrum, output=tmpfile)

        self.assertGreaterEqual(len(ws), 2)
        n_warned = 0
        for w in ws:
            n_warned += (issubclass(w.category, UserWarning) and "discouraged" in str(w.message))
        self.assertGreaterEqual(n_warned, 2)
        tmpfile.seek(0)
        tmpreader = mgf.read(tmpfile)
        self.assertEqual(data.mgf_spectra_long, list(tmpreader))

    def test_read_dtype(self):
        dtypes = {'m/z array': np.float32, 'intensity array': np.int32}
        with mgf.read(self.path, dtype=dtypes) as f:
            for spec in f:
                for k, v in dtypes.items():
                    self.assertEqual(spec[k].dtype, v)

    def test_get_spectrum(self):
        key = 'Spectrum 2'
        for klass in [mgf.MGF, mgf.IndexedMGF]:
            f = klass(self.path)
            self.assertEqual(data.mgf_spectra_long[1], f[key])
            self.assertEqual(data.mgf_spectra_long[1], f.get_spectrum(key))
        self.assertEqual(data.mgf_spectra_long[1], mgf.get_spectrum(self.path, key))

    def test_key_access_ions(self):
        with mgf.IndexedMGF(self.path_annotated, read_ions=True) as f:
            np.testing.assert_array_equal(f['RAEYWENYPPAH||3']['ion array'], self.spectra_annotated[1]['ion array'])

    def test_read_list(self):
        key = ['Spectrum 2', 'Spectrum 1']
        with mgf.IndexedMGF(self.path) as f:
            self.assertEqual(data.mgf_spectra_long[::-1], f[key])

    def test_indexedmgf_picklable(self):
        with mgf.IndexedMGF(self.path, block_size=12345) as reader:
            spec = pickle.dumps(reader)
        with pickle.loads(spec) as reader:
            self.assertEqual(data.mgf_spectra_long[0], next(reader))
            self.assertEqual(reader.block_size, 12345)

    def test_mgf_picklable(self):
        with mgf.MGF(self.path, convert_arrays=0) as reader:
            spec = pickle.dumps(reader)
        with pickle.loads(spec) as reader:
            self.assertEqual(data.mgf_spectra_lists[0], next(reader))

    def test_map(self):
        with mgf.IndexedMGF(self.path) as reader:
            spectra = sorted(list(reader.map()), key=lambda s: s['params']['title'])
        self.assertEqual(data.mgf_spectra_long, spectra)

    def test_prebuild_index(self):
        test_dir = tempfile.mkdtemp()
        work_path = os.path.join(test_dir, self.path)
        with open(work_path, 'w') as dest, open(self.path) as source:
            dest.write(source.read())
        assert dest.closed
        with mgf.IndexedMGF(work_path) as inst:
            offsets_exist = os.path.exists(inst._byte_offset_filename)
            self.assertEqual(offsets_exist, inst._check_has_byte_offset_file())
            self.assertTrue(isinstance(inst._offset_index, aux.OffsetIndex))
        self.assertTrue(inst._source.closed)
        with mgf.IndexedMGF(work_path) as inst:
            inst._offset_index.pop('Spectrum 1')
            inst.write_byte_offsets()
        with mgf.IndexedMGF(work_path) as inst:
            offsets_exist = os.path.exists(inst._byte_offset_filename)
            self.assertTrue(offsets_exist)
            self.assertEqual(offsets_exist, inst._check_has_byte_offset_file())
            self.assertTrue(isinstance(inst._offset_index, aux.OffsetIndex))
            self.assertEqual(len(inst), 1)
        self.assertTrue(inst._source.closed)
        os.remove(inst._byte_offset_filename)
        with mgf.IndexedMGF(work_path) as inst:
            offsets_exist = os.path.exists(inst._byte_offset_filename)
            self.assertEqual(offsets_exist, inst._check_has_byte_offset_file())
            self.assertTrue(isinstance(inst._offset_index, aux.OffsetIndex))
        self.assertTrue(inst._source.closed)
        shutil.rmtree(test_dir, True)

    def test_write_index_keys(self):
        test_dir = tempfile.mkdtemp()
        work_path = os.path.join(test_dir, self.path)
        with open(work_path, 'wb') as dest, open(self.path, 'rb') as source:
            dest.write(source.read())
        assert dest.closed
        mgf.IndexedMGF.prebuild_byte_offset_file(work_path)
        with mgf.IndexedMGF(work_path) as inst:
            ipath = inst._byte_offset_filename
        with open(ipath) as ifp:
            container = json.load(ifp, object_hook=OrderedDict)
        tag_key = mgf.IndexedMGF._index_class._schema_version_tag_key
        self.assertEqual(set(container.keys()), {tag_key, 'index'})
        self.assertEqual(tuple(container[tag_key]), mgf.IndexedMGF._index_class.schema_version)
        self.assertEqual(container['index'], [['Spectrum 1', [217, 343]], ['Spectrum 2', [343, 504]]])


class UtilityTest(unittest.TestCase):
    def test_charge_repr_single(self):
        self.assertEqual(mgf._charge_repr('charge', 2), 'CHARGE=2+')
        self.assertEqual(mgf._charge_repr('charge', '2'), 'CHARGE=2+')
        self.assertEqual(mgf._charge_repr('charge', [2]), 'CHARGE=2+')
        self.assertEqual(mgf._charge_repr('charge', aux.Charge(2)), 'CHARGE=2+')
        self.assertEqual(mgf._charge_repr('charge', aux.ChargeList([2])), 'CHARGE=2+')
        self.assertEqual(mgf._charge_repr('charge', np.int64(2)), 'CHARGE=2+')

    def test_charge_repr_multiple(self):
        self.assertEqual(mgf._charge_repr('charge', [2, 3]), 'CHARGE=2+ and 3+')
        self.assertEqual(mgf._charge_repr('charge', aux.ChargeList([2, 3])), 'CHARGE=2+ and 3+')
        self.assertEqual(mgf._charge_repr('charge', '2+, 3+'), 'CHARGE=2+ and 3+')
        self.assertEqual(mgf._charge_repr('charge', np.array([2, 3])), 'CHARGE=2+ and 3+')


if __name__ == "__main__":
    unittest.main()
