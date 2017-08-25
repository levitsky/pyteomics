import os
import shutil
from os import path
import tempfile
import pyteomics
pyteomics.__path__ = [path.abspath(path.join(path.dirname(__file__), path.pardir, 'pyteomics'))]
from itertools import product
import unittest
from pyteomics.mzml import *
from data import mzml_spectra
import numpy as np

class MzmlTest(unittest.TestCase):
    maxDiff = None
    path = 'test.mzML'

    def testReadSpectrum(self):
        for rs, it, ui in product([True, False], repeat=3):
            for func in [MzML, read, chain,
                    lambda x, **kw: chain.from_iterable([x], **kw), PreIndexedMzML]:
                with func(self.path, read_schema=rs, iterative=it, use_index=ui) as r:
                    # http://stackoverflow.com/q/14246983/1258041
                    self.assertEqual(mzml_spectra, list(r))

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
        # inst._source.close()
        self.assertTrue(inst._source.closed)
        MzML.prebuild_byte_offset_file(work_path)
        with MzML(work_path, use_index=True) as inst:
            offsets_exist = path.exists(inst._byte_offset_filename)
            self.assertTrue(offsets_exist)
            self.assertEqual(offsets_exist, inst._check_has_byte_offset_file())
            self.assertTrue(isinstance(inst._offset_index, xml.PrebuiltOffsetIndex))
        # inst._source.close()
        self.assertTrue(inst._source.closed)
        os.remove(inst._byte_offset_filename)
        with MzML(work_path, use_index=True) as inst:
            offsets_exist = path.exists(inst._byte_offset_filename)
            self.assertEqual(offsets_exist, inst._check_has_byte_offset_file())
            self.assertTrue(isinstance(inst._offset_index, xml.FlatTagSpecificXMLByteIndex))
            self.assertTrue(not isinstance(inst._offset_index, xml.PrebuiltOffsetIndex))
        # inst._source.close()
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

if __name__ == '__main__':
    unittest.main()
