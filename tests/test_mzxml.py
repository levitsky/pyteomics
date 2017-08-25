import os
from os import path
import pyteomics
pyteomics.__path__ = [path.abspath(path.join(path.dirname(__file__), path.pardir, 'pyteomics'))]
from itertools import product
import unittest
from pyteomics.mzxml import *
from data import mzxml_spectra
import tempfile
import shutil
import numpy as np

class MzXMLTest(unittest.TestCase):
    maxDiff = None
    path = 'test.mzXML'

    def testReadSpectrum(self):
        for rs, it, ui in product([True, False], repeat=3):
            for func in [MzXML, read, chain,
                    lambda x, **kw: chain.from_iterable([x], **kw)]:
                with func(self.path, read_schema=rs, iterative=it, use_index=ui) as r:
                    # http://stackoverflow.com/q/14246983/1258041
                    self.assertEqual(mzxml_spectra, list(r))

    def test_prebuild_index(self):
        test_dir = tempfile.mkdtemp()
        work_path = path.join(test_dir, self.path)
        with open(work_path, 'w') as dest, open(self.path) as source:
            dest.write(source.read())
        assert dest.closed
        with MzXML(work_path, use_index=True) as inst:
            offsets_exist = path.exists(inst._byte_offset_filename)
            self.assertEqual(offsets_exist, inst._check_has_byte_offset_file())
            self.assertTrue(isinstance(inst._offset_index, xml.FlatTagSpecificXMLByteIndex))
            self.assertTrue(not isinstance(inst._offset_index, xml.PrebuiltOffsetIndex))
        inst._source.close()
        self.assertTrue(inst._source.closed)
        MzXML.prebuild_byte_offset_file(work_path)
        with MzXML(work_path, use_index=True) as inst:
            offsets_exist = path.exists(inst._byte_offset_filename)
            self.assertTrue(offsets_exist)
            self.assertEqual(offsets_exist, inst._check_has_byte_offset_file())
            self.assertTrue(isinstance(inst._offset_index, xml.PrebuiltOffsetIndex))
        inst._source.close()
        self.assertTrue(inst._source.closed)
        os.remove(inst._byte_offset_filename)
        with MzXML(work_path, use_index=True) as inst:
            offsets_exist = path.exists(inst._byte_offset_filename)
            self.assertEqual(offsets_exist, inst._check_has_byte_offset_file())
            self.assertTrue(isinstance(inst._offset_index, xml.FlatTagSpecificXMLByteIndex))
            self.assertTrue(not isinstance(inst._offset_index, xml.PrebuiltOffsetIndex))
        inst._source.close()
        self.assertTrue(inst._source.closed)
        shutil.rmtree(test_dir, True)

    def test_coerce_duration_type(self):
        with MzXML(self.path) as handle:
            scan = next(handle)
            time = scan['retentionTime']
            self.assertEqual(time.unit_info, 'minute')

    def test_read_dtype(self):
        dtypes = {'m/z array': np.float32, 'intensity array': np.int32}
        with read(self.path, dtype=dtypes) as f:
            for spec in f:
                for k, v in dtypes.items():
                    self.assertEqual(spec[k].dtype, v)

if __name__ == '__main__':
    unittest.main()
