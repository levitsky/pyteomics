import unittest
import os
import pickle
import pyteomics
from io import BytesIO
pyteomics.__path__ = [os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir, 'pyteomics'))]
from data import mzml_spectra
try:
    from pyteomics.mzmlb import MzMLb, read, chain
    reason = None
except ImportError as err:
    MzMLb = read = chain = None
    reason = err

from pyteomics.auxiliary import FileReader


class MzMLbTest(unittest.TestCase):
    maxDiff = None
    path = 'test.mzMLb'

    def test_read(self):
        for func in [MzMLb, read, chain]:
            with func(self.path) as r:
                # http://stackoverflow.com/q/14246983/1258041
                self.assertEqual(mzml_spectra, list(r))

    def test_picklable(self):
        with MzMLb(self.path) as reader:
            expected_data = next(reader)
            spec = pickle.dumps(reader)
        with pickle.loads(spec) as reader:
            self.assertEqual(next(reader)['id'], expected_data['id'])

    def test_in_memory_buffer(self):
        with open(self.path, 'rb') as fh:
            data_buffer = BytesIO(fh.read())
        with MzMLb(data_buffer) as reader:
            spectrum = next(reader)
            self.assertEqual(
                spectrum['id'], 'controllerType=0 controllerNumber=1 scan=1')
        data_buffer.seek(0)
        with MzMLb(data_buffer, use_index=True) as reader:
            spectrum = next(reader)
            self.assertEqual(
                spectrum['id'], 'controllerType=0 controllerNumber=1 scan=1')

    def test_registered_filereader(self):
        self.assertTrue(issubclass(MzMLb, FileReader))


if __name__ == '__main__':
    unittest.main()
