import unittest
import platform
import os
import pyteomics
pyteomics.__path__ = [os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir, 'pyteomics'))]
from pyteomics import auxiliary as aux


class UtilTest(unittest.TestCase):
    def test_ensure_prefix(self):
        pairs = [
            ('file:///home/test/unimod.xml', 'file:///home/test/unimod.xml'),
            ('https://example.org/test/unimod.xml', 'https://example.org/test/unimod.xml'),
            ('ftp://example.org/test/unimod.xml', 'ftp://example.org/test/unimod.xml'),
            ('http://example.org/test/unimod.xml', 'http://example.org/test/unimod.xml'),
        ]
        pairs_windows = [
            ('C:/Data folder/unimod.xml', 'file:///C:/Data%20folder/unimod.xml'),
            ('file:///C:/Data folder/unimod.xml', 'file:///C:/Data folder/unimod.xml'),
        ]
        pairs_other = [('/home/test/unimod.xml', 'file:///home/test/unimod.xml'),]
        system = platform.system()
        print('Testing on', system)
        if system == 'Windows':
            pairs.extend(pairs_windows)
        else:
            pairs.extend(pairs_other)
        for inp, out in pairs:
            try:
                self.assertEqual(aux.ensure_url_prefix(inp), out)
            except Exception:
                print('Failed with:', inp, out)
                raise


if __name__ == '__main__':
    unittest.main()
