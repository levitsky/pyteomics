from os import path
import doctest
import unittest
import pyteomics
pyteomics.__path__ = [path.abspath(path.join(path.dirname(__file__), path.pardir, 'pyteomics'))]
from pyteomics import achrom
# doctest.testmod(achrom, verbose=True)


def load_tests(loader, tests, ignore):
    tests.addTests(doctest.DocTestSuite(achrom))
    return tests


if __name__ == "__main__":
    unittest.main()
