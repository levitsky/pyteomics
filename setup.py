#!/usr/bin/env python

'''
setup.py file for pyteomics
'''

from setuptools import setup
import re
import os


# from https://packaging.python.org/guides/single-sourcing-package-version/

def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with open(os.path.join(here, rel_path), 'r') as fp:
        return fp.read()


def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith('__version__'):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")


with open('README.rst') as r, open('INSTALL') as i:
    long_description = re.sub(r':py:\w+:`([^`]+)`',
            lambda m: '**{}**'.format(m.group(1)),
            ''.join(r) + '\n' + ''.join(i))


extras_require = {'XML': ['lxml', 'numpy'],
                  'TDA': ['numpy'],
                  'graphics': ['matplotlib'],
                  'DF': ['pandas>=0.17'],
                  'Unimod': ['lxml', 'sqlalchemy'],
                  'numpress': ['pynumpress'],
                  'mzMLb': ['h5py', 'hdf5plugin'],
                  'proforma': ['psims > v0.1.42']}
extras_require['all'] = sum(extras_require.values(), ['scikit-learn'])


setup(
    name               = 'pyteomics',
    version            = get_version('pyteomics/version.py'),
    description        = 'A framework for proteomics data analysis.',
    long_description   = long_description,
    author             = 'Anton Goloborodko & Lev Levitsky',
    author_email       = 'pyteomics@googlegroups.com',
    url                = 'http://pyteomics.readthedocs.io',
    packages           = ['pyteomics', 'pyteomics.mass', 'pyteomics.openms', 'pyteomics.auxiliary'],
    project_urls       = {
                          'Documentation': 'http://pyteomics.readthedocs.io',
                          'Source Code'  : 'https://github.com/levitsky/pyteomics',
                          'Issue Tracker': 'https://github.com/levitsky/pyteomics/issues',
                          'Mailing List' : 'https://groups.google.com/group/pyteomics',
                         },
    namespace_packages = ['pyteomics'],
    extras_require     = extras_require,
    classifiers        = ['Intended Audience :: Science/Research',
                          'Programming Language :: Python :: 2.7',
                          'Programming Language :: Python :: 3',
                          'Topic :: Education',
                          'Topic :: Scientific/Engineering :: Bio-Informatics',
                          'Topic :: Scientific/Engineering :: Chemistry',
                          'Topic :: Scientific/Engineering :: Physics',
                          'Topic :: Software Development :: Libraries'],
    license            = 'License :: OSI Approved :: Apache Software License',
    zip_safe           = False,
    )
