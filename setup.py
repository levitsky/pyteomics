#!/usr/bin/env python
import ez_setup
ez_setup.use_setuptools()

'''
setup.py file for pyteomics
'''

from setuptools import setup
import re

with open('VERSION') as v:
    version = next(v).strip()
with open('README.rst') as r, open('INSTALL') as i:
    long_description = re.sub(r':py:\w+:`([^`]+)`',
            lambda m: '**{}**'.format(m.group(1)),
            ''.join(r) + '\n' + ''.join(i))

setup(
    name               = 'pyteomics',
    version            = version,
    description        = 'A framework for proteomics data analysis.',
    long_description   = long_description,
    author             = 'Anton Goloborodko & Lev Levitsky',
    author_email       = 'pyteomics@googlegroups.com',
    url                = 'http://pyteomics.readthedocs.io',
    packages           = ['pyteomics', 'pyteomics.mass', 'pyteomics.openms', 'pyteomics.auxiliary'],
    project_urls       = {
                          'Documentation': 'http://pyteomics.readthedocs.io',
                          'Source Code'  : 'https://bitbucket.org/levitsky/pyteomics',
                          'Issue Tracker': 'http://bitbucket.org/levitsky/pyteomics/issues',
                          'Mailing List' : 'https://groups.google.com/group/pyteomics',
                         },
    namespace_packages = ['pyteomics'],
    extras_require     = {'XML': ['lxml', 'numpy'],
                          'TDA': ['numpy'],
                          'graphics': ['matplotlib'],
                          'DF': ['pandas'],
                          'Unimod': ['lxml', 'sqlalchemy'],
                          'numpress': ['pynumpress']},
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
