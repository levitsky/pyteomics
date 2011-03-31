#!/usr/bin/env python

'''
setup.py file for pyteomics
'''

from distutils.core import setup
import os.path
import glob

setup(
    name = 'pyteomics',
    version = version,
    description      = '''A framework for proteomics data analysis.''',
    long_description = ''.join(open('README').readlines()),
    author           = 'Anton Goloborodko',
    author_email     = 'goloborodko.anton@gmail.com',
    url              = 'http://pyteomics.theorchromo.ru',
    classifiers      = ['Intended Audience :: Science/Research',
                        'Programming Language :: Python :: 2',
                        'Topic :: Education',
                        'Topic :: Scientific/Engineering :: Bio-Informatics',
                        'Topic :: Scientific/Engineering :: Chemistry',
                        'Topic :: Scientific/Engineering :: Physics',
                        'Topic :: Software Development :: Libraries'],
    license          = 'License :: OSI Approved :: MIT License',
    )
