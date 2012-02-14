#!/usr/bin/env python

'''
setup.py file for pyteomics
'''

from distutils.core import setup

version = open('VERSION').readline().strip()

setup(
    name = 'pyteomics',
    version = version,
    description      = '''A framework for proteomics data analysis.''',
    long_description = (''.join(open('README').readlines()) + '\n'
                        + ''.join(open('INSTALL').readlines())),
    author           = 'Anton Goloborodko & Lev Levitskiy',
    author_email     = 'goloborodko.anton@gmail.com',
    url              = 'http://hg.theorchromo.ru/pyteomics',
    packages         = ['pyteomics', ],
    classifiers      = ['Intended Audience :: Science/Research',
                        'Programming Language :: Python :: 2.7',
                        'Topic :: Education',
                        'Topic :: Scientific/Engineering :: Bio-Informatics',
                        'Topic :: Scientific/Engineering :: Chemistry',
                        'Topic :: Scientific/Engineering :: Physics',
                        'Topic :: Software Development :: Libraries'],
    license          = 'License :: OSI Approved :: MIT License',
    )
