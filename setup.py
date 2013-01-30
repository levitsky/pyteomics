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
    author           = 'Anton Goloborodko & Lev Levitsky',
    author_email     = 'pyteomics@googlegroups.com',
    url              = 'http://hg.theorchromo.ru/pyteomics',
    packages         = ['pyteomics', ],
    requires         = [line.strip() for line in open('requirements.txt')],
    classifiers      = ['Intended Audience :: Science/Research',
                        'Programming Language :: Python :: 2.7',
                        'Programming Language :: Python :: 3.2',
                        'Topic :: Education',
                        'Topic :: Scientific/Engineering :: Bio-Informatics',
                        'Topic :: Scientific/Engineering :: Chemistry',
                        'Topic :: Scientific/Engineering :: Physics',
                        'Topic :: Software Development :: Libraries'],
    license          = 'License :: OSI Approved :: Apache Software License',
    )
