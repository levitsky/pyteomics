#!/usr/bin/env python

"""
setup.py file for pyteomics
"""

from distutils.core import setup, Extension
import os.path
import glob

sources = (glob.glob(os.path.join('biolccc', 'src', 'core', '*.cpp'))
           + [os.path.join('biolccc', 'src', 'bindings', 'biolccc_wrap.cc')])
version = open('VERSION').readline().strip()

# Compile a shared library for biolccc module.
pyBioLCCC_ext = Extension(
    '_biolccc',
    sources=sources,
    include_dirs=[os.path.join('biolccc','include')],
    extra_compile_args=['-DVERSION=\"%s\"' % version,]
    )

setup(
    name = 'pyteomics',
    version = version,
    description      = """Collection of tools for proteomics data analysis""",
    long_description = ''.join(open('README').readlines()),
    author           = 'Anton Goloborodko',
    author_email     = 'goloborodko.anton@gmail.com',
    url              = 'http://theorchromo.ru',
    ext_package      = 'pyteomics.biolccc', # Store _biolccc.so in biolccc dir.
    ext_modules      = [pyBioLCCC_ext],     # Add _biolccc.so to package
    classifiers      = ['Intended Audience :: Science/Research',
                        'Topic :: Scientific/Engineering :: Bio-Informatics',
                        'Topic :: Scientific/Engineering :: Chemistry',
                        'Topic :: Scientific/Engineering :: Physics',],
    license          = 'License :: Free for non-commercial use',
    packages         = ['pyteomics', 'pyteomics.biolccc'],
    package_dir      = {'pyteomics': '',
                        'pyteomics.biolccc': './biolccc/src/bindings'},
    )
