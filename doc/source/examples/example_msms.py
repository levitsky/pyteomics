# This is written in Python 2 for simplicity
# Can be done forward-compatible easily, though
from pyteomics import mgf, pepxml, mass
import os
from urllib import urlretrieve
import pylab

# get the files
for fname in ('mgf', 'pep.xml'):
    if not os.path.isfile('example.' + fname):
        urlretrieve('http://packages.python.org/pyteomics/_downloads/example.'
                + fname, 'example.' + fname)

def fragments(peptide, types=('b', 'y'), maxcharge=1):
    """
    The function generates all possible m/z for fragments of types 
    `types` and of charges from 1 to `maxharge`.
    """
    for i in xrange(1, len(peptide)-1):
        for ion_type in types:
            for charge in xrange(1, maxcharge+1):
                if ion_type[0] in 'abc':
                    yield mass.fast_mass(
                            peptide[:i], ion_type=ion_type, charge=charge)
                else:
                    yield mass.fast_mass(
                            peptide[i:], ion_type=ion_type, charge=charge)

with mgf.read('example.mgf') as spectra, pepxml.read('example.pep.xml') as psms:
    spectrum = next(spectra)
    psm = next(psms)
pylab.figure()
pylab.title('Theoretical and experimental spectra for '
        + psm['search_hit'][0]['peptide'])
pylab.xlabel('m/z, Th')
pylab.ylabel('Intensity, rel. units')
pylab.bar(spectrum['m/z array'], spectrum['intensity array'], width=0.1, linewidth=2,
        edgecolor='black')
theor_spectrum = list(fragments(psm['search_hit'][0]['peptide'],
    maxcharge=psm['assumed_charge']))
pylab.bar(theor_spectrum,
        [spectrum['intensity array'].max()]*len(theor_spectrum),
        width=0.1, edgecolor='red', alpha=0.7)
pylab.show()
