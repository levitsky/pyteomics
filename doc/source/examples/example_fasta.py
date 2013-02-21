import os
from urllib import urlretrieve
# In python 3, use 'from urllib.request import urlretrieve' instead
import gzip
import matplotlib.pyplot as plt
import numpy as np
from pyteomics import fasta, parser, mass, achrom, electrochem, auxiliary

if not os.path.isfile('yeast.fasta.gz'):
    print 'Downloading the FASTA file for Saccharomyces cerevisiae...'
    urlretrieve(
        'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/'
        'knowledgebase/proteomes/YEAST.fasta.gz',
        'yeast.fasta.gz')
    print 'Done!'

print 'Cleaving the proteins with trypsin...'
unique_peptides = set()
for description, sequence in fasta.read(gzip.open('yeast.fasta.gz')):
    new_peptides = parser.cleave(sequence, parser.expasy_rules['trypsin'])
    unique_peptides.update(new_peptides)
print 'Done, {0} sequences obtained!'.format(len(unique_peptides))

peptides = [{'sequence': i} for i in unique_peptides]

print 'Parsing peptide sequences...'
for peptide in peptides:
    peptide['parsed_sequence'] = parser.parse(
        peptide['sequence'],
        show_unmodified_termini=True)
    peptide['length'] = parser.length(peptide['parsed_sequence']) 
print 'Done!'

peptides = [peptide for peptide in peptides if peptide['length'] <= 100]

print 'Calculating the mass, charge and m/z...'
for peptide in peptides:
    peptide['charge'] = int(round(
        electrochem.charge(peptide['parsed_sequence'], pH=2.0)))
    peptide['mass'] = mass.calculate_mass(peptide['parsed_sequence'])
    peptide['m/z'] = mass.calculate_mass(peptide['parsed_sequence'], 
        charge=peptide['charge'])
print 'Done!'

print 'Calculating the retention time...'
for peptide in peptides:
    peptide['RT_RP'] = achrom.calculate_RT(
        peptide['parsed_sequence'],
        achrom.RCs_zubarev)
    peptide['RT_normal'] = achrom.calculate_RT(
        peptide['parsed_sequence'],
        achrom.RCs_yoshida_lc)
print 'Done!'

plt.figure()
plt.hist([peptide['m/z'] for peptide in peptides], 
    bins = 2000,
    range=(0,4000))
plt.xlabel('m/z, Th')
plt.ylabel('# of peptides within 2 Th bin')

plt.figure()
plt.hist([peptide['charge'] for peptide in peptides],
    bins = 20,
    range=(0,10))
plt.xlabel('charge, e')
plt.ylabel('# of peptides')

x = [peptide['RT_RP'] for peptide in peptides]
y = [peptide['RT_normal'] for peptide in peptides]
heatmap, xbins, ybins = np.histogram2d(x, y, bins=100)
heatmap[heatmap == 0] = np.nan
a, b, r, stderr = auxiliary.linear_regression(x,y)

plt.figure()
plt.imshow(heatmap)
plt.xlabel('RT on RP, min')
plt.ylabel('RT on normal phase, min')
plt.title('All tryptic peptides, RT correlation = {0}'.format(r))

x = [peptide['m/z'] for peptide in peptides]
y = [peptide['RT_RP'] for peptide in peptides]
heatmap, xbins, ybins = np.histogram2d(x, y, 
    bins=[150, 2000],
    range=[[0, 4000], [0, 150]])
heatmap[heatmap == 0] = np.nan
a, b, r, stderr = auxiliary.linear_regression(x,y)

plt.figure()
plt.imshow(heatmap, 
    aspect='auto',
    origin='lower')
plt.xlabel('m/z, Th')
plt.ylabel('RT on RP, min')
plt.title('All tryptic peptides, correlation = {0}'.format(r))

close_mass_peptides = [peptide for peptide in peptides 
                       if 700.0 <= peptide['m/z'] <= 701.0]
x = [peptide['RT_RP'] for peptide in close_mass_peptides]
y = [peptide['RT_normal'] for peptide in close_mass_peptides]
a, b, r, stderr = auxiliary.linear_regression(x, y)

plt.figure()
plt.scatter(x, y)
plt.xlabel('RT on RP, min')
plt.ylabel('RT on normal phase, min')
plt.title('Tryptic peptides with m/z=700-701 Th\nRT correlation = {0}'.format(r))

plt.show()
