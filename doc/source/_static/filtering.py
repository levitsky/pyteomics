# place this file in the same directory as example files
import pylab
from pyteomics import tandem, pepxml, mzid, auxiliary as aux, pylab_aux as pa
import pandas as pd
import numpy as np

pylab.figure()
with tandem.read('example.t.xml') as tf:
    h = pylab.hist([psm['rt'] / 60 for psm in tf], bins=25, label='All IDs')
with tandem.filter('example.t.xml', fdr=0.01, full_output=False) as ftf:
    pylab.hist([psm['rt'] / 60 for psm in ftf], bins=h[1], label='1% FDR')
pylab.xlabel('RT, min')
pylab.legend()

q1 = pepxml.qvalues('example.pep.xml', read_schema=False,
                    key=lambda x: x['search_hit'][0]['search_score']['Morpheus Score'], reverse=True)
q2 = tandem.qvalues('example.t.xml')
pylab.figure()
pa.plot_qvalue_curve(q1['q'], label='Morpheus')
pa.plot_qvalue_curve(q2['q'], label='X!Tandem')
pylab.legend()

msgf = mzid.filter('example.mzid', retrieve_refs=True,
                   key=lambda x: x['SpectrumIdentificationItem'][0]['MS-GF:EValue'], fdr=0.01)
pylab.figure()
pylab.hist([psm['SpectrumIdentificationItem'][0]['chargeState'] for psm in msgf], bins=np.arange(5), align='left')
pylab.xticks(np.linspace(0, 4, 5))
pylab.xlabel('charge state')

morpheus = pd.read_table('example.PSMs.tsv')
amanda = pd.read_table('example_output.csv', skiprows=1)

morph_filt = aux.filter(morpheus, fdr=0.01, key='Morpheus Score', reverse=True,
                       is_decoy='Decoy?')

morph_filt.plot(x='Retention Time (minutes)' , y='Precursor Mass (Da)', kind='scatter')

amanda['isDecoy'] = [all(s.startswith('DECOY') for s in prot.split(';')) for prot in amanda['Protein Accessions']]
amanda_filt = aux.filter(amanda[amanda['Rank'] == 1], key='Weighted Probability', is_decoy='isDecoy', fdr=0.01)

amanda_pep = amanda_filt.sort_values('Weighted Probability').groupby('Sequence').first()
morph_pep = morph_filt.sort_values('Q-Value (%)').groupby('Base Peptide Sequence').first()

inter = amanda_pep.join(morph_pep, how='inner', lsuffix='[amanda]', rsuffix='[morpheus]')
inter.plot('Amanda Score', 'Morpheus Score', kind='hexbin', gridsize=10)
pylab.show()
