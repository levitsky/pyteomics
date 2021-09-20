from pyteomics import pylab_aux as pa, usi
import matplotlib.pyplot as plt
spectrum = usi.proxi(
    'mzspec:PXD004732:01650b_BC2-TUM_first_pool_53_01_01-3xHCD-1h-R2:scan:41840',
    'massive')
peptide = 'WNQLQAFWGTGK'

pa.annotate_spectrum(spectrum, peptide, precursor_charge=2, backend='spectrum_utils',
    ion_types='aby', title=peptide)

from pyteomics import mass
peptide = 'DLTDYLoxMK'  # oxidized methionine
aa_mass = mass.std_aa_mass.copy()
aa_mass['ox'] = 15.9949  # define the mass of the label

usi_top = 'mzspec:MSV000079960:DY_HS_Exp7-Ad1:scan:30372'
usi_bottom = 'mzspec:MSV000080679:j11962_C1orf144:scan:10671'

spectrum_top = usi.proxi(usi_top, 'massive')
spectrum_bottom = usi.proxi(usi_bottom, 'massive')

fig, ax = plt.subplots(figsize=(12, 6))
pa.mirror(spectrum_top, spectrum_bottom, peptide=peptide, precursor_charge=2,
    aa_mass=aa_mass, ion_types='aby', ax=ax, title=peptide, ftol=0.5, scaling='root',
    remove_precursor_peak=True)

plt.show()
