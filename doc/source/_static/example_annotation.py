from pyteomics import pylab_aux as pa, usi
spectrum = usi.proxi('mzspec:PXD004732:01650b_BC2-TUM_first_pool_53_01_01-3xHCD-1h-R2:scan:41840', 'massive')
peptide = 'WNQLQAFWGTGK'

pa.annotate_spectrum(spectrum, peptide, precursor_charge=2, backend='spectrum_utils',
     scaling='root', min_intensity=0.05, ion_types='aby', remove_precursor_peak=True, colors={'a': 'yellow'})
