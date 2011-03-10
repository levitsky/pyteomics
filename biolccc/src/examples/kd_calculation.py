import pyBioLCCC

peptide = 'Ac-PEPTIDE-NH2'

kd = pyBioLCCC.calculateKd(
    peptide, # the peptide sequence 
    15.0,    # the concentration of the second solvent, %
    pyBioLCCC.rpAcnFaRod, # the chemical basis
    100.0)   # the size of the pores, angstroms

print 'The coefficient of distribution of', peptide, 'is', kd
