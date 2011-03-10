import pyBioLCCC

peptide = 'Ac-PEPTIDE-NH2'
RT = pyBioLCCC.calculateRT(peptide,
    pyBioLCCC.rpAcnFaRod,
    pyBioLCCC.standardChromoConditions)
print 'The retention time of', peptide, 'is', RT
