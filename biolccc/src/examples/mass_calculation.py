import pyBioLCCC

peptide = 'Ac-PEPTIDE-NH2'

averageMass = pyBioLCCC.calculateAverageMass(
    peptide, pyBioLCCC.rpAcnFaRod)
monoisotopicMass = pyBioLCCC.calculateMonoisotopicMass(
    peptide, pyBioLCCC.rpAcnFaRod)

print 'The average mass of', peptide, 'is', averageMass, 'Da'
print 'The monoisotopic mass of', peptide, 'is', monoisotopicMass, 'Da'
