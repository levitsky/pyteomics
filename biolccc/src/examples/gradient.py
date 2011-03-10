import pyBioLCCC

peptide = 'Ac-PEPTIDE-NH2'
myChromoConditions = pyBioLCCC.ChromoConditions()
myGradient = pyBioLCCC.Gradient()
myGradient.addPoint(0.0, 5.0)
myGradient.addPoint(20.0, 5.0)
myGradient.addPoint(60.0, 45.0)
myGradient.addPoint(65.0, 100.0)
myGradient.addPoint(85.0, 100.0)
myChromoConditions.setGradient(myGradient)

RT = pyBioLCCC.calculateRT(peptide,
         pyBioLCCC.rpAcnFaRod,
         myChromoConditions)
print 'The retention time of', peptide, 'in the custom gradient is',RT
