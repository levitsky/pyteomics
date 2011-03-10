import pyBioLCCC

# Deriving a new ChemicalBasis instance from a predefined one.
myChemicalBasis = pyBioLCCC.ChemicalBasis(
    pyBioLCCC.RP_ACN_FA_ROD)

# Changing the bind energy of a chemical group.
myChemicalBasis.chemicalGroups()['E'].setBindEnergy(0.0)
myChemicalBasis.chemicalGroups()['-NH2'].setBindEnergy(0.0)

print "The bind energy of E is", \
    myChemicalBasis.chemicalGroups()['E'].bindEnergy()
print "The bind energy of -NH2 is", \
    myChemicalBasis.chemicalGroups()['-NH2'].bindEnergy()

# Adding a new chemical group. The energy is not valid.
myChemicalBasis.addChemicalGroup(
    pyBioLCCC.ChemicalGroup(
        'Hydroxyproline',      # full name
        'hoP',                 # label
        0.40,                  # bind energy
        97.1167+15.9994,       # average mass
        97.05276+15.9994915))  # monoisotopic mass

# Setting a new type of model. Without a massive recalibration
# it will ruin the accuracy of prediction.
myChemicalBasis.setModel(pyBioLCCC.CHAIN);

peptide = "Ac-PEhoPTIDE-NH2"
RT = pyBioLCCC.calculateRT(peptide,
    myChemicalBasis,
    pyBioLCCC.standardChromoConditions)

monoisotopicMass = pyBioLCCC.calculateMonoisotopicMass(
    peptide, myChemicalBasis)

print 'The retention time of', peptide, 'is', RT
print 'The monoisotopic mass of', peptide, 'is', monoisotopicMass,'Da'
