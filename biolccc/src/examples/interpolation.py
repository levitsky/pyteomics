import pyBioLCCC

print pyBioLCCC.calculateRT('QWERTYIPASDFGHKLCVNM', pyBioLCCC.rpAcnTfaChain,
    pyBioLCCC.standardChromoConditions)

# Using 21 interpolating points.
print pyBioLCCC.calculateRT('QWERTYIPASDFGHKLCVNM', pyBioLCCC.rpAcnTfaChain,
    pyBioLCCC.standardChromoConditions, 21)
