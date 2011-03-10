import pyBioLCCC
peptide = 'PEPTIDE'

parsedSequence = pyBioLCCC.parseSequence(peptide)

for chemicalGroup in parsedSequence:
    print chemicalGroup.name()
