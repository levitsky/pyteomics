import pylab
import numpy as np
from pyteomics import electrochem

pHs = np.arange(1, 14, 0.5) # list of values of pH
charges = electrochem.charge('PEPTIDE', pHs) # charge function accepts lists of pHs

pylab.figure()
pylab.plot(pHs, charges)
pylab.title("Charge of peptide 'PEPTIDE' vs pH")
pylab.xlabel('pH')
pylab.ylabel('Charge')
pylab.show()
