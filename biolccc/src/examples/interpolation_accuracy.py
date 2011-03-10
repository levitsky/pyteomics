import math
import random
import pylab
import pyBioLCCC

peptides = []
random.seed()
for i in [5] * 5 + [10] * 5 + [20] * 5 + [40] * 5:
    peptides.append(
        ''.join([random.choice("QWERTYIPASDFGHKLCVNM") for j in range(i)]))
pylab.figure(figsize=(9,4))
pylab.subplots_adjust(top=0.85, hspace=0.30, wspace=0.30, right=0.96)
pylab.suptitle('The difference of RT calculated by the fast and standard '
               'algorithm for 20 random peptides')
for chembasis, subplot_num, title in [
(pyBioLCCC.rpAcnTfaChain, 121, 'rpAcnTfaChain'),
(pyBioLCCC.rpAcnFaRod, 122, 'rpAcnFaRod')]:
    pylab.subplot(subplot_num)
    for peptide in peptides:
        x = range(5, 31, 1)
        RTs = [pyBioLCCC.calculateRT(peptide, chembasis,
                 pyBioLCCC.standardChromoConditions, i)
               for i in x]
        ref_RT = pyBioLCCC.calculateRT(peptide, chembasis,
                     pyBioLCCC.standardChromoConditions)
        y = [(RT / ref_RT ) * 100.0 - 100.0 for RT in RTs] 
        pylab.plot(x, y, label=peptide)
    pylab.title(title)
    pylab.xlabel('Number of interpolating points')
    pylab.ylabel('$\Delta RT, \; \%$')
    pylab.rcParams['legend.fontsize'] = 10
    #pylab.legend(loc='upper right')
    #pylab.legend(loc='lower right')
pylab.show()
