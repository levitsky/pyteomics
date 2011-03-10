import math
import random
import pylab
import pyBioLCCC

peptides = []
random.seed()
for i in [5, 10, 20, 30, 40]:
    peptides.append(
        ''.join([random.choice("QWERTYIPASDFGHKLCVNM") for j in range(i)]))
pylab.figure(figsize=(9,4))
pylab.subplots_adjust(top=0.8, hspace=0.30, wspace=0.30, right=0.96)
pylab.suptitle('The dependency of log(Kd) on the second solvent concentration '
               'for five random peptides')
for chembasis, subplot_num, title in [
(pyBioLCCC.rpAcnTfaChain, 121, 'rpAcnTfaChain'),
(pyBioLCCC.rpAcnFaRod, 122, 'rpAcnFaRod')]:
    pylab.subplot(subplot_num)
    for peptide in peptides:
        x = range(0, 101, 1)
        y = [math.log(pyBioLCCC.calculateKd(peptide, i, chembasis))
             for i in x] 
        pylab.plot(x, y, label='%d aa residues' % (len(peptide)))
    pylab.rcParams['legend.fontsize'] = 10
    pylab.legend(loc='upper right')
    pylab.xlim((0,100))
    pylab.xlabel('ACN concentration, %')
    pylab.ylabel('log(Kd)')
    pylab.title(title)
pylab.show()
