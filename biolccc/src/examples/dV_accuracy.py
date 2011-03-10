import math
import random
import pylab
import pyBioLCCC

nanocolumn = pyBioLCCC.ChromoConditions()
macrocolumn = pyBioLCCC.ChromoConditions()
macrocolumn.update({
    'columnLength'   : 250.0,
    'columnDiameter' : 4.6,
    'columnPoreSize' : 100.0,
    'flowRate'       : 0.5,
    'dV'             : 0.0,
    'secondSolventConcentrationA' : 0.0,
    'secondSolventConcentrationB' : 100.0,
})

peptides = []
random.seed()
for i in [5] * 3 + [10] * 2 + [20] * 3 + [40] * 2:
    peptides.append(
        ''.join([random.choice("QWERTYIPASDFGHKLCVNM") for j in range(i)]))
pylab.figure(figsize=(9,8))
pylab.subplots_adjust(top=0.9, hspace=0.30, wspace=0.30, right=0.96)
pylab.suptitle('The effect of changing dV for ten random peptides'
               ' of different length')
for chembasis, chromoconditions, subplot_num, title in [
(pyBioLCCC.rpAcnTfaChain, nanocolumn, 221, 'Nanocolumn, rpAcnTfaChain'),
(pyBioLCCC.rpAcnFaRod, nanocolumn, 222, 'Nanocolumn, rpAcnFaRod'),
(pyBioLCCC.rpAcnTfaChain, macrocolumn, 223, 'Macrocolumn, rpAcnTfaChain'),
(pyBioLCCC.rpAcnFaRod, macrocolumn, 224, 'Macrocolumn, rpAcnFaRod')]:
    pylab.subplot(subplot_num)
    for peptide in peptides:
        divisors = [(2.0 ** i) for i in range(11)]
        chromoconditions['dV'] = chromoconditions['flowRate'] / divisors[-1]
        reference_time = ( 
                pyBioLCCC.calculateRT(peptide, chembasis, chromoconditions))
        y = []
        for divisor in divisors:
            chromoconditions['dV'] = chromoconditions['flowRate'] / divisor
            y.append(
                (pyBioLCCC.calculateRT(peptide, chembasis, chromoconditions) /
                reference_time - 1.0) * 100.0)
        pylab.plot(divisors, y, label=peptide)
    pylab.title(title)
    pylab.xlabel('Flow rate divisor')
    pylab.ylabel('$\Delta RT, \; \%$')
    pylab.gca().set_xscale('log')
    pylab.rcParams['legend.fontsize'] = 10
    #pylab.legend(loc='upper right')
    #pylab.legend(loc='lower right')
pylab.show()
