"""
electrochem - calculate charge and pI of a polypeptide
======================================================

The theory behind this module is based on the Henderson-Hasselbalch
equation and was thoroughly described in a number of sources [#Aronson]_,
[#Moore]_.

Briefly, the formula is the following:

.. math::

   Q_{peptide} = \Sigma{\\frac{Q_i}{1+10^{Q_i(pH-pK_i)}}},

where the sum is taken over all ionizable groups of a polypeptide, and
:math:`Q_i` is -1 and +1 for acidic and basic functional groups,
respectively.

**Computation:**


.. [#Aronson] John N Aronson, The Henderson-Hasselbalch equation
   revisited.  Biochemical Education, vol. 11, issue 2, pp. 68, April
   1983.  `link <http://dx.doi.org/10.1016/0307-4412(83)90046-8>`_

.. [#Moore] Dexter S Moore, Amino acid and peptide net charges: A
   simple calculational procedure. Biochemical Education, vol. 13,
   issue 1, pp. 10-12, January 1986.  `link
   <http://dx.doi.org/10.1016/0307-4412(85)90114-1>`_

-------------------------------------------------------------------------------

"""
# Licensed under the MIT license:
# http://www.opensource.org/licenses/mit-license.php 

import parser

def charge(sequence, pH, **kwargs):
    """Calculate charge of a polypeptide in given pH of list of pHs using
    a given list of amino acid electrochemical properties.

    Keyword arguments:
    sequence -- a polypeptide sequence;
    pH -- a single value of pH or a list of pHs;
    data -- a set of pK of amino acids' charged groups.

    Return a single value of charge or a list of charges.
    """
    data = kwargs.get('data', data_lehniger)
    
    amino_acids = data.keys()
    peptide_dict = parser.amino_acid_composition(sequence, True, False,
                                                 labels=amino_acids)

    # Processing the case what pH is a single float.
    pH_list = pH if isinstance(pH, list) else [pH,]

    # Calculating charge for each value of pH.
    charge_list = []
    for pH_value in pH_list:
        charge = 0
        for aa in peptide_dict:
            for charged_group in data[aa]:
                charge += peptide_dict[aa] * charged_group[1] * (
                    1.0 / (1.0 + 10 ** (
                            charged_group[1] * (pH_value - charged_group[0]))))
        charge_list.append(charge)
    
    return charge_list[0] if len(charge_list) == 1 else charge_list

def pI(sequence, pI_range=(0.0, 14.0), precision_pI=0.01, **kwargs):
    """Calculate isoelectric point of a polypeptide using a given set
    of amino acids' electrochemical properties.

    Parameters
    ----------
    sequence : str
    pI_range: tuple of float
    data : dict, optional
        a set of pK of amino acids' charged groups;

    """
    #     Blahblah
    # precision_pI -- a precision of the calculated pI.

    # Return a value of pI.

    data = kwargs.get('data', data_lehniger)

    left_x, right_x = pI_range
    left_y = charge(sequence, left_x, data)
    right_y = charge(sequence, right_x, data)

    while (right_x-left_x) > precision_pI:
        if left_y * right_y > 0:
            return left_x if abs(left_y) < abs(right_y) else right_x
        
        middle_x = (left_x + right_x) / 2.0
        middle_y = charge(sequence, middle_x, data)
        
        if middle_y * left_y < 0:
            right_x = middle_x
            right_y = middle_y
        else:
            left_x = middle_x
            left_y = middle_y
    return (left_x + right_x) / 2.0

# Nelson, D. L.; Cox, M. M. Lehninger Principles of Biochemistry,
# Fourth Edition; W. H. Freeman, 2004; p. 1100.
data_lehniger = {
    'Q': [],
    'W': [],
    'E': [(4.25, -1),],
    'R': [(12.48, +1),],
    'T': [],
    'Y': [(10.07, -1),],
    'I': [],
    'P': [],
    'A': [],
    'S': [],
    'D': [(3.65, -1),],
    'F': [],
    'G': [],
    'H': [(6.00, +1),],
    'K': [(10.53, +1),],
    'L': [],
    'C': [(8.18, -1),],
    'V': [],
    'N': [],
    'M': [],
    'NH2-': [(9.69, +1),],
    '-OH':  [(2.34, -1),],
    }

# Sillero, A.; Ribeiro, J. Isoelectric points of proteins: Theoretical
# determination. Analytical Biochemistry. 1989, 179, 319-325.
data_sillero_ribeiro = {
    'Q': [],
    'W': [],
    'E': [(4.5, -1),],
    'R': [(12.0, +1),],
    'T': [],
    'Y': [(10.0, -1),],
    'I': [],
    'P': [],
    'A': [],
    'S': [],
    'D': [(4.0, -1),],
    'F': [],
    'G': [],
    'H': [(6.4, +1),],
    'K': [(10.4, +1),],
    'L': [],
    'C': [(9.0, -1),],
    'V': [],
    'N': [],
    'M': [],
    'NH2-': [(8.2, +1),],
    '-OH':  [(3.2, -1),],
    }

# Dawson, R. M. C.; Elliot, D. C.; Elliot, W. H.; Jones, K. M.  
# Data for biochemical research; Oxford University Press, 1989; p. 592.
# pKs for NH2- and -OH are taken from data_sillero_ribeiro
data_dawson = {
    'Q': [],
    'W': [],
    'E': [(4.3, -1),],
    'R': [(12.0, +1),],
    'T': [],
    'Y': [(10.1, -1),],
    'I': [],
    'P': [],
    'A': [],
    'S': [],
    'D': [(3.9, -1),],
    'F': [],
    'G': [],
    'H': [(6.0, +1),],
    'K': [(10.5, +1),],
    'L': [],
    'C': [(8.3, -1),],
    'V': [],
    'N': [],
    'M': [],
    'NH2-': [(8.2, +1),],
    '-OH':  [(3.2, -1),],
    }

# Rodwell, J. Analytical Biochemistry. 1982, 119, 440-449.
data_rodwell = {
    'Q': [],
    'W': [],
    'E': [(4.25, -1),],
    'R': [(11.5, +1),],
    'T': [],
    'Y': [(10.7, -1),],
    'I': [],
    'P': [],
    'A': [],
    'S': [],
    'D': [(3.86, -1),],
    'F': [],
    'G': [],
    'H': [(6.0, +1),],
    'K': [(11.5, +1),],
    'L': [],
    'C': [(8.33, -1),],
    'V': [],
    'N': [],
    'M': [],
    'NH2-': [(8.0, +1),],
    '-OH':  [(3.1, -1),],
    }

