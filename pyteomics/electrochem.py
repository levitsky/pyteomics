"""
electrochem - electrochemical properties of polypeptides
========================================================

Summary
-------

This module is used to calculate the
electrochemical properties of polypeptide molecules.

The theory behind this module is based on the Henderson-Hasselbalch
equation and was thoroughly described in a number of sources [#Aronson]_,
[#Moore]_.

Briefly, the formula for the charge of a polypeptide in given pH is the following:

.. math::

   Q_{peptide} = \sum{\\frac{Q_i}{1+10^{Q_i(pH-pK_i)}}},

where the sum is taken over all ionizable groups of the polypeptide, and
:math:`Q_i` is -1 and +1 for acidic and basic functional groups,
respectively.

Main functions:
---------------

  :py:func:`charge` - calculate the charge of a polypeptide
  
  :py:func:`pI` - calculate the isoelectric point of a polypeptide

Data:
-----

  :py:data:`pK_lehninger` - a set of pK from [#Lehninger]_.

  :py:data:`pK_sillero` - a set of pK from [#Sillero]_.

  :py:data:`pK_dawson` - a set of pK from [#Dawson]_, the pK values for NH2-
  and -OH are taken from [#Sillero]_.

  :py:data:`pK_rodwell` - a set of pK from [#Rodwell]_.

References
----------

.. [#Aronson] Aronson, J. N. The Henderson-Hasselbalch equation
   revisited.  Biochemical Education, 1983, 11 (2), 68.
   `Link. <http://dx.doi.org/10.1016/0307-4412(83)90046-8>`_

.. [#Moore] Moore, D. S.. Amino acid and peptide net charges: A
   simple calculational procedure. Biochemical Education, 1986, 13 (1), 10-12.
   `Link. <http://dx.doi.org/10.1016/0307-4412(85)90114-1>`_

.. [#Lehninger] Nelson, D. L.; Cox, M. M. Lehninger Principles of
   Biochemistry, Fourth Edition; W. H. Freeman, 2004; p. 1100.

.. [#Sillero] Sillero, A.; Ribeiro, J. Isoelectric points of proteins:
   Theoretical determination. Analytical Biochemistry, 1989, 179 (2), 319-325.
   `Link. <http://dx.doi.org/10.1016/0003-2697(89)90136-X>`_

.. [#Dawson] Dawson, R. M. C.; Elliot, D. C.; Elliot, W. H.; Jones, K. M.
   Data for biochemical research. Oxford University Press, 1989; p. 592.

.. [#Rodwell] Rodwell, J. Heterogeneity of component bands in isoelectric
   focusing patterns. Analytical Biochemistry, 1982, 119 (2), 440-449.
   `Link. <http://dx.doi.org/10.1016/0003-2697(82)90611-X>`_

.. ipython::
   :suppress:

   In [1]: import pyteomics.parser; from pprint import pprint

-------------------------------------------------------------------------------

"""
# Licensed under the MIT license:
# http://www.opensource.org/licenses/mit-license.php 

import parser
from auxiliary import PyteomicsError

def charge(sequence, pH, **kwargs):
    """Calculate the charge of a polypeptide in given pH or list of pHs using
    a given list of amino acid electrochemical properties.

    .. warning::

        Be cafeful when supplying a list with a parsed sequence or a dict with
        amino acid composition as `sequence`. Such values must be obtained
        with enabled `show_unmodified_termini` option.

    Parameters
    ----------
    sequence : str or list or dict    
        A string with a polypeptide sequence, a list with a parsed
        sequence or a dict of amino acid composition.
    pH : float or list of floats
        pH or list of pHs for which the charge is calculated.
    pK : dict {str: [(float, int),]}, optional
        A set of pK of amino acids' ionizable groups. It is a dict, where keys
        are amino acid labels and the values are lists of tuples (pK,
        charge_in_ionized_state), a tuple per ionizable group. The default
        value is `pK_lehninger`.

    Returns
    -------
    out : float or list of floats or None    
        A single value of charge or a list of charges. Returns None if
        `sequence` is not of supported type.
    """

    # Get the list of valid modX labels.
    pK = kwargs.get('pK', pK_lehninger)
    labels = list(parser.std_labels)
    for label in pK:
        if label not in labels:
            labels.append(label)

    # Parse the sequence.
    if isinstance(sequence, basestring) or isinstance(sequence, list):
        peptide_dict = parser.amino_acid_composition(sequence, True, False,
                                                     labels=labels)
    elif isinstance(sequence, dict):
        peptide_dict = sequence
    else:
        raise PyteomicsError('Unsupported type of a sequence.')

    # Check if a sequence was parsed with `show_unmodified_termini` enabled.
    num_term_mod = 0
    for aa in peptide_dict:
        if parser.is_term_mod(aa):
            num_term_mod += 1
    if num_term_mod != 2:
        raise PyteomicsError('Parsed sequences must contain unmodified termini.')
        
    # Process the case when pH is a single float.
    pH_list = pH if isinstance(pH, list) else [pH,]

    # Calculate the charge for each value of pH.
    charge_list = []
    for pH_value in pH_list:
        charge = 0
        for aa in peptide_dict:
            for ionizable_group in pK.get(aa, []):
                charge += peptide_dict[aa] * ionizable_group[1] * (
                    1.0 
                    / (1.0 + 10 ** (ionizable_group[1]
                                  * (pH_value - ionizable_group[0]))))
        charge_list.append(charge)
    
    return charge_list[0] if len(charge_list) == 1 else charge_list

def pI(sequence, pI_range=(0.0, 14.0), precision_pI=0.01, **kwargs):
    """Calculate the isoelectric point of a polypeptide using a given set
    of amino acids' electrochemical properties.

    .. warning::

        Be cafeful when supplying a list with a parsed sequence or a dict with
        amino acid composition as `sequence`. Such values must be obtained
        with enabled `show_unmodified_termini` option.

    Parameters
    ----------
    sequence : str or list or dict    
        A string with a polypeptide sequence, a list with a parsed
        sequence or a dict of amino acid composition.
    pI_range : tuple (float, float)
        The range of allowable pI values. Default is (0.0, 14.0).
    precision_pI : float
        The precision of the calculated pI. Default is 0.01.
    pK : dict {str: [(float, int),]}, optional
        A set of pK of amino acids' ionizable groups. It is a dict, where keys
        are amino acid labels and the values are lists of tuples (pK,
        charge_in_ionized_state), a tuple per ionizable group. The default
        value is `pK_lehninger`.
        
    Returns
    -------
    out : float
    """

    pK = kwargs.get('pK', pK_lehninger)

    # The algorithm is based on the fact that charge(pH) is a monotonic function.
    left_x, right_x = pI_range
    left_y = charge(sequence, left_x, pK=pK)
    right_y = charge(sequence, right_x, pK=pK)

    while (right_x - left_x) > precision_pI:
        if left_y * right_y > 0:
            return left_x if abs(left_y) < abs(right_y) else right_x
        
        middle_x = (left_x + right_x) / 2.0
        middle_y = charge(sequence, middle_x, pK=pK)
        
        if middle_y * left_y < 0:
            right_x = middle_x
            right_y = middle_y
        else:
            left_x = middle_x
            left_y = middle_y
    return (left_x + right_x) / 2.0

pK_lehninger = {
    'E': [(4.25, -1),],
    'R': [(12.48, +1),],
    'Y': [(10.07, -1),],
    'D': [(3.65, -1),],
    'H': [(6.00, +1),],
    'K': [(10.53, +1),],
    'C': [(8.18, -1),],
    'H-': [(9.69, +1),],
    '-OH':  [(2.34, -1),],
    }
"""A set of pK from Nelson, D. L.; Cox, M. M. Lehninger Principles of
Biochemistry, Fourth Edition; W. H. Freeman, 2004; p. 1100.

.. ipython::
   
   In [2]: pprint(pyteomics.electrochem.pK_lehninger)
"""

pK_sillero = {
    'E': [(4.5, -1),],
    'R': [(12.0, +1),],
    'Y': [(10.0, -1),],
    'D': [(4.0, -1),],
    'H': [(6.4, +1),],
    'K': [(10.4, +1),],
    'C': [(9.0, -1),],
    'H-': [(8.2, +1),],
    '-OH':  [(3.2, -1),],
    }
"""A set of pK from Sillero, A.; Ribeiro, J. Isoelectric points of proteins:
Theoretical determination. Analytical Biochemistry, vol. 179 (2), pp. 319-325,
1989.
   
.. ipython::
   
   In [2]: pprint(pyteomics.electrochem.pK_sillero)
"""

pK_dawson = {
    'E': [(4.3, -1),],
    'R': [(12.0, +1),],
    'Y': [(10.1, -1),],
    'D': [(3.9, -1),],
    'H': [(6.0, +1),],
    'K': [(10.5, +1),],
    'C': [(8.3, -1),],
    'H-': [(8.2, +1),],
    '-OH':  [(3.2, -1),],
    }
"""A set of pK from Dawson, R. M. C.; Elliot, D. C.; Elliot, W. H.; Jones,
K. M.  Data for biochemical research. Oxford University Press, 1989; p. 592.
pKs for NH2- and -OH are taken from `pK_sillero`.

.. ipython::
   
   In [2]: pprint(pyteomics.electrochem.pK_dawson)
"""

pK_rodwell = {
    'E': [(4.25, -1),],
    'R': [(11.5, +1),],
    'Y': [(10.7, -1),],
    'D': [(3.86, -1),],
    'H': [(6.0, +1),],
    'K': [(11.5, +1),],
    'C': [(8.33, -1),],
    'H-': [(8.0, +1),],
    '-OH':  [(3.1, -1),],
    }
"""A set of pK from Rodwell, J. Heterogeneity of component bands in
isoelectric focusing patterns. Analytical Biochemistry, vol. 119 (2),
pp. 440-449, 1982.

.. ipython::
   
   In [2]: pprint(pyteomics.electrochem.pK_rodwell)
"""

if __name__ == "__main__":
    import doctest
    doctest.testmod()
