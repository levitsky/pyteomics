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

Main functions
--------------

  :py:func:`charge` - calculate the charge of a polypeptide

  :py:func:`pI` - calculate the isoelectric point of a polypeptide

Data
----

  :py:data:`pK_lehninger` - a set of pK from [#Lehninger]_.

  :py:data:`pK_sillero` - a set of pK from [#Sillero]_.

  :py:data:`pK_dawson` - a set of pK from [#Dawson]_, the pK values for NH2-
  and -OH are taken from [#Sillero]_.

  :py:data:`pK_rodwell` - a set of pK from [#Rodwell]_.

  :py:data:`pK_bjellqvist` - a set of pK from [#Bjellqvist]_.

  :py:data:`pK_nterm_bjellqvist` - a set of N-terminal pK from [#Bjellqvist]_.

  :py:data:`pK_cterm_bjellqvist` - a set of C-terminal pK from [#Bjellqvist]_.

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

.. [#Bjellqvist] Bjellqvist, B., Basse, B., Olsen, E. and Celis, J.E.
    Reference points for comparisons of two-dimensional maps of proteins from
    different human cell types defined in a pH scale where isoelectric points
    correlate with polypeptide compositions. Electrophoresis 1994, 15, 529-539.
    `Link. <http://dx.doi.org/10.1002/elps.1150150171>`_

-------------------------------------------------------------------------------

"""

#   Copyright 2012 Anton Goloborodko, Lev Levitsky
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

from . import parser
from .auxiliary import PyteomicsError
from collections import Iterable

def charge(sequence, pH, **kwargs):
    """Calculate the charge of a polypeptide in given pH or list of pHs using
    a given list of amino acid electrochemical properties.

    .. warning::

        Be cafeful when supplying a list with a parsed sequence or a dict with
        amino acid composition as `sequence`. Such values must be obtained
        with enabled `show_unmodified_termini` option.

    .. warning::

        If you provide `pK_nterm` or `pK_cterm` and provide `sequence` as a dict,
        it is assumed that it was obtained with ``term_aa=True`` (see
        :py:func:`pyteomics.parser.amino_acid_composition` for details).

    Parameters
    ----------
    sequence : str or list or dict
        A string with a polypeptide sequence, a list with a parsed
        sequence or a dict of amino acid composition.
    pH : float or iterable of floats
        pH or iterable of pHs for which the charge is calculated.
    pK : dict {str: [(float, int), ...]}, optional
        A set of pK of amino acids' ionizable groups. It is a dict, where keys
        are amino acid labels and the values are lists of tuples (pK,
        charge_in_ionized_state), a tuple per ionizable group. The default
        value is `pK_lehninger`.

    pK_nterm : dict {str: [(float, int),]}, optional
    pK_cterm : dict {str: [(float, int),]}, optional
        Sets of pK of N-terminal and C-terminal (respectively) amino acids'
        ionizable groups. Dicts with the same structure as ``pK``. These
        values (if present) are used for N-terminal and C-terminal residues,
        respectively. If given, `sequence` must be a :py:class:`str` or a
        :py:class:`list`. The default value is an empty dict.

    Returns
    -------
    out : float or list of floats
        A single value of charge or a list of charges.
    """

    nterm = cterm = n_aa = c_aa = None
    pK = kwargs.get('pK', pK_lehninger).copy()
    pK_nterm = kwargs.get('pK_nterm', {})
    pK_cterm = kwargs.get('pK_cterm', {})

    if isinstance(sequence, dict):
        peptide_dict = sequence.copy()
        for k, v in sequence.items():
            if k[-1] == '-':
                if v > 1 or nterm:
                    raise PyteomicsError(
                            'More that one N-terminal group in {}'.format(
                                sequence))
                nterm = k
            if k[0] == '-':
                if v > 1 or cterm:
                    raise PyteomicsError(
                            'More that one C-terminal group in {}'.format(
                                sequence))
                cterm = k
            if k[:5] == 'nterm':
                if v > 1 or n_aa:
                    raise PyteomicsError(
                            'More that one N-terminal residue in {}'.format(
                                sequence))
                n_aa = k[5:]
                peptide_dict[n_aa] = peptide_dict.get(n_aa, 0) + 1
            if k[:5] == 'cterm':
                if v > 1 or c_aa:
                    raise PyteomicsError(
                            'More that one C-terminal residue in {}'.format(
                                sequence))
                c_aa = k[5:]
                peptide_dict[c_aa] = peptide_dict.get(c_aa, 0) + 1

        if nterm is None or cterm is None:
            raise PyteomicsError('Peptide must have two explicit terminal groups')
        if (n_aa is None or c_aa is None) and (pK_nterm or pK_cterm):
            raise PyteomicsError('Two terminal residues must be present in '
                    'peptide (designated as "ntermX" and "ctermX", where "X" is '
                    'the one-letter residue label). Use '
                    '``term_aa=True`` when calling '
                    '`parser.amino_acid_composition`.')

    elif isinstance(sequence, (str, list)):
        if isinstance(sequence, str):
            parsed_sequence = parser.parse(sequence,
                    show_unmodified_termini=True)
        elif isinstance(sequence, list):
            if sequence[0][-1] != '-' or sequence[-1][0] != '-':
                raise PyteomicsError('Parsed sequences must contain terminal '
                                     'groups at 0-th and last positions.')
            parsed_sequence = sequence

        n_aa = parsed_sequence[1]
        c_aa = parsed_sequence[-2]
        nterm = parsed_sequence[0]
        cterm = parsed_sequence[-1]
        peptide_dict = parser.amino_acid_composition(parsed_sequence)

    else:
        raise PyteomicsError('Unsupported type of sequence: %s' % type(sequence))

    if nterm in pK_nterm:
        if n_aa in pK_nterm[nterm]:
            pK[nterm] = pK_nterm[nterm][n_aa]
    if cterm in pK_cterm:
        if c_aa in pK_cterm[cterm]:
            pK[cterm] = pK_cterm[cterm][c_aa]

    # Process the case when pH is a single float.
    pH_list = pH if isinstance(pH, Iterable) else [pH,]

    # Check if a sequence was parsed with `show_unmodified_termini` enabled.

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


def pI(sequence, pI_range=(0.0, 14.0), precision_pI=0.0001, **kwargs):
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
    pK : dict {str: [(float, int), ...]}, optional
        A set of pK of amino acids' ionizable groups. It is a dict, where keys
        are amino acid labels and the values are lists of tuples (pK,
        charge_in_ionized_state), a tuple per ionizable group. The default
        value is `pK_lehninger`.
    pK_nterm : dict {str: [(float, int),]}, optional
    pK_cterm : dict {str: [(float, int),]}, optional
        Sets of pK of N-terminal and C-terminal (respectively) amino acids'
        ionizable groups. Dicts with the same structure as ``pK``. These
        values (if present) are used for N-terminal and C-terminal residues,
        respectively. If given, `sequence` must be a :py:class:`str` or a
        :py:class:`list`. The default value is an empty dict.

    Returns
    -------
    out : float
    """

    pK = kwargs.get('pK', pK_lehninger.copy())
    pK_nterm = {}
    pK_cterm = {}
    if isinstance(sequence, str) or isinstance(sequence, list):
        pK_nterm = kwargs.get('pK_nterm', {})
        pK_cterm = kwargs.get('pK_cterm', {})
    elif isinstance(sequence,dict) and (("pK_nterm" in kwargs) or ("pK_cterm" in kwargs)):
        raise PyteomicsError('Can not use terminal features for %s' % type(sequence))
    # The algorithm is based on the fact that charge(pH) is a monotonic function.
    left_x, right_x = pI_range
    left_y = charge(sequence, left_x, pK=pK, pK_cterm=pK_cterm, pK_nterm=pK_nterm)
    right_y = charge(sequence, right_x, pK=pK, pK_cterm=pK_cterm, pK_nterm=pK_nterm)
    while (right_x - left_x) > precision_pI:
        if left_y * right_y > 0:
            return left_x if abs(left_y) < abs(right_y) else right_x
        middle_x = (left_x + right_x) / 2.0
        middle_y = charge(sequence, middle_x, pK=pK, pK_cterm=pK_cterm, pK_nterm=pK_nterm)
        if middle_y * left_y < 0:
            right_x = middle_x
            right_y = middle_y
        else:
            left_x = middle_x
            left_y = middle_y
    return (left_x + right_x) / 2.0


pK_lehninger = {
    'E':   [(4.25,  -1)],
    'R':   [(12.48,  1)],
    'Y':   [(10.07, -1)],
    'D':   [(3.65,  -1)],
    'H':   [(6.00,  +1)],
    'K':   [(10.53, +1)],
    'C':   [(8.18,  -1)],
    'H-':  [(9.69,  +1)],
    '-OH': [(2.34,  -1)],
    }
"""A set of pK from Nelson, D. L.; Cox, M. M. Lehninger Principles of
Biochemistry, Fourth Edition; W. H. Freeman, 2004; p. 1100.
"""

pK_sillero = {
    'E':   [(4.5,  -1)],
    'R':   [(12.0, +1)],
    'Y':   [(10.0, -1)],
    'D':   [(4.0,  -1)],
    'H':   [(6.4,  +1)],
    'K':   [(10.4, +1)],
    'C':   [(9.0,  -1)],
    'H-':  [(8.2,  +1)],
    '-OH': [(3.2,  -1)],
    }
"""A set of pK from Sillero, A.; Ribeiro, J. Isoelectric points of proteins:
Theoretical determination. Analytical Biochemistry, vol. 179 (2), pp. 319-325,
1989.
"""

pK_dawson = {
    'E':   [(4.3,  -1)],
    'R':   [(12.0, +1)],
    'Y':   [(10.1, -1)],
    'D':   [(3.9,  -1)],
    'H':   [(6.0,  +1)],
    'K':   [(10.5, +1)],
    'C':   [(8.3,  -1)],
    'H-':  [(8.2,  +1)],
    '-OH': [(3.2,  -1)],
    }
"""A set of pK from Dawson, R. M. C.; Elliot, D. C.; Elliot, W. H.; Jones,
K. M.  Data for biochemical research. Oxford University Press, 1989; p. 592.
pKs for NH2- and -OH are taken from `pK_sillero`.
"""

pK_rodwell = {
    'E':   [(4.25, -1)],
    'R':   [(11.5, +1)],
    'Y':   [(10.7, -1)],
    'D':   [(3.86, -1)],
    'H':   [(6.0,  +1)],
    'K':   [(11.5, +1)],
    'C':   [(8.33, -1)],
    'H-':  [(8.0,  +1)],
    '-OH': [(3.1,  -1)],
}
"""A set of pK from Rodwell, J. Heterogeneity of component bands in
isoelectric focusing patterns. Analytical Biochemistry, vol. 119 (2),
pp. 440-449, 1982.
"""

pK_bjellqvist = {
    'E':   [(4.45, -1)],
    'R':   [(12.0, +1)],
    'Y':   [(10.0, -1)],
    'D':   [(4.05, -1)],
    'H':   [(5.98, +1)],
    'K':   [(10.0, +1)],
    'C':   [(9.0,  -1)],
    'H-':  [(7.5,  +1)],
    '-OH': [(3.55, -1)],
}
"""
A set of pK from Bjellqvist, B., Basse, B., Olsen, E. and Celis, J.E.
Reference points for comparisons of two-dimensional maps of proteins from
different human cell types defined in a pH scale where isoelectric points
correlate with polypeptide compositions. Electrophoresis 1994, 15, 529-539.
"""

pK_nterm_bjellqvist = {
    'H-': {
        'A': [(7.59, +1)],
        'M': [(7.0,  +1)],
        'S': [(6.93, +1)],
        'P': [(8.36, +1)],
        'T': [(6.82, +1)],
        'V': [(7.44, +1)],
        'E': [(7.7,  +1)]
        }
    }
"""
A set of N-terminal pK from Bjellqvist, B., Basse, B., Olsen, E. and Celis, J.E.
Reference points for comparisons of two-dimensional maps of proteins from
different human cell types defined in a pH scale where isoelectric points
correlate with polypeptide compositions. Electrophoresis 1994, 15, 529-539.
"""

pK_cterm_bjellqvist = {
    '-OH': {
        'D': [(4.55, -1)],
        'E': [(4.75, -1)]
        }
    }
"""
A set of C-terminal pK from Bjellqvist, B., Basse, B., Olsen, E. and Celis, J.E.
Reference points for comparisons of two-dimensional maps of proteins from
different human cell types defined in a pH scale where isoelectric points
correlate with polypeptide compositions. Electrophoresis 1994, 15, 529-539.
"""

if __name__ == "__main__":
    import doctest

    doctest.testmod()
