"""
achrom - additive model of polypeptide chromatography
=====================================================

Summary
-------

The additive model of polypeptide chromatography or the model of
retention coefficients was the earliest attempt to describe the dependence of
retention time of a polypeptide in liquid chromatography on its sequence
[#Meek]_, [#Guo1]_. In this model, each amino acid is assigned a number or
a *retention coefficient* (RC) describing its retention properties. The
retention time (RT) during a gradient elution is then calculated as:

.. math::

    RT = \sum_{i=1}^{i=N}{RC_i} + RT_0,

which is a sum of retention coefficients of all amino acid residues in a
polypeptide. This equation can also be expressed in the terms of linear algebra:

.. math::
    
    RT = \overline{aa} \cdot \overline{RC} + RT_0,
    
where :math:`\overline{aa}` is a vector of amino acid composition,
i.e. :math:`\overline{aa}_i` is the number of amino acid residues of i-th type
in a polypeptide; :math:`\overline{RC}` is a vector of respective retention
coefficients.

In this formulation, it is clear that additive model give the same results for
any two peptides that have different sequences but the same amino acid
composition. In other words, **additive model is not sequence-specific**.

However, the additive model has

Several attempts were made in order to improve the accuracy of prediction by
the additive model (for a review of the field we suggest to read [#Baczek]_
and [#Babushok]_). The two implemented in this module are the logarithmic
length correction term described in [#MantLogLen]_ and additional sets of
retention coefficients for terminal amino acid residues [#Tripet]_.

Logarithmic length correction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This enhancement was firstly described in [#MantLogLen]_. Briefly, it was
found that the equation better describing the dependency of RT on peptide
sequence has the following form:

.. math::

    RT = \sum_{i=1}^{i=N}{RC_i} + m\,ln(N) \sum_{i=1}^{i=N}{RC_i} + RT_0

We would call the second term :math:`m\,ln(N) \sum_{i=1}^{i=N}{RC_i}` *the
length correction term* and m - *the length correction factor*. The simplified
and vectorized form of this equation would be:

.. math::
    
    RT = (1 + m\,ln(N)) \, \overline{RC} \cdot \overline{aa} + RT_0

    



References
----------

.. [#Meek] Meek, J. L. Prediction of peptide retention times in high-pressure
   liquid chromatography on the basis of amino acid composition. PNAS, 1980,
   77 (3), 1632-1636.
   `Link <http://www.ncbi.nlm.nih.gov/pubmed/6929513>`_

.. [#Guo1] Guo, D.; Mant, C. T.; Taneja, A. K.; Parker, J. M. R.; Rodges, R. S.
   Prediction of peptide retention times in reversed-phase high-performance
   liquid chromatography I. Determination of retention coefficients of amino
   acid residues of model synthetic peptides. Journal of Chromatography A,
   1986, 359, 499-518.
   `Link. <http://dx.doi.org/10.1016/0021-9673(86)80102-9>`_
   
.. [#Baczek] Baczek, T.; Kaliszan, R. Predictions of peptides' retention times
   in reversed-phase liquid chromatography as a new supportive tool to improve
   protein identification in proteomics. Proteomics, 2009, 9 (4), 835-47.
   `Link. <http://dx.doi.org/10.1002/pmic.200800544>`_

.. [#Babushok] Babushok, V. I.; Zenkevich, I. G. Retention Characteristics of
   Peptides in RP-LC: Peptide Retention Prediction. Chromatographia, 2010, 72
   (9-10), 781-797.
   `Link. <http://dx.doi.org/10.1365/s10337-010-1721-8>`_
   
.. [#MantLogLen] Mant, C. T.; Zhou, N. E.; Hodges, R. S. Correlation of
   protein retention times in reversed-phase chromatography with polypeptide
   chain length and hydrophobicity. Journal of Chromatography A, 1989, 476,
   363-375. `Link. <http://dx.doi.org/10.1016/S0021-9673(01)93882-8>`_

.. [#Tripet] Tripet, B.; Cepeniene, D.; Kovacs, J. M.; Mant, C. T.; Krokhin,
   O. V.; Hodges, R. S. Requirements for prediction of peptide retention time
   in reversed-phase high-performance liquid chromatography:
   hydrophilicity/hydrophobicity of side-chains at the N- and C-termini of
   peptides are dramatically affected by the end-groups and location. Journal
   of chromatography A, 2007, 1141 (2), 212-25.
   `Link. <http://dx.doi.org/10.1016/j.chroma.2006.12.024>`_

"""

# Licensed under the MIT license:
# http://www.opensource.org/licenses/mit-license.php 

import operator
import numpy
from parser import std_labels, peptide_length, amino_acid_composition

def get_RCs(peptides, RTs, length_correction_factor = -0.21,
            labels = std_labels,
            term_aa = False):
    """Calculate the retention coefficients of amino acids given
    retention time of a peptide sample.

    Keyword arguments:

    peptides -- a list of peptide sequences;
    RTs -- a list of retention times of the peptides;
    length_correction_factor -- a multiplier before ln(L) term in the
    equation for the retention time of a peptide;
    labels -- a list of all possible amino acids and terminal
    groups (default 20 standard amino acids, N-terminal NH2- and
    C-terminal -OH);                   
    term_aa -- if True than terminal amino acids are treated as being
    modified with 'ntermX'/'ctermX' modifications. False by default.

    Return a dictionary RC_dict containing the calculated retention
    coefficients.
    RC_dict['aa'] -- the retention coeffitients of amino acids
    RC_dict['const'] -- value of the constant term 
    RC_dict['lcf'] -- the length correction factor

    >>> RCs = get_RCs(['A','AA'], [1.0, 2.0], 0.0, ['A'])
    >>> RCs['const'] = round(RCs['const'], 4) # Rounding for comparison
    >>> RCs == {'aa': {'A': 1.0}, 'lcf': 0.0, 'const': 0.0}
    True
    >>> RCs = get_RCs(['A','AA','B'], [1.0, 2.0, 2.0], 0.0, ['A','B'])
    >>> RCs['aa']['A'] = round(RCs['aa']['A'], 4)
    >>> RCs['aa']['B'] = round(RCs['aa']['B'], 4)
    >>> RCs == {'aa':{'A': 1.0, 'B': 2.0},'const': 0.0, 'lcf': 0.0}
    True
    """

    # Make a list of all amino acids present in the sample.
    peptide_dicts = [
        amino_acid_composition(peptide, False, term_aa, labels=labels)
        for peptide in peptides]

    detected_amino_acids = set([aa for peptide_dict in peptide_dicts
                                for aa in peptide_dict])

    # Determine retention coeffitients using multidimensional fitting.
    composition_array = numpy.array([
        [peptide_dicts[i].get(aa, 0.0) 
         * (1.0 + length_correction_factor
            * numpy.log(peptide_length(peptide_dicts[i])))
           for aa in detected_amino_acids] + [1.0]
        for i in range(len(peptides))])
    
    RCs, res, rank, s = numpy.linalg.lstsq(composition_array,
                                           numpy.array(RTs))
    RC_dict = {}
    RC_dict['aa'] = dict(
        zip(list(detected_amino_acids),
            RCs[:len(detected_amino_acids)]))
    RC_dict['const'] = RCs[len(detected_amino_acids)]
    RC_dict['lcf'] = length_correction_factor

    return RC_dict

def get_RCs_vary_lcf(peptides, RTs,
                labels = std_labels,
                term_aa = False,
                min_lcf = -1.0,
                max_lcf = 1.0):
    """Finds best combination of a length correction factor and
    retention coefficients for given peptide sample.
    Keyword arguments:

    peptides -- a list of peptide sequences;
    RTs -- a list of retention times of the peptides;    
    labels -- a list of all possible amino acids and terminal
    groups (default 20 standard amino acids, N-terminal NH2- and
    C-terminal -OH);                   
    min_lcf -- the minimal value of the length correction factor;
    max_lcf -- the maximal value of the length correction factor.

    Return a dictionary RC_dict containing the calculated retention
    coefficients.
    RC_dict['aa'] -- the retention coeffitients of amino acids
    RC_dict['const'] -- value of the constant term 
    RC_dict['lcf'] -- the length correction factor

    >>> RC_dict = get_RCs_vary_lcf(['A', 'AA', 'AAA'], [1.0, 2.0, 3.0], ['A'])
    >>> RC_dict['aa']['A'] = round(RC_dict['aa']['A'], 4)
    >>> RC_dict['lcf'] = round(RC_dict['lcf'], 4)
    >>> RC_dict['const'] = round(RC_dict['const'], 4)
    >>> RC_dict == {'aa': {'A': 1.0}, 'lcf': 0.0, 'const': 0.0}
    True
    """

    best_r = -1.1
    best_RC_dict = {}
    
    step = (max_lcf - min_lcf) / 10.0
    while step > 0.1:
        lcf_range = numpy.arange(min_lcf, max_lcf,
                                 (max_lcf - min_lcf) / 10.0)
        for lcf in lcf_range:
            RC_dict = get_RCs(peptides, RTs, lcf, labels, term_aa)
            regression_coeffs = linear_regression(
                RTs, 
                [calculate_RT(peptide, RC_dict) for peptide in peptides])
            if regression_coeffs[2] > best_r:
                best_r = regression_coeffs[2]
                best_RC_dict = dict(RC_dict)
        min_lcf = best_RC_dict['lcf'] - step
        max_lcf = best_RC_dict['lcf'] + step
        step = (max_lcf - min_lcf) / 10.0

    return best_RC_dict

def calculate_RT(peptide, RC_dict):
    """Calculate retention time of a peptide using a predetermined set
    of retention coefficients.

    Keyword arguments:
    peptide -- a peptide sequence
    RC_dict -- a set of retention coefficients.

    >>> RT = calculate_RT('AA', {'aa':{'A':1.1},'lcf':0.0,'const':0.1})
    >>> abs(RT - 2.3) < 1e-6      # Float comparison
    True
    >>> RT = calculate_RT('AAA', {'aa': {'ntermA':1.0, 'A':1.1, 'ctermA':1.2},\
                                  'lcf':0.0,\
                                  'const':0.1})
    >>> abs(RT - 3.4) < 1e-6      # Float comparison
    True
    """
    
    amino_acids = [aa for aa in RC_dict['aa']
                   if not (aa.startswith('cterm') or aa.startswith('nterm'))]

    # Check if there are retention coefficients for terminal amino acids.
    term_aa = False
    for aa in RC_dict['aa']:
        if aa.startswith('nterm') or aa.startswith('cterm'):
            term_aa = True
            break

    # Calculate retention time.
    peptide_dict = amino_acid_composition(peptide, False, term_aa,
                                          labels=amino_acids)
    length_correction_term = (
        1.0 + RC_dict['lcf'] * numpy.log(peptide_length(peptide_dict)))
    RT = reduce(operator.add, 
                [peptide_dict[aa] * length_correction_term * RC_dict['aa'][aa]
                 for aa in peptide_dict],
                0.0)
    RT += RC_dict['const']

    return RT

RCs_guo_ph2_0= {'aa':{'K': -2.1,
                      'G': -0.2,
                      'L':  8.1,
                      'A':  2.0,
                      'C':  2.6,
                      'E':  1.1,
                      'D':  0.2,
                      'F':  8.1,
                      'I':  7.4,
                      'H': -2.1,
                      'M':  5.5,
                      'N': -0.6,
                      'Q':  0.0,
                      'P':  2.0,
                      'S': -0.2,
                      'R': -0.6,
                      'T':  0.6,
                      'W':  8.8,
                      'V':  5.0,
                      'Y':  4.5},
                'lcf': 0.0,
                'const': 0.0}


RCs_guo_ph7_0= {'aa':{'K': -0.2,
                      'G': -0.2,
                      'L':  9.0,
                      'A':  2.2,
                      'C':  2.6,
                      'E': -1.3,
                      'D': -2.6,
                      'F':  9.0,
                      'I':  8.3,
                      'H':  2.2,
                      'M':  6.0,
                      'N': -0.8,
                      'Q':  0.0,
                      'P':  2.2,
                      'S': -0.5,
                      'R':  0.9,
                      'T':  0.3,
                      'W':  9.5,
                      'V':  5.7,
                      'Y':  4.6},
                'lcf': 0.0,
                'const': 0.0}

RCs_meek_ph2_1 = {'aa':{'K': -3.2,
                        'G': -0.5,
                        'L': 10.0,
                        'A': -0.1,
                        'C': -2.2,
                        'E': -7.5,
                        'D': -2.8,
                        'F': 13.9,
                        'I': 11.8,
                        'H':  0.8,
                        'M':  7.1,
                        'N': -1.6,
                        'Q': -2.5,
                        'P':  8.0,
                        'S': -3.7,
                        'R': -4.5,
                        'T':  1.5,
                        'W': 18.1,
                        'V':  3.3,
                        'Y':  8.2},
                  'lcf': 0.0,
                  'const': 0.0}

RCs_meek_ph7_4 = {'aa':{'K':  0.1,
                        'G':  0.0,
                        'L':  8.8,
                        'A':  0.5,
                        'C': -6.8,
                        'E':-16.9,
                        'D': -8.2,
                        'F': 13.2,
                        'I': 13.9,
                        'H': -3.5,
                        'M':  4.8,
                        'N':  0.8,
                        'Q': -4.8,
                        'P':  6.1,
                        'S':  1.2,
                        'R':  0.8,
                        'T':  2.7,
                        'W': 14.9,
                        'V':  2.7,
                        'Y':  6.1},
                  'lcf': 0.0,
                  'const': 0.0}

RCs_browne_tfa = {'aa':{'K': -3.7,
                        'G': -1.2,
                        'L': 20.0,
                        'A':  7.3,
                        'C': -9.2,
                        'E': -7.1,
                        'D': -2.9,
                        'F': 19.2,
                        'I':  6.6,
                        'H': -2.1,
                        'M':  5.6,
                        'N': -5.7,
                        'Q': -0.3,
                        'P':  5.1,
                        'S': -4.1,
                        'pS':-6.5,
                        'R': -3.6,
                        'T':  0.8,
                        'pT':-1.6,
                        'W': 16.3,
                        'V':  3.5,
                        'Y':  5.9,
                        'pY': 3.5},
                  'lcf': 0.0,
                  'const': 0.0}

RCs_browne_hfba = {'aa':{'K': -2.5,
                         'G': -2.3,
                         'L': 15.0,
                         'A':  3.9,
                         'C':-14.3,
                         'E': -7.5,
                         'D': -2.8,
                         'F': 14.7,
                         'I': 11.0,
                         'H':  2.0,
                         'M':  4.1,
                         'N': -2.8,
                         'Q':  1.8,
                         'P':  5.6,
                         'S': -3.5,
                         'pS':-7.6,
                         'R':  3.2,
                         'T':  1.1,
                         'pT':-3.0,
                         'W': 17.8,
                         'V':  2.1,
                         'Y':  3.8,
                         'pY':-0.3},
                   'lcf': 0.0,
                   'const': 0.0}

RCs_palmblad = {'aa':{'K': -0.66,
                      'G': -0.29,
                      'L':  2.28,
                      'A':  0.41,
                      'C': -1.32,
                      'E': -0.26,
                      'D':  0.04,
                      'F':  2.68,
                      'I':  2.70,
                      'H':  0.57,
                      'M':  0.98,
                      'N': -0.54,
                      'Q':  1.02,
                      'P':  0.97,
                      'S': -0.71,
                      'R': -0.76,
                      'T':  0.37,
                      'W':  4.68,
                      'V':  2.44,
                      'Y':  2.78},
                'lcf': 0.0,
                'const': 0.0}

RCs_yoshida = {'aa':{'K':  2.77,
                     'G': -0.16,
                     'L': -2.31,
                     'A':  0.28,
                     'C':  0.80,
                     'E':  1.58,
                     'D':  2.45,
                     'F': -2.94,
                     'I': -1.34,
                     'H':  3.44,
                     'M': -0.14,
                     'N':  3.25,
                     'Q':  2.35,
                     'P':  0.77,
                     'S':  2.53,
                     'R':  3.90,
                     'T':  1.73,
                     'W': -1.80,
                     'V': -2.19,
                     'Y': -0.11},
               'lcf': 0.0,
               'const': 0.0}

if __name__ == "__main__":
    import doctest
    doctest.testmod()
