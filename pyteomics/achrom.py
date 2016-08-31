"""
achrom - additive model of polypeptide chromatography
=====================================================

Summary
-------

The additive model of polypeptide chromatography, or achrom, is the most basic
model for peptide retention time prediction. The main equation behind
achrom has the following form:

.. math::

    RT = (1 + m\,ln N) \sum_{i=1}^{i=N}{RC_i n_i} + RT_0


Here, :math:`RC_i` is the retention coefficient of the amino acid
residues of the i-th type, :math:`n_i` corresponds to the number of amino acid
residues of type :math:`i` in the peptide sequence, N is the total number of
different *types* of amino acid residues present,
and :math:`RT_0` is a constant retention time shift.

In order to use achrom, one needs to find the retention
coeffcients, using experimentally determined retention times for a training set
of peptide retention times, i.e. to *calibrate* the model.

Calibration
-----------

  :py:func:`get_RCs` - find a set of retention coefficients using a
  given set of peptides with known retention times and a fixed value of
  length correction parameter.

  :py:func:`get_RCs_vary_lcp` - find the best length correction parameter
  and a set of retention coefficients for a given peptide sample.

Retention time calculation
--------------------------

  :py:func:`calculate_RT` - calculate the retention time of a peptide
  using a given set of retention coefficients.

Data
----

  :py:data:`RCs_guo_ph2_0` - a set of retention coefficients (RCs)
  from [#Guo1]_. Conditions: Synchropak RP-P C18 column (250 x 4.1 mm
  I.D.), gradient (A = 0.1% aq. TFA, pH 2.0; B = 0.1% TFA in acetonitrile) at
  1% B/min, flow rate 1 ml/min, 26 centigrades.

  :py:data:`RCs_guo_ph7_0` - a set of retention coefficients (RCs)
  from [#Guo1]_. Conditions: Synchropak RP-P C18 column (250 x 4.1 mm
  I.D.), gradient (A = aq. 10 mM (NH4)2HPO4 - 0.1 M NaClO4, pH 7.0; B
  = 0.1 M NaClO4 in 60% aq. acetonitrile) at 1.67% B/min, flow rate 1
  ml/min, 26 centigrades.

  :py:data:`RCs_meek_ph2_1` - a set of RCs from [#Meek]_. Conditions: Bio-Rad
  "ODS" column, gradient (A = 0.1 M NaClO4, 0.1% phosphoric acid in
  water; B = 0.1 M NaClO4, 0.1% phosphoric acid in 60%
  aq. acetonitrile) at 1.25% B/min, room temperature.

  :py:data:`RCs_meek_ph7_4` - a set of RCs from [#Meek]_. Conditions: Bio-Rad
  "ODS" column, gradient (A = 0.1 M NaClO4, 5 mM phosphate buffer in
  water; B = 0.1 M NaClO4, 5 mM phosphate buffer in 60%
  aq. acetonitrile) at 1.25% B/min, room temperature.

  :py:data:`RCs_browne_tfa` - a set of RCs found in
  [#Browne]_. Conditions: Waters mjuBondapak C18 column, gradient (A =
  0.1% aq. TFA, B = 0.1% TFA in acetonitrile) at 0.33% B/min, flow
  rate 1.5 ml/min.

  :py:data:`RCs_browne_hfba` - a set of RCs found in
  [#Browne]_. Conditions: Waters mjuBondapak C18 column, gradient (A =
  0.13% aq. HFBA, B = 0.13% HFBA in acetonitrile) at 0.33% B/min, flow
  rate 1.5 ml/min.

  :py:data:`RCs_palmblad` - a set of RCs from
  [#Palmblad]_. Conditions: a fused silica column (80-100 x 0.200 mm
  I.D.) packed in-house with C18 ODS-AQ; solvent A = 0.5% aq. HAc,
  B = 0.5% HAc in acetonitrile.

  :py:data:`RCs_yoshida` - a set of RCs for normal phase chromatography
  from [#Yoshida]_. Conditions:
  TSK gel Amide-80 column (250 x 4.6 mm I.D.), gradient (A = 0.1% TFA
  in ACN-water (90:10); B = 0.1% TFA in ACN-water (55:45)) at 0.6%
  water/min, flow rate 1.0 ml/min, 40 centigrades.

  :py:data:`RCs_yoshida_lc` - a set of length-corrected RCs for normal phase
  chromatography. The set was calculated in [#Moskovets]_ for the data from
  [#Yoshida]_.
  Conditions:
  TSK gel Amide-80 column (250 x 4.6 mm I.D.), gradient (A = 0.1% TFA
  in ACN-water (90:10); B = 0.1% TFA in ACN-water (55:45)) at 0.6%
  water/min, flow rate 1.0 ml/min, 40 centigrades.

  :py:data:`RCs_zubarev` - a set of length-corrected RCs calculated
  on a dataset used in [#Goloborodko]_.
  Conditions: Reprosil-Pur C18-AQ column (150 x 0.075 mm I.D.), gradient (A =
  0.5% AA in water; B = 0.5% AA in ACN-water (90:10)) at
  0.5% water/min, flow rate 200.0 nl/min, room temperature.

  :py:data:`RCs_gilar_atlantis_ph3_0` - a set of retention coefficients obtained
  in [#Gilar]_.
  Conditions: Atlantis HILIC silica column, (150 x 2.1 mm I.D.), 3 um, 100 A,
  gradient (A = water, B = ACN, C = 200 mM ammonium formate):
  0 min, 5% A, 90% B, 5% C; 62.5 min, 55% A, 40% B, 5% C
  at 0.2 ml/min, temperature 40 C, pH 3.0

  :py:data:`RCs_gilar_atlantis_ph4_5` - a set of retention coefficients obtained
  in [#Gilar]_.
  Conditions: Atlantis HILIC silica column, (150 x 2.1 mm I.D.), 3 um, 100 A,
  gradient (A = water, B = ACN, C = 200 mM ammonium formate):
  0 min, 5% A, 90% B, 5% C; 62.5 min, 55% A, 40% B, 5% C
  at 0.2 ml/min, temperature 40 C, pH 4.5

  :py:data:`RCs_gilar_atlantis_ph10_0` - a set of retention coefficients
  obtained in [#Gilar]_.
  Conditions: Atlantis HILIC silica column, (150 x 2.1 mm I.D.), 3 um, 100 A,
  gradient (A = water, B = ACN, C = 200 mM ammonium formate):
  0 min, 5% A, 90% B, 5% C; 62.5 min, 55% A, 40% B, 5% C
  at 0.2 ml/min, temperature 40 C, pH 10.0

  :py:data:`RCs_gilar_beh` - a set of retention coefficients obtained in
  [#Gilar]_.
  Conditions: ACQUITY UPLC BEH HILIC column (150 x 2.1 mm I.D.), 1.7 um, 130 A,
  Mobile phase A: 10 mM ammonium formate buffer, pH 4.5 prepared by
  titrating 10 mM solution of FA with ammonium hydroxide. Mobile phase B:
  90% ACN, 10% mobile phase A (v:v).
  Gradient: 90-60% B in 50 min.

  :py:data:`RCs_gilar_beh_amide` - a set of retention coefficients obtained in
  [#Gilar]_.
  Conditions: ACQUITY UPLC BEH glycan column (150 x 2.1 mm I.D.), 1.7 um, 130 A,
  Mobile phase A: 10 mM ammonium formate buffer, pH 4.5 prepared by
  titrating 10 mM solution of FA with ammonium hydroxide. Mobile phase B:
  90% ACN, 10% mobile phase A (v:v).
  Gradient: 90-60% B in 50 min.

  :py:data:`RCs_gilar_rp` - a set of retention coefficients obtained in
  [#Gilar]_.
  Conditions: ACQUITY UPLC BEH C18 column (100 mm x 2.1 mm I.D.), 1.7 um, 130 A.
  Mobile phase A: 0.02% TFA in water, mobile phase B: 0.018% TFA in ACN.
  Gradient: 0 to 50% B in 50 min, flow rate 0.2 ml/min, temperature 40 C.,
  pH 2.6.

  :py:data:`RCs_krokhin_100A_fa` - a set of retention coefficients obtained in
  [#Krokhin]_.
  Conditions: 300 um x 150mm PepMap100 (Dionex, 0.1% FA), packed with
  5-um Luna C18(2) (Phenomenex, Torrance, CA), pH=2.0.
  Both eluents A (2% ACN in water) and B (98% ACN) contained
  0.1% FA as ion-pairing modifier. 0.33% ACN/min
  linear gradient (0-30% B).

  :py:data:`RCs_krokhin_100A_tfa` - a set of retention coefficients obtained in
  [#Krokhin]_.
  Conditions: 300 um x 150mm PepMap100 (Dionex, 0.1% TFA), packed with
  5-um Luna C18(2) (Phenomenex, Torrance, CA), pH=2.0.
  Both eluents A (2% ACN in water) and B (98% ACN) contained
  0.1% TFA as ion-pairing modifier. 0.33% ACN/min
  linear gradient (0-30% B).

Theory
------

The additive model of polypeptide chromatography, or the model of
retention coefficients was the earliest attempt to describe the dependence of
retention time of a polypeptide in liquid chromatography on its sequence
[#Meek]_, [#Guo1]_. In this model, each amino acid is assigned a number, or
a *retention coefficient* (RC) describing its retention properties. The
retention time (RT) during a gradient elution is then calculated as:

.. math::

    RT = \sum_{i=1}^{i=N}{RC_i \cdot n_i} + RT_0,

which is the sum of retention coefficients of all amino acid residues in a
polypeptide. This equation can also be expressed in terms of linear
algebra:

.. math::

    RT = \\bar{aa} \cdot \\bar{RC} + RT_0,

where :math:`\\bar{aa}` is a vector of amino acid composition,
i.e. :math:`\\bar{aa}_i` is the number of amino acid residues of i-th
type in a polypeptide; :math:`\\bar{RC}` is a vector of respective
retention coefficients.

In this formulation, it is clear that additive model gives the same results for
any two peptides with different sequences but the same amino acid
composition. In other words, **additive model is not sequence-specific**.

The additive model has two advantages over all other models of chromatography
- it is easy to understand and use. The rule behind the additive model is as
simple as it could be: **each amino acid residue shifts retention time by a
fixed value, depending only on its type**. This rule allows geometrical
interpretation. Each peptide may be represented by a point in 21-dimensional
space, with first 20 coordinates equal to the amounts of corresponding amino
acid residues in the peptide and 21-st coordinate equal to RT. The additive
model assumes that a line may be drawn through these points. Of course, this
assumption is valid only partially, and most points would not lie on the
line. But the line would describe the main trend and could be used to estimate
retention time for peptides with known amino acid composition.

This best fit line is described by retention coefficients and :math:`RT_0`.
The procedure of finding these coefficients is called *calibration*. There is
`an analytical solution to calibration of linear models
<http://en.wikipedia.org/wiki/Linear_regression>`_, which makes them
especially useful in real applications.

Several attempts were made in order to improve the accuracy of prediction by
the additive model (for a review of the field we suggest to read [#Baczek]_
and [#Babushok]_). The two implemented in this module are the logarithmic
length correction term described in [#MantLogLen]_ and additional sets of
retention coefficients for terminal amino acid residues [#Tripet]_.

Logarithmic length correction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This enhancement was firstly described in [#MantLogLen]_. Briefly, it was
found that the following equation better describes the dependence of RT on the
peptide sequence:

.. math::

    RT = \sum_{i=1}^{i=N}{RC_i} + m\,ln N \sum_{i=1}^{i=N}{RC_i} + RT_0

We would call the second term :math:`m\,ln N \sum_{i=1}^{i=N}{RC_i}` *the
length correction term* and m - *the length correction parameter*. The
simplified and vectorized form of this equation would be:

.. math::

    RT = (1 + m\,ln N) \, \\bar{RC} \cdot \\bar{aa} + RT_0

This equation may be reduced to a linear form and solved by the standard
methods.

Terminal retention coefficients
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Another significant improvement may be obtained through introduction of
separate sets of retention coefficients for terminal amino acid residues
[#Tripet]_.

References
----------

.. [#Meek] Meek, J. L. `Prediction of peptide retention times in high-pressure
   liquid chromatography on the basis of amino acid composition.
   <http://www.ncbi.nlm.nih.gov/pubmed/6929513>`_
   PNAS, 1980, 77 (3), 1632-1636.

.. [#Guo1] Guo, D.; Mant, C. T.; Taneja, A. K.; Parker, J. M. R.; Hodges,
   R. S.  `Prediction of peptide retention times in reversed-phase
   high-performance liquid chromatography I. Determination of retention
   coefficients of amino acid residues of model synthetic peptides.
   <http://dx.doi.org/10.1016/0021-9673(86)80102-9>`_
   Journal of Chromatography A, 1986, 359, 499-518.

.. [#Baczek] Baczek, T.; Kaliszan, R. `Predictions of peptides' retention times
   in reversed-phase liquid chromatography as a new supportive tool to improve
   protein identification in proteomics.
   <http://dx.doi.org/10.1002/pmic.200800544>`_
   Proteomics, 2009, 9 (4), 835-47.

.. [#Babushok] Babushok, V. I.; Zenkevich, I. G. `Retention Characteristics of
   Peptides in RP-LC: Peptide Retention Prediction.
   <http://dx.doi.org/10.1365/s10337-010-1721-8>`_
   Chromatographia, 2010, 72 (9-10), 781-797.

.. [#MantLogLen] Mant, C. T.; Zhou, N. E.; Hodges, R. S. `Correlation of
   protein retention times in reversed-phase chromatography with polypeptide
   chain length and hydrophobicity.
   <http://dx.doi.org/10.1016/S0021-9673(01)93882-8>`_
   Journal of Chromatography A, 1989, 476, 363-375.

.. [#Tripet] Tripet, B.; Cepeniene, D.; Kovacs, J. M.; Mant, C. T.; Krokhin,
   O. V.; Hodges, R. S. `Requirements for prediction of peptide retention time
   in reversed-phase high-performance liquid chromatography:
   hydrophilicity/hydrophobicity of side-chains at the N- and C-termini of
   peptides are dramatically affected by the end-groups and location.
   <http://dx.doi.org/10.1016/j.chroma.2006.12.024>`_
   Journal of chromatography A, 2007, 1141 (2), 212-25.

.. [#Browne] Browne, C. A.; Bennett, H. P. J.; Solomon, S. `The
   isolation of peptides by high-performance liquid chromatography
   using predicted elution positions
   <http://www.sciencedirect.com/science/article/pii/000326978290238X>`_.
   Analytical Biochemistry, 1982, 124 (1), 201-208.

.. [#Palmblad] Palmblad, M.; Ramstrom, M.; Markides, K. E.; Hakansson,
   P.; Bergquist, J. `Prediction of Chromatographic Retention and
   Protein Identification in Liquid Chromatography/Mass
   Spectrometry
   <http://pubs.acs.org/doi/abs/10.1021/ac0256890>`_.
   Analytical Chemistry, 2002, 74 (22), 5826-5830.

.. [#Yoshida] Yoshida, T. Calculation of peptide retention
   coefficients in normal-phase liquid chromatography. Journal of
   Chromatography A, 1998, 808 (1-2), 105-112.

.. [#Moskovets] Moskovets, E.; Goloborodko A. A.; Gorshkov A. V.; Gorshkov M.V.
   `Limitation of predictive 2-D liquid chromatography in reducing the database
   search space in shotgun proteomics: In silico studies.
   <http://dx.doi.org/10.1002/jssc.201100798>`_
   Journal of Separation Science, 2012, 35 (14), 1771-1778.

.. [#Goloborodko] Goloborodko A. A.; Mayerhofer C.; Zubarev A. R.;
   Tarasova I. A.; Gorshkov A. V.; Zubarev, R. A.; Gorshkov, M. V.
   `Empirical approach to false discovery rate
   estimation in shotgun proteomics. <http://dx.doi.org/10.1002/rcm.4417>`_
   Rapid communications in mass spectrometry, 2010, 24(4), 454-62.

.. [#Gilar] Gilar, M., & Jaworski, A. (2011). `Retention behavior of peptides in
    hydrophilic-interaction chromatography.
    <http://dx.doi.org/10.1016/j.chroma.2011.04.005>`_
    Journal of chromatography A, 1218(49), 8890-6.

.. [#Krokhin] Dwivedi, R. C.; Spicer, V.; Harder, M.; Antonovici, M.; Ens, W.;
    Standing, K. G.; Wilkins, J. A.; Krokhin, O. V. (2008). `Practical
    implementation of 2D HPLC scheme with accurate peptide retention prediction
    in both dimensions for high-throughput bottom-up proteomics
    <http://pubs.acs.org/doi/abs/10.1021/ac800984n>`_.
    Analytical Chemistry, 80(18), 7036-42.

Dependencies
------------

This module requires :py:mod:`numpy`.

--------------------------------------------------------------------------------
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

import numpy as np
from .auxiliary import linear_regression, PyteomicsError
from . import parser

def get_RCs(sequences, RTs, lcp = -0.21,
            term_aa = False, **kwargs):
    """Calculate the retention coefficients of amino acids using
    retention times of a peptide sample and a fixed value of length
    correction parameter.

    Parameters
    ----------
    sequences : list of str
        List of peptide sequences.
    RTs: list of float
        List of corresponding retention times.
    lcp : float, optional
        A multiplier before ln(L) term in the equation for the retention
        time of a peptide. Set to -0.21 by default.
    term_aa : bool, optional
        If :py:const:`True`, terminal amino acids are treated as being
        modified with 'ntermX'/'ctermX' modifications. :py:const:`False`
        by default.
    labels : list of str, optional
        List of all possible amino acids and terminal groups
        If not given, any modX labels are allowed.

    Returns
    -------
    RC_dict : dict
        Dictionary with the calculated retention coefficients.

        - RC_dict['aa'] -- amino acid retention coefficients.

        - RC_dict['const'] -- constant retention time shift.

        - RC_dict['lcp'] -- length correction parameter.

    Examples
    --------
    >>> RCs = get_RCs(['A','AA'], [1.0, 2.0], 0.0, labels=['A'])
    >>> abs(RCs['aa']['A'] - 1) < 1e-6 and abs(RCs['const']) < 1e-6
    True
    >>> RCs = get_RCs(['A','AA','B'], [1.0, 2.0, 2.0], 0.0, labels=['A','B'])
    >>> abs(RCs['aa']['A'] - 1) + abs(RCs['aa']['B'] - 2) + \
            abs(RCs['const']) < 1e-6
    True
    """

    labels = kwargs.get('labels')

    # Make a list of all amino acids present in the sample.
    peptide_dicts = [
            parser.amino_acid_composition(peptide, False, term_aa,
                               allow_unknown_modifications=True,
                               labels=labels)
            if not isinstance(peptide, dict) else peptide
        for peptide in sequences]

    detected_amino_acids = {aa for peptide_dict in peptide_dicts
                                for aa in peptide_dict}

    # Determine retention coefficients using multidimensional linear
    # regression.
    composition_array = []
    for pdict in peptide_dicts:
        loglen = np.log(parser.length(pdict))
        composition_array.append([pdict.get(aa, 0.)
             * (1. + lcp * loglen)
               for aa in detected_amino_acids] + [1.])

    # Add normalizing conditions for terminal retention coefficients. The
    # condition we are using here is quite arbitrary. It implies that the sum
    # of N- or C-terminal RCs minus the sum of corresponding internal RCs must
    # be equal to zero.
    if term_aa:
        for term_label in ['nterm', 'cterm']:
            normalizing_peptide = []
            for aa in detected_amino_acids:
                if aa.startswith(term_label):
                    normalizing_peptide.append(1.0)
                elif (term_label+aa) in detected_amino_acids:
                    normalizing_peptide.append(-1.0)
                else:
                    normalizing_peptide.append(0.0)
            normalizing_peptide.append(0.0)
            composition_array.append(normalizing_peptide)
            RTs.append(0.0)

    # Use least square linear regression.
    RCs, res, rank, s = np.linalg.lstsq(np.array(composition_array),
                                           np.array(RTs))

    # Remove normalizing elements from the RTs vector.
    if term_aa:
        for term_label in ['nterm', 'cterm']:
            RTs.pop()

    # Form output.
    RC_dict = {}
    RC_dict['aa'] = dict(
        zip(list(detected_amino_acids),
            RCs[:len(detected_amino_acids)]))
    RC_dict['aa'][parser.std_nterm] = 0.0
    RC_dict['aa'][parser.std_cterm] = 0.0
    RC_dict['const'] = RCs[len(detected_amino_acids)]
    RC_dict['lcp'] = lcp

    # Find remaining terminal RCs.
    if term_aa:
        for term_label in ['nterm', 'cterm']:
            # Check if there are terminal RCs remaining undefined.
            undefined_term_RCs = [aa for aa in RC_dict['aa']
                                if aa[1:5] != 'term'
                                and term_label + aa not in RC_dict['aa']]
            if not undefined_term_RCs:
                continue

            # Find a linear relationship between internal and terminal RCs.
            defined_term_RCs = [aa for aa in RC_dict['aa']
                              if aa[1:5] != 'term'
                              and term_label + aa in RC_dict['aa']]

            a, b, r, stderr = linear_regression(
                [RC_dict['aa'][aa] for aa in defined_term_RCs],
                [RC_dict['aa'][term_label+aa] for aa in defined_term_RCs])

            # Define missing terminal RCs using this linear equation.
            for aa in undefined_term_RCs:
                RC_dict['aa'][term_label + aa] = a * RC_dict['aa'][aa] + b

    return RC_dict

def get_RCs_vary_lcp(sequences, RTs,
                term_aa = False,
                lcp_range = (-1.0, 1.0),
                **kwargs):
    """Find the best combination of a length correction parameter and
    retention coefficients for a given peptide sample.

    Parameters
    ----------
    sequences : list of str
        List of peptide sequences.
    RTs : list of float
        List of corresponding retention times.
    term_aa : bool, optional
        If True, terminal amino acids are treated as being
        modified with 'ntermX'/'ctermX' modifications. False by default.
    lcp_range : 2-tuple of float, optional
        Range of possible values of the length correction parameter.
    labels : list of str, optional
        List of labels for all possible amino acids and terminal groups
        If not given, any modX labels are allowed.
    lcp_accuracy : float, optional
        The accuracy of the length correction parameter calculation.

    Returns
    -------
    RC_dict : dict
        Dictionary with the calculated retention coefficients.

        - RC_dict['aa'] -- amino acid retention coefficients.

        - RC_dict['const'] -- constant retention time shift.

        - RC_dict['lcp'] -- length correction parameter.

    Examples
    --------
    >>> RCs = get_RCs_vary_lcp(['A', 'AA', 'AAA'], \
        [1.0, 2.0, 3.0], \
        labels=['A'])
    >>> abs(RCs['aa']['A'] - 1) + abs(RCs['lcp']) + abs(RCs['const']) < 1e-6
    True
    """
    labels = kwargs.get('labels')

    best_r = -1.1
    best_RC_dict = {}
    lcp_accuracy = kwargs.get('lcp_accuracy', 0.1)

    min_lcp = lcp_range[0]
    max_lcp = lcp_range[1]
    step = (max_lcp - min_lcp) / 10.0
    peptide_dicts = [
            parser.amino_acid_composition(peptide, False, term_aa,
                                      allow_unknown_modifications=True,
                                      labels=labels)
            if not isinstance(peptide, dict) else peptide
        for peptide in sequences]
    while step > lcp_accuracy:
        lcp_grid = np.arange(min_lcp, max_lcp,
                                (max_lcp - min_lcp) / 10.0)
        for lcp in lcp_grid:
            RC_dict = get_RCs(peptide_dicts, RTs, lcp, term_aa, labels=labels)
            regression_coeffs = linear_regression(
                RTs,
                [calculate_RT(peptide, RC_dict) for peptide in peptide_dicts])
            if regression_coeffs[2] > best_r:
                best_r = regression_coeffs[2]
                best_RC_dict = dict(RC_dict)
        min_lcp = best_RC_dict['lcp'] - step
        max_lcp = best_RC_dict['lcp'] + step
        step = (max_lcp - min_lcp) / 10.0

    return best_RC_dict

def calculate_RT(peptide, RC_dict, raise_no_mod=True):
    """Calculate the retention time of a peptide using a given set
    of retention coefficients.

    Parameters
    ----------
    peptide : str or dict
        A peptide sequence or amino acid composition.
    RC_dict : dict
        A set of retention coefficients, length correction parameter and
        a fixed retention time shift. Keys are: 'aa', 'lcp' and 'const'.
    raise_no_mod : bool, optional
        If :py:const:`True` then an exception is raised when a modified amino
        acid from `peptides` is not found in `RC_dict`. If :py:const:`False`,
        then the retention coefficient for the non-modified amino acid residue
        is used instead. :py:const:`True` by default.

    Returns
    -------
    RT : float
        Calculated retention time.

    Examples
    --------
    >>> RT = calculate_RT('AA', {'aa': {'A': 1.1}, 'lcp':0.0, 'const': 0.1})
    >>> abs(RT - 2.3) < 1e-6      # Float comparison
    True
    >>> RT = calculate_RT('AAA', {'aa': {'ntermA': 1.0, 'A': 1.1, 'ctermA': 1.2},\
        'lcp': 0.0, 'const':0.1})
    >>> abs(RT - 3.4) < 1e-6      # Float comparison
    True
    >>> RT = calculate_RT({'A': 3}, {'aa': {'ntermA': 1.0, 'A': 1.1, 'ctermA': 1.2},\
        'lcp': 0.0, 'const':0.1})
    >>> abs(RT - 3.4) < 1e-6      # Float comparison
    True
    """

    amino_acids = [aa for aa in RC_dict['aa']
                   if not (aa[:5] == 'nterm' or aa[:5] == 'cterm')]

    # Check if there are retention coefficients for terminal amino acids.
    term_aa = False
    for aa in RC_dict['aa']:
        if aa[:5] == 'nterm' or aa[:5] == 'cterm':
            term_aa = True
            break

    # Calculate retention time.
    if isinstance(peptide, dict):
        peptide_dict = peptide
    else:
        peptide_dict = parser.amino_acid_composition(peptide, False, term_aa,
            allow_unknown_modifications=True, labels=amino_acids)
    RT = 0.0
    for aa in peptide_dict:
        if aa not in RC_dict['aa']:
            if len(aa) == 1:
                raise PyteomicsError('No RC for residue "{}".'.format(aa))
            if (not raise_no_mod) and aa[-1] in RC_dict['aa']:
                RT += peptide_dict[aa] * RC_dict['aa'][aa[-1]]
            else:
                raise PyteomicsError(
                    'Residue "{0}" not found in RC_dict. '.format(aa) +
                    'Set raise_no_mod=False to ignore this error ' +
                    'and use the RC for "{0}"" instead.'.format(aa[-1]))
        else:
            RT += peptide_dict[aa] * RC_dict['aa'][aa]

    length_correction_term = (
        1.0 + RC_dict.get('lcp', 0) * np.log(parser.length(peptide_dict)))
    RT *= length_correction_term

    RT += RC_dict.get('const', 0)

    return RT

RCs_guo_ph2_0 = {'aa':{'K': -2.1,
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
                       'Y':  4.5,
                       'H-': 0.0,
                       '-OH':0.0},
                 'lcp': 0.0,
                 'const': 0.0}
"""A set of retention coefficients from Guo, D.; Mant, C. T.; Taneja,
A. K.; Parker, J. M. R.; Hodges, R. S.  Prediction of peptide
retention times in reversed-phase high-performance liquid
chromatography I. Determination of retention coefficients of amino
acid residues of model synthetic peptides. Journal of Chromatography
A, 1986, 359, 499-518.

Conditions: Synchropak RP-P C18 column (250 x 4.1 mm I.D.), gradient
(A = 0.1% aq. TFA, pH 2.0; B = 0.1% TFA in acetonitrile) at 1% B/min,
flow rate 1 ml/min, 26 centigrades.
"""

RCs_guo_ph7_0 = {'aa':{'K': -0.2,
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
                       'Y':  4.6,
                       'H-': 0.0,
                       '-OH':0.0},
                 'lcp': 0.0,
                 'const': 0.0}
"""A set of retention coefficients from Guo, D.; Mant, C. T.; Taneja,
A. K.; Parker, J. M. R.; Hodges, R. S.  Prediction of peptide
retention times in reversed-phase high-performance liquid
chromatography I. Determination of retention coefficients of amino
acid residues of model synthetic peptides. Journal of Chromatography
A, 1986, 359, 499-518.

Conditions: Synchropak RP-P C18 column (250 x 4.1 mm I.D.), gradient
(A = aq. 10 mM (NH4)2HPO4 - 0.1 M NaClO4, pH 7.0; B = 0.1 M NaClO4 in
60% aq. acetonitrile) at 1.67% B/min, flow rate 1 ml/min, 26
centigrades.
"""

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
                        'Y':  8.2,
                        'H-': 0.0,
                        '-OH':0.0},
                  'lcp': 0.0,
                  'const': 0.0}
"""A set of retention coefficients determined in Meek,
J. L. Prediction of peptide retention times in high-pressure liquid
chromatography on the basis of amino acid composition. PNAS, 1980, 77
(3), 1632-1636.

.. note :: C stands for Cystine.

Conditions: Bio-Rad "ODS" column, gradient (A = 0.1 M NaClO4,
0.1% phosphoric acid in water; B = 0.1 M NaClO4, 0.1% phosphoric acid
in 60% aq. acetonitrile) at 1.25% B/min, room temperature.
"""

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
                        'Y':  6.1,
                        'H-': 0.0,
                        '-OH':0.0},
                  'lcp': 0.0,
                  'const': 0.0}
"""A set of retention coefficients determined in Meek,
J. L. Prediction of peptide retention times in high-pressure liquid
chromatography on the basis of amino acid composition. PNAS, 1980, 77
(3), 1632-1636.

.. note :: C stands for Cystine.

Conditions: Bio-Rad "ODS" column, gradient (A = 0.1 M NaClO4,
5 mM phosphate buffer in water; B = 0.1 M NaClO4, 5 mM phosphate buffer
in 60% aq. acetonitrile) at 1.25% B/min, room temperature.
"""

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
                        'pY': 3.5,
                        'H-': 0.0,
                        '-OH':0.0},
                  'lcp': 0.0,
                  'const': 0.0}
"""A set of retention coefficients determined in Browne, C. A.;
Bennett, H. P. J.; Solomon, S. The isolation of peptides by
high-performance liquid chromatography using predicted elution
positions. Analytical Biochemistry, 1982, 124 (1), 201-208.

Conditions: Waters mjuBondapak C18 column, gradient (A = 0.1% aq. TFA,
B = 0.1% TFA in acetonitrile) at 0.33% B/min, flow rate 1.5 ml/min.
"""

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
                         'pY':-0.3,
                         'H-': 0.0,
                         '-OH':0.0},
                   'lcp': 0.0,
                   'const': 0.0}
"""A set of retention coefficients determined in Browne, C. A.;
Bennett, H. P. J.; Solomon, S. The isolation of peptides by
high-performance liquid chromatography using predicted elution
positions. Analytical Biochemistry, 1982, 124 (1), 201-208.

Conditions: Waters mjuBondapak C18 column, gradient (A = 0.13% aq. HFBA,
B = 0.13% HFBA in acetonitrile) at 0.33% B/min, flow rate 1.5 ml/min.
"""

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
                      'Y':  2.78,
                      'H-': 0.0,
                      '-OH':0.0},
                'lcp': 0.0,
                'const': 0.0}
"""A set of retention coefficients determined in Palmblad, M.;
Ramstrom, M.; Markides, K. E.; Hakansson, P.; Bergquist, J. Prediction
of Chromatographic Retention and Protein Identification in Liquid
Chromatography/Mass Spectrometry. Analytical Chemistry, 2002, 74 (22),
5826-5830.

Conditions: a fused silica column (80-100 x 0.200 mm I.D.) packed
in-house with C18 ODS-AQ; solvent A = 0.5% aq. HAc, B = 0.5% HAc in
acetonitrile.
"""

RCs_yoshida = {'aa':{'K':  2.77,
                     'G': -0.16,
                     'L': -2.31,
                     'A':  0.28,
                     'C':  0.80,
                  'camC':  0.80,
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
                     'Y': -0.11,
                     'H-': 0.0,
                     '-OH':0.0},
               'lcp': 0.0,
               'const': 0.0}
"""A set of retention coefficients determined in Yoshida,
T. Calculation of peptide retention coefficients in normal-phase
liquid chromatography. Journal of Chromatography A, 1998, 808 (1-2),
105-112.

.. note::  Cysteine is Carboxymethylated.

Conditions: TSK gel Amide-80 column (250 x 4.6 mm I.D.), gradient (A =
0.1% TFA in ACN-water (90:10); B = 0.1% TFA in ACN-water (55:45)) at
0.6% water/min, flow rate 1.0 ml/min, 40 centigrades.
"""

RCs_yoshida_lc = {'aa': {'A': 1.29,
                         'C': 0.94,
                      'camC': 0.94,
                         'D': 3.89,
                         'E': 4.40,
                         'F': -4.18,
                         'G': 1.29,
                         'H': 7.57,
                         'I': -2.65,
                         'K': 7.33,
                         'L': -3.93,
                         'M': -1.48,
                         'N': 6.65,
                         'P': 1.03,
                         'Q': 6.68,
                         'R': 7.08,
                         'S': 5.09,
                         'T': 3.46,
                         'V': -2.52,
                         'W': -1.87,
                         'Y': -0.46,
                         'H-': 0.0,
                         '-OH': 0.0},
                  'const': 0.0,
                  'lcp': -0.2}
"""A set of retention coefficients from the length-corrected model
of normal-phase peptide chromatography. The dataset comes from Yoshida, T.
Calculation of peptide retention coefficients in normal-phase
liquid chromatography. Journal of Chromatography A, 1998, 808 (1-2),
105-112. The RCs were calculated in Moskovets, E.; Goloborodko A. A.;
Gorshkov A. V.; Gorshkov M.V. Limitation of predictive 2-D liquid chromatography
in reducing the database search space in shotgun proteomics: In silico studies.
Journal of Separation Science, 2012, 35 (14), 1771-1778.

.. note::  Cysteine is Carboxymethylated.

Conditions: TSK gel Amide-80 column (250 x 4.6 mm I.D.), gradient (A =
0.1% TFA in ACN-water (90:10); B = 0.1% TFA in ACN-water (55:45)) at
0.6% water/min, flow rate 1.0 ml/min, 40 centigrades.
"""

RCs_zubarev = {'aa': {'A': 6.73,
                      'E': 5.66,
                      'C': 3.25,
                      'D': 5.64,
                      'G': 2.35,
                      'F': 27.43,
                      'I': 20.50,
                      'H': -0.66,
                      'K': -4.47,
                      'M': 17.39,
                      'L': 23.38,
                      'N': 2.57,
                      'Q': 2.93,
                      'P': 5.66,
                      'S': 3.58,
                      'R': -2.55,
                      'T': 4.88,
                      'Y': 13.22,
                      'W': 31.27,
                      'V': 13.05,
                   'camC': 3.25,
                      'C': 3.25,
                    'oxM': -7.61,
                    '-OH': 0.0,
                     'H-': 0.0},
            'const': 0.53,
              'lcp': -0.21}
"""A set of retention coefficients from the length-corrected model
of reversed-phase peptide chromatography. The dataset was taken from
Goloborodko A. A.; Mayerhofer C.; Zubarev A. R.; Tarasova I. A.; Gorshkov A. V.;
Zubarev, R. A.; Gorshkov, M. V. Empirical approach to false discovery rate
estimation in shotgun proteomics. Rapid communications in mass spectrometry,
2010, 24(4), 454-62.

.. note::  Cysteine is Carbamidomethylated.

Conditions: Reprosil-Pur C18-AQ column (150 x 0.075 mm I.D.), gradient (A =
0.5% AA in water; B = 0.5% AA in ACN-water (90:10)) at
0.5% water/min, flow rate 200.0 nl/min, room temperature.
"""

RCs_gilar_atlantis_ph3_0 = {'aa': {'K': 15.90,
    'R': 13.64,
    'H': 12.94,
    'E': 2.97,
    'P': 4.77,
    'Q': 5.43,
    'D': 3.20,
   'C*': 4.87,
    'C': 4.87,
    'N': 3.91,
    'A': 3.34,
    'G': 3.33,
    'S': 3.04,
    'T': 2.71,
    'V': 1.75,
    'I': 0.65,
    'M': 1.13,
    'L': 0.13,
    'F': -1.17,
    'Y': -0.22,
    'W': -2.47},
        'lcp': 0.0,
        'const': 21.33}
"""A set of retention coefficients for normal phase chromatography obtained in
Gilar, M., & Jaworski, A. (2011). Retention behavior of peptides in
hydrophilic-interaction chromatography. Journal of chromatography A, 1218(49),
8890-6.

.. note::  Cysteine is Carbamidomethylated.

Conditions: Atlantis HILIC silica column (150 x 2.1 mm I.D.), 3 um, 100 A,
gradient (A = water, B = ACN, C = 200 mM ammonium formate):
0 min, 5% A, 90% B, 5% C; 62.5 min, 55% A, 40% B, 5% C
at 0.2 ml/min, temperature 40 C, pH 3.0"""

RCs_gilar_atlantis_ph4_5 = {'aa': {'K': 15.49,
    'R': 13.33,
    'H': 12.19,
    'E': 6.93,
    'P': 5.89,
    'Q': 5.68,
    'D': 5.31,
   'C*': 5.23,
    'C': 5.23,
    'N': 4.07,
    'A': 3.6,
    'G': 3.46,
    'S': 2.62,
    'T': 2.33,
    'V': 1.42,
    'I': 0.84,
    'M': 0.34,
    'L': 0.29,
    'F': -1.21,
    'Y': -1.62,
    'W': -2.08},
        'lcp': 0.0,
        'const': 23.95}
"""A set of retention coefficients for normal phase chromatography obtained in
Gilar, M., & Jaworski, A. (2011). Retention behavior of peptides in
hydrophilic-interaction chromatography. Journal of chromatography A, 1218(49),
8890-6.

.. note::  Cysteine is Carbamidomethylated.

Conditions: Atlantis HILIC silica column (150 x 2.1 mm I.D.), 3 um, 100 A,
gradient (A = water, B = ACN, C = 200 mM ammonium formate):
0 min, 5% A, 90% B, 5% C; 62.5 min, 55% A, 40% B, 5% C
at 0.2 ml/min, temperature 40 C, pH 4.5"""

RCs_gilar_atlantis_ph10_0 = {'aa': {'K': 25.23,
    'R': 23.38,
    'H': 5.94,
    'E': 0.59,
    'P': 4.00,
    'Q': 3.53,
    'D': -0.84,
   'C*': 3.52,
    'C': 3.52,
    'N': 3.26,
    'A': 3.64,
    'G': 3.02,
    'S': 2.28,
    'T': 1.74,
    'V': 1.05,
    'I': 1.51,
    'M': -0.61,
    'L': 0.25,
    'F': -0.17,
    'Y': -0.79,
    'W': 0.23},
        'lcp': 0.0,
        'const': 13.78}
"""A set of retention coefficients for normal phase chromatography obtained in
Gilar, M., & Jaworski, A. (2011). Retention behavior of peptides in
hydrophilic-interaction chromatography. Journal of chromatography A, 1218(49),
8890-6.

.. note::  Cysteine is Carbamidomethylated.

Conditions: Atlantis HILIC silica column (150 x 2.1 mm I.D.), 3 um, 100 A,
gradient (A = water, B = ACN, C = 200 mM ammonium formate):
0 min, 5% A, 90% B, 5% C; 62.5 min, 55% A, 40% B, 5% C
at 0.2 ml/min, temperature 40 C, pH 10.0"""

RCs_gilar_beh = {'aa': {'K': 9.49,
    'R': 8.56,
    'H': 8.40,
    'E': 5.95,
    'P': 4.73,
    'Q': 4.65,
    'D': 4.97,
    'C': 3.47,
   'C*': 3.47,
    'N': 3.50,
    'A': 2.90,
    'G': 2.63,
    'S': 2.14,
    'T': 2.19,
    'V': 1.71,
    'I': 1.30,
    'M': 1.40,
    'L': 0.73,
    'F': -0.09,
    'Y': -0.40,
    'W': 0.11},
        'lcp': 0.0,
        'const': 18.41}
"""A set of retention coefficients for normal phase chromatography obtained in
Gilar, M., & Jaworski, A. (2011). Retention behavior of peptides in
hydrophilic-interaction chromatography. Journal of chromatography A, 1218(49),
8890-6.

.. note::  Cysteine is Carbamidomethylated.

Conditions: ACQUITY UPLC BEH HILIC column (150 x 2.1 mm I.D.), 1.7 um, 130 A,
Mobile phase A: 10 mM ammonium formate buffer, pH 4.5 prepared by
titrating 10 mM solution of FA with ammonium hydroxide. Mobile phase B:
90% ACN, 10% mobile phase A (v:v).
Gradient: 90-60% B in 50 min."""

RCs_gilar_beh_amide = {'aa': {'K': 7.19,
    'R': 6.68,
    'H': 6.16,
    'E': 6.11,
    'P': 3.18,
    'Q': 5.19,
    'D': 6.02,
   'C*': 3.71,
    'C': 3.71,
    'N': 4.16,
    'A': 2.64,
    'G': 3.12,
    'S': 3.17,
    'T': 3.41,
    'V': 0.83,
    'I': -0.69,
    'M': -0.12,
    'L': -1.24,
    'F': -1.93,
    'Y': 0.46,
    'W': -2.11},
        'lcp': 0.0,
        'const': 24.26}
"""A set of retention coefficients for normal phase chromatography obtained in
Gilar, M., & Jaworski, A. (2011). Retention behavior of peptides in
hydrophilic-interaction chromatography. Journal of chromatography A, 1218(49),
8890-6.

.. note::  Cysteine is Carbamidomethylated.

Conditions: ACQUITY UPLC BEH glycan column (150 x 2.1 mm I.D.), 1.7 um, 130 A,
Mobile phase A: 10 mM ammonium formate buffer, pH 4.5 prepared by
titrating 10 mM solution of FA with ammonium hydroxide. Mobile phase B:
90% ACN, 10% mobile phase A (v:v).
Gradient: 90-60% B in 50 min."""

RCs_gilar_rp = {'aa': {'K': -1.015,
    'R': -0.681,
    'H': -1.937,
    'E': 1.475,
    'P': 3.496,
    'Q': 1.228,
    'D': 1.326,
    'C': 1.832,
   'C*': 1.832,
    'N': 0.299,
    'A': 2.322,
    'G': 1.172,
    'S': 1.165,
    'T': 1.894,
    'V': 5.695,
    'I': 8.343,
    'M': 5.128,
    'L': 9.069,
    'F': 10.877,
    'Y': 5.603,
    'W': 12.183},
        'lcp': 0.0,
        'const': -3.696}
"""A set of retention coefficients for normal phase chromatography obtained in
Gilar, M., & Jaworski, A. (2011). Retention behavior of peptides in
hydrophilic-interaction chromatography. Journal of chromatography A, 1218(49),
8890-6.

.. note::  Cysteine is Carbamidomethylated.

Conditions: ACQUITY UPLC BEH C18 column (100 mm x 2.1 mm I.D.), 1.7 um, 130 A.
Mobile phase A: 0.02% TFA in water, mobile phase B: 0.018% TFA in ACN.
Gradient: 0 to 50% B in 50 min, flow rate 0.2 ml/min, temperature 40 C., pH 2.6.
"""

RCs_krokhin_100A_fa = {'aa':{'K': -5.08,
                       'G': -0.07,
                       'L':  9.89,
                       'A':  1.63,
                       'C':  0.7,
                    'camC':  0.7,
                       'E':  1.75,
                       'D':  0.95,
                       'F':  11.92,
                       'I':  9.06,
                       'H': -5.05,
                       'M':  6.96,
                       'N': -0.59,
                       'Q':  0.2,
                       'P':  1.98,
                       'S': 0.27,
                       'R': -3.55,
                       'T':  1.37,
                       'W':  13.67,
                       'V':  5.72,
                       'Y':  5.97},
                 'lcf': 0.0,
                 'const': 0.0}
"""A set of retention coefficients from R.C. Dwivedi, V. Spicer,
M. Harder, M. Antonovici, W. Ens, K.G. Standing, J.A. Wilkins, and O.V. Krokhin;
Analytical Chemistry 2008 80 (18), 7036-7042.
Practical Implementation of 2D HPLC Scheme with Accurate Peptide
Retention Prediction in Both Dimensions for High-Throughput Bottom-Up Proteomics.

.. note::  Cysteine is Carbamidomethylated.

Conditions: 300 um x 150mm PepMap100 (Dionex, 0.1% FA), packed with
5-um Luna C18(2) (Phenomenex, Torrance, CA), pore size 100A,  pH=2.0.
Both eluents A (2% ACN in water) and B (98% ACN) contained
0.1% FA as ion-pairing modifier. 0.33% ACN/min
linear gradient (0-30% B).
"""

RCs_krokhin_100A_tfa = {'aa':{'K': -3.53,
                       'G': -0.35,
                       'L':  9.44,
                       'A':  1.11,
                       'C':  0.04,
                    'camC':  0.04,
                       'E':  1.08,
                       'D':  -0.22,
                       'F':  11.34,
                       'I':  7.86,
                       'H': -3.04,
                       'M':  6.57,
                       'N': -1.44,
                       'Q':  -0.53,
                       'P':  1.62,
                       'S': -0.33,
                       'R': -2.58,
                       'T':  0.48,
                       'W':  13.12,
                       'V':  4.86,
                       'Y':  5.4},
                 'lcf': 0.0,
                 'const': 0.0}
"""A set of retention coefficients from R.C. Dwivedi, V. Spicer,
M. Harder, M. Antonovici, W. Ens, K.G. Standing, J.A. Wilkins, and O.V. Krokhin;
Analytical Chemistry 2008 80 (18), 7036-7042.
Practical Implementation of 2D HPLC Scheme with Accurate Peptide
Retention Prediction in Both Dimensions for High-Throughput Bottom-Up Proteomics.

.. note :: Cysteine is Carbamidomethylated.

Conditions: 300 um x 150mm PepMap100 (Dionex, 0.1% TFA), packed with
5-um Luna C18(2) (Phenomenex, Torrance, CA), pore size 100 A,  pH=2.0.
Both eluents A (2% ACN in water) and B (98% ACN) contained
0.1% TFA as ion-pairing modifier. 0.33% ACN/min
linear gradient (0-30% B).
"""


if __name__ == "__main__":
    import doctest
    doctest.testmod()
