Retention time prediction
=========================

Pyteomics has two modules for prediction of retention times (RTs) of peptides 
and proteins in liquid chromatography.

BioLCCC
-------

The first module is :py:mod:`pyteomics.biolccc`. This module implements
the BioLCCC model of liquid chromatography of polypeptides. 
:py:mod:`pyteomics.biolccc` is not distributed with the main package and has 
to be installed separately. :py:mod:`pyteomics.biolccc` can be downloaded from 
http://pypi.python.org/pypi/pyteomics.biolccc, and the project documentation
is hosted at http://packages.python.org/pyteomics.biolccc.

Additive model of peptide chromatography
----------------------------------------

Another option for retention time prediction is the :py:mod:`achrom` module
distributed with Pyteomics. It implements the additive model of polypeptide
chromatography. Briefly, in the additive model each amino acid residue changes 
retention time by a fixed value, depending only on its type (e.g. an alanine r
esidue add 2.0 mins to RT, while an arginine decreases it by 1.1 min). The module 
documentation contains the complete description of this model and the references. 
In this tutorial we will focus on the basic usage.

Retention time prediction
.........................

Retention time prediction with :py:mod:`achrom` is done by the
:py:func:`calculate_RT` function:

.. code-block:: python

    >>> from pyteomics import achrom
    >>> achrom.calculate_RT('PEPTIDE', achrom.RCs_guo_ph7_0)
    7.8000000000000025
    
The first argument of the function is the sequence of a peptide in modX 
notation.

The second argument is the set parameters called 'retention coefficients' which
describe chromatographic properties of individual amino acid residues in
a polypeptide chain. :py:mod:`achrom` has a number of predefined sets of 
retention coefficients obtained from publications. The list, detailed 
descriptions and references related to these sets can be found in the module
documentation.

Calibration
...........

The main advantage of the additive model is that it gives more accurate 
predictions if adjusted to specific chromatographic setups and conditions. 
This adjustment, or 'calibration' requires a set of known peptide 
sequences and corresponding retention times (a 'training set'( and returns
a set of new retention coefficients. The following code illustrates the 
calibration procedure in Pyteomics.
    
.. code-block:: python

    >>> from pyteomics import achrom
    >>> RCs = achrom.get_RCs(sequences, RTs)
    >>> achrom.calculate_RT('PEPTIDE', RCs)
    
The first argument of :py:func:`get_RCs` should be a list of modX sequences, 
the second - a list of float-point retention times.
 
As in :py:func:`parser.parse_sequence`:, all non-standard amino modX
acid labels used in the training set should be supplied to `labels` keyword 
argument of :py:func:`get_RCs` along with the standard ones:

.. code-block:: python

    >>> RCs = achrom.get_RCs(sequences, RTs, labels=achrom.std_labels + ['pS', 'pT'])

Advanced calibration
....................

The standard additive model allows a couple of improvements. Firstly, an 
explicit dependency on the length of a peptide may be introduced by multiplying
the retention time by :math:`(1.0 + m * log(L))`, where L is the number of amino
acid residues in a peptide and m is the length correction factor, typically ~ -0.2.

The value of the length correction factor is set at the calibration and stored along
with the retention coefficients. By default, length correction is enabled in
:py:func:`get_RCs` and the factor equals -0.21. You can change
the value of the length correction factor by supplying the 'lcf' keyword argument, 
or you can disable length correction completely by setting lcf=0:

.. code-block:: python

    >>> RCs = achrom.get_RCs(sequences, RTs, lcf=-0.18) # A new value of the length correction factor

    >>> RCs = achrom.get_RCs(sequences, RTs, lcf=0) # Disable length correction.
    
Another considerable improvement over the standard additive model is to treat
terminal amino acid residues as separate chemical entities. This behavior
is disabled by default, but can be enabled by setting term_add=True:

.. code-block:: python

    >>> RCs = achrom.get_RCs(sequences, RTs, term_aa = True) 

This correction is implemented by addition of the 'nterm' and 'cterm' prefixes
to the labels of terminal amino acid residues of the training peptides. In order 
for this correction to work, the training peptides should represent all possible
variations of terminal amino acid residues.
