Charge and pI
=============

Electrochemical properties of polypeptides can be assessed via the
:py:mod:`electrochem` module. For now, it allows to calculate:

*  the charge of a polypeptide molecule at given pH;
*  the isoelectric point.

The :py:mod:`electrochem` module is based on the Henderson-Hasselbalch 
equation.


Examples
--------

Both functions in the module accepts input in form of a modX sequence, 
a parsed sequence or a dict with amino acid composition.

.. code-block:: python

    >>> from pyteomics import electrochem
    >>> electrochem.charge('PEPTIDE', 7)
    -2.9980189709606284
    >>> from pyteomics import parser
    >>> parsed_seq = parser.parse_sequence('PEPTIDE', show_unmodified_termini=True)
    >>> electrochem.charge(parsed_seq, 7)
    -2.9980189709606284
    >>> aa_composition = parser.amino_acid_composition('PEPTIDE', show_unmodified_termini=True)
    >>> electrochem.charge(aa_composition, 7)
    -2.9980189709606284
    >>> electrochem.pI('PEPTIDE')
    2.87451171875
    >>> electrochem.pI('PEPTIDE', precision_pI=0.0001)
    2.876354217529297

.. plot:: source/charge_vs_ph.py

Customization
-------------

The pKas of individual amino acids are stored in dicts in the following format:
{`modX label` : (`pKa`, `charge`)}. The module contains several datasets 
published in scientific journals: :py:data:`pK_lehninger` (used by default), 
:py:data:`pK_sillero`, :py:data:`pK_dawson`, :py:data:`pK_rodwell`. 

