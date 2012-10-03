Example 1: in silico cleavage, property calculation and visualization
=====================================================================

In this example, a FASTA file with yeast protein sequences is retrieved,
opened, and used to generate peptides *in silico*. For the peptide sequences,
probable charge states, as well as masses, m/z and retention times for normal
and reversed phase are calculated (additive model is used for retention time
prediction).
The code is shown in two versions, for Python 2.7 and Python 3.x.
(Note: currently, there's no official matplotlib release for Python 3.)

Python 2.7
----------

.. literalinclude:: fasta.py2.7.py
    :linenos:
    :language: python

Python 3.x
----------

.. literalinclude:: fasta.py3.x.py
    :linenos:
    :language: python

