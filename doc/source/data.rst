===========
Data Access
===========

The following section is dedicated to data manipulation. **Pyteomics** aims to
support the most common formats of (LC-)MS/MS data, peptide identification
results and protein databases.

.. include :: data/notes.rst

.. contents:: Document contents
    :backlinks: top

.. include :: data/text.rst

.. include :: data/xml.rst


FDR estimation and filtering
============================

Three modules for reading proteomics search engine output (:py:mod:`tandem`,
:py:mod:`pepxml` and :py:mod:`mzid`) expose similar functions
:py:func:`is_decoy`, :py:func:`fdr` and :py:func:`!filter`. These functions
implement the widely used
Target-Decoy Approach (TDA) to estimation of False Discovery Rate (FDR).

The :py:func:`is_decoy` function is supposed to determine if a particular
spectrum identification is coming from the decoy database. In :py:mod:`tandem`
and :py:mod:`pepxml` this is done by checking if the protein description/name
starts with a certain prefix. In :py:mod:`mzid`, a boolean value that stores
this information in the PSM dict is used.

.. warning ::
     Because of the variety of the software producing files in pepXML and
     mzIdentML formats, the :py:func:`is_decoy` function provided in the
     corresponding modules may not work for your specific files. In this case
     you will have to refer to the source of
     :py:func:`pyteomics.pepxml.is_decoy` and
     :py:func:`pyteomics.mzid.is_decoy` and create your own function in a
     similar manner.

The :py:func:`fdr` function estimates the FDR in a set of PSMs by counting
the decoy matches. Since it is using the :py:func:`is_decoy` function, the
warning above applies. You can supply a custom function so that :py:func:`fdr`
works for your data. :py:func:`fdr` can also be imported from
:py:mod:`auxiliary`, where it has no default for :py:func:`is_decoy`.

The :py:func:`!filter` function works like :py:func:`chain`, but instead of
yielding all PSMs, it filters them to a certain level of FDR. PSM filtering
requires counting decoy matches, too (see above), but it also implies sorting
the PSMs by some kind of a score. This score cannot be universal due to the
above-mentioned reasons, and it can be specified as a user-defined function.
For instance, the default sorting key in :py:func:`pyteomics.mzid.filter` is
only expected to work with mzIdentML files created with Mascot.
So once again,

.. warning ::
     The default parameters of :py:func:`!filter` may not work for your files.

There are also :py:func:`filter.chain` and
:py:func:`filter.chain.from_iterable`. These are different from
:py:func:`!filter` in that they apply FDR filtering to all files separately
and then provide a reader over top PSMs of all files, whereas
:py:func:`!filter` pools all PSMs together and applies a single threshold.

If you want to filter a list representing PSMs in arbitrary format, you can
use :py:func:`pyteomics.auxiliary.filter`. Instead of files it takes lists
(or other iterables) of PSMs. The rest is the same as for other
:py:func:`!filter` functions.

NumPy and Pandas support, etc.
------------------------------

:py:func:`pyteomics.auxiliary.filter` supports structured :py:mod:`numpy` arrays and
:py:class:`pandas.DataFrames` of PSMs. This makes it easy to filter search results
stored as CSV files (see :ref:`example-3` for more info).

Generally, PSMs can be provided as iterators, lists, arrays, and :py:class:`DataFrames`,
and `key` and `is_decoy` parameters to :py:func:`!filter` can be functions, strings,
lists, arrays, or iterators. If a string is given, it is used as a key in a structured
array, :py:class:`DataFrame` or an iterable of :py:class:`dicts`.