Peptide sequence formats
========================

Two common formats for representing peptide sequences are supported in Pyteomics: *modX* and *ProForma*.

*modX* is a custom format historically used in Pyteomics to represent peptide sequences with modifications.
This format is supported by all sequence-related functions in Pyteomics, including functions for digestion, mass and composition calculations,
spectrum annotation, etc. The core functions working with *modX* sequences (sequence parsing, peptidoform generation, etc.) are located in the :py:doc:`api/parser` module.
See :doc:`modX format and the parser module <parser>` for details.

*ProForma* is a more recent format that provides a more structured way to represent peptide sequences and their modifications.
It is supported by most sequence-related functions in Pyteomics, including mass and composition calculations and spectrum annotation.
See :doc:`ProForma format <proforma>` for details.

.. contents::
   :local:
   :depth: 2

.. include:: parser.rst

.. include:: proforma.rst