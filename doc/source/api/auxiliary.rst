auxiliary - common functions and objects
========================================

Math
----

  :py:func:`linear_regression_vertical` - a wrapper for NumPy linear regression,
  minimizes the sum of squares of *y* errors.

  :py:func:`linear_regression` - alias for :py:func:`linear_regression_vertical`.

  :py:func:`linear_regression_perpendicular` - a wrapper for NumPy linear regression,
  minimizes the sum of squares of (perpendicular) distances between the points and the line.


Target-Decoy Approach
---------------------

  :py:func:`qvalues` - estimate q-values for a set of PSMs.

  :py:func:`!filter` - filter PSMs to specified FDR level using TDA or given PEPs.

  :py:func:`filter.chain` - a chained version of :py:func:`!filter`.

  :py:func:`fdr` - estimate FDR in a set of PSMs using TDA or given PEPs.

Project infrastructure
----------------------

  :py:class:`PyteomicsError` - a pyteomics-specific exception.

Helpers
-------

  :py:class:`Charge` - a subclass of :py:class:`int` for charge states.

  :py:class:`ChargeList` - a subclass of :py:class:`list` for lists of charges.

  :py:func:`print_tree` - display the structure of a complex nested
  :py:class:`dict`.

  :py:func:`memoize` - makes a
  `memoization <http://stackoverflow.com/a/1988826/1258041>`_
  `function decorator <http://stackoverflow.com/a/1594484/1258041>`_.

  :py:func:`cvquery` - traverse an arbitrarily nested dictionary looking
  for keys which are :py:class:`cvstr` instances, or objects
  with an attribute called ``accession``.

-------------------------------------------------------------------------------


.. automodule :: pyteomics.auxiliary.math

.. automodule :: pyteomics.auxiliary.target_decoy

.. automodule :: pyteomics.auxiliary.utils

.. automodule :: pyteomics.auxiliary.structures

.. automodule :: pyteomics.auxiliary.file_helpers