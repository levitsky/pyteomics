mzML
====

**mzML** is an XML-based format for experimental data obtained on MS/MS or LC-MS
setups. **Pyteomics** offers you the functionality of :py:mod:`mzml.py` module
to gain access to the information contained in .mzML files from Python.

Examples
--------

The main (and the top-most-level) function in this module is
:py:func:`iter_spectrum`::

    >>> from pyteomics.mzml import iter_spectrum
    >>> s = iter_spectrum('my.mzml')
    >>> for spectrum in s:
        …

.. note:: :py:func:`iter_spectrum` returns a *generator object* instead of a 
    *list* to prevent unnecessary memory use.
    Type::

        spectrum_list = [x for x in iter_spectrum(‘my.mzml’)]
       
    to get a list.
