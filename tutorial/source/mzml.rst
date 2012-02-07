mzML
====

**mzML** is an XML-based format for experimental data obtained on MS/MS or LC-MS
setups. **Pyteomics** offers you the functionality of :py:mod:`mzml.py` module
to gain access to the information contained in .mzML files from Python.

Examples
--------

The main (and the top-most-level) function in this module is
:py:func:`iter_spectrum`. It allows the user to iterate through spectra
contained in the file. In the following example a list of parent m/z is created::

    >>> from pyteomics.mzml import iter_spectrum
    >>> s = iter_spectrum('my.mzml')
    >>> parent_mz = [spectrum['base peak m/z'] for spectrum in s]

The list can now be used for all kinds of statistical analysis. For example, one
can calculate the corresponding histogram with pylab::

    >>> import pylab
    >>> pylab.figure()
    >>> pylab.hist(parent_mz, bins=100)
    >>> pylab.show()

The data accessible through :py:func:`iter_spectrum` represent the data
contained in the corresponding mzML file and can thus differ from one file
to another. However, one can easily inspect the format just by printing any
of the entries::

    >>> print iter_spectrum('my.mzml').next()

Then, you can develop more complex scripts relying on the exact format of the 
files you work with.

.. note::

    :py:func:`iter_spectrum` returns a generator object instead of a 
    list to prevent unnecessary memory use.

If you actually **need** a list, type::

    spectrum_list = [x for x in iter_spectrum(‘my.mzml’)]
      

