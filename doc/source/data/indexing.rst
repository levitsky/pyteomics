.. _indexing:

Indexed Parsers
===============

Most of the parsers implement indexing: MGF, mzML, mzXML, FASTA, PEFF, pepXML, mzIdentML, ms1, TraML, featureXML.
Some formats do not have indexing parsers, because there is no unique ID field in the files to identify entries.

XML parser classes are called according to the format, e.g. :py:class:`pyteomics.mzml.MzML`. Text format parsers
that implement indexing are called with the word "Indexed", e.g. :py:class:`pyteomics.fasta.IndexedFASTA`,
as opposed to :py:class:`pyteomics.fasta.FASTA`, which does not implement indexing.
This distinction is due to the fact that indexed parsers need to open the files in binary mode.
This may affect performance for text-based formats and is not always backwards-compatible
(you cannot instantiate an indexed parser class using a previously opened file if it is in text mode).
XML files, on the other hand, are always meant to be opened in binary mode.
So, there is no duplication of classes for XML formats, but indexing can still be disabled by passing
``use_index=False`` to the class constructor or the :py:func:`read` function.

Basic usage
-----------

Indexed parsers can be instantiated using the class name or the :py:func:`read` function::

    In [1]: from pyteomics import mgf

    In [2]: f = mgf.IndexedMGF('tests/test.mgf')

    In [3]: f
    Out[3]: <pyteomics.mgf.IndexedMGF at 0x7fc983cbaeb8>

    In [4]: f.close()

    In [5]: f = mgf.read('tests/test.mgf', use_index=True)

    In [6]: f
    Out[6]: <pyteomics.mgf.IndexedMGF at 0x7fc980c63898>


They support direct assignment and iteration or the `with` syntax, the same way as the older, iterative parsers.

Parser objects can be used as dictionaries mapping entry IDs to entries, or as lists::

    In [7]: f['Spectrum 2']
    Out[7]:
    {'params': {'com': 'Based on http://www.matrixscience.com/help/data_file_help.html',
      'itol': '1',
      'itolu': 'Da',
      'mods': 'Carbamidomethyl (C)',
      'it_mods': 'Oxidation (M)',
      'mass': 'Monoisotopic',
      'username': 'Lou Scene',
      'useremail': 'leu@altered-state.edu',
      'charge': [2, 3],
      'title': 'Spectrum 2',
      'pepmass': (1084.9, 1234.0),
      'scans': '3',
      'rtinseconds': 25.0 second},
     'm/z array': array([ 345.1,  370.2,  460.2, 1673.3, 1674. , 1675.3]),
     'intensity array': array([ 237.,  128.,  108., 1007.,  974.,   79.]),
     'charge array': masked_array(data=[3, 2, 1, 1, 1, 1],
                  mask=False,
            fill_value=0)}

    In [8]: f[1]['params']['title'] # positional indexing
    Out[8]: 'Spectrum 2'

Like dictionaries, indexed parsers support membership testing and :py:func:`len`::

    In [9]: 'Spectrum 1' in f
    Out[9]: True

    In [10]: len(f)
    Out[10]: 2

Rich Indexing
-------------

Indexed parsers also support positional indexing, slices of IDs and integers. ID-based slices include both
endpoints; integer-based slices exclude the right edge of the interval. With integer indexing, *step*
is also supported. Here is a self-explanatory demo of indexing functionality using a test file of two spectra::

    In [11]: len(f['Spectrum 1':'Spectrum 2'])
    Out[11]: 2

    In [12]: len(f['Spectrum 2':'Spectrum 1'])
    Out[12]: 2

    In [13]: len(f[:])
    Out[13]: 2

    In [14]: len(f[:1])
    Out[14]: 1

    In [15]: len(f[1:0])
    Out[15]: 0

    In [16]: len(f[1:0:-1])
    Out[16]: 1

    In [17]: len(f[::2])
    Out[17]: 1

RT-based indexing
.................

In MGF, mzML and mzXML the spectra are usually time-ordered. The corresponding indexed parsers allow accessing the
spectra by retention time, including slices::

    In [18]: f = mzxml.MzXML('tests/test.mzXML')

    In [19]: spec = f.time[5.5] # get the spectrum closest to this retention time

    In [20]: len(f.time[5.5:6.0]) # get spectra from a range
    Out[20]: 2


RT lookup is performed using binary search.
When retrieving ranges, the closest spectra to the start and end of the range
are used as endpoints, so it is possible that they are slightly outside the range.

Multiprocessing
---------------

Indexed parsers provide a unified interface for multiprocessing: :py:meth:`map`.
The method applies a user-defined function to entries from the file, calling it in different processes.
If the function is not provided, the parsing itself is parallelized. Depending on the format,
this may speed up or slow down the parsing overall.
:py:meth:`map` is a generator and yields items as they become available, not preserving the original order::

    In [1]: from pyteomics import mzml

    In [2]: f = mzml.MzML('tests/test.mzML')

    In [3]: for spec in f.map():
       ...:     print(spec['id'])
       ...:
    controllerType=0 controllerNumber=1 scan=2
    controllerType=0 controllerNumber=1 scan=1

    In [4]: for item in f.map(lambda spec: spec['id']):
       ...:     print(item)
       ...:
    controllerType=0 controllerNumber=1 scan=1
    controllerType=0 controllerNumber=1 scan=2


.. note ::
  To use :py:meth:`map` with lambda functions (and in some other corner cases, like
  parsers instantiated with pre-opened file objects), the :py:mod:`dill` package is required.
  This is because the target callable and the parser itself need to be pickled for multiprocessing to work.

Apart from parser objects, :py:meth:`map` is available on objects returned by :py:func:`chain` functions
and :py:meth:`iterfind`::

    In [5]: for c in f.iterfind('chromatogram').map():
       ...:     print(c['id'])
       ...:
    TIC

    In [6]: for spec in mzml.chain('tests/test.mzML', 'tests/test.mzML').map():
       ...:     print(spec['id'])
       ...:
    controllerType=0 controllerNumber=1 scan=1
    controllerType=0 controllerNumber=1 scan=2
    controllerType=0 controllerNumber=1 scan=1
    controllerType=0 controllerNumber=1 scan=2