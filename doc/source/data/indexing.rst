Indexed Parsers
===============

Most of the parsers implement indexing: MGF, mzML, mzXML, FASTA, PEFF, pepXML, mzIdentML, ms1, TraML, featureXML.
Some formats do not have indexing parsers, because there is no unique ID field in the files to identify entries.

XML parser classes are called according to the format, e.g. :py:class:`pyteomics.mzml.MzML`. Text format parsers
that implement indexing are called with the word "Indexed", e.g. `:py:class:`pyteomics.fasta.IndexedFASTA`,
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

Parser objects can be used as dictionaries mapping entry IDs to entries::

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

Rich Indexing
-------------

Indexed parsers also support positional indexing, slices of IDs and integers. ID-based slices include both
endpoints; integer-based slices exclude the right edge of the interval. With integer indexing, __step__
is also supported. Here is a self-explanatory demo of indexing functionality using a test file of two spectra::

    In [9]: len(f['Spectrum 1':'Spectrum 2'])
    Out[9]: 2

    In [10]: len(f['Spectrum 2':'Spectrum 1'])
    Out[10]: 2

    In [11]: len(f[:])
    Out[11]: 2

    In [12]: len(f[:1])
    Out[12]: 1

    In [13]: len(f[1:0])
    Out[13]: 0

    In [14]: len(f[1:0:-1])
    Out[14]: 1

    In [15]: len(f[::2])
    Out[15]: 1

RT-based indexing
.................

In MGF, mzML and mzXML the spectra are usually time-ordered. The corresponding indexed parsers allow accessing the
spectra by retention time, including slices::

Multiprocessing
---------------