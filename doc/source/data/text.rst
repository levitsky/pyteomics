Text-based formats
==================

MGF
---

Mascot Generic Format
(`MGF <http://www.matrixscience.com/help/data_file_help.html>`_) is a simple
human-readable format for MS/MS data. It allows storing MS/MS peak lists and
exprimental parameters. :py:mod:`pyteomics.mgf` is a module that implements
reading and writing MGF files.

Reading
.......

:py:func:`pyteomics.mgf.read` function allows iterating through spectrum entries.
Spectra are represented as :py:class:`dicts`. By default, MS/MS peak lists are stored
as :py:class:`numpy.ndarray` objects `m/z array` and `intensity array`.
Fragment charges will be stored in a masked array under the `charge array` key.
Parameters are stored as a :py:class:`dict` under `params` key.

Here is an example of use::

    >>> from pyteomics import mgf, auxiliary
    >>> with mgf.read('tests/test.mgf') as reader:
    >>>     auxiliary.print_tree(next(reader))
    m/z array
    params
     -> username
     -> useremail
     -> mods
     -> pepmass
     -> title
     -> itol
     -> charge
     -> mass
     -> itolu
     -> it_mods
     -> com
    intensity array
    charge array

To speed up parsing, or if you want to avoid using :py:mod:`numpy`, you can tweak the
behaviour of :py:func:`pyteomics.mgf.read` with parameters `convert_arrays` and `read_charges`.


Reading file headers
....................

Also, :py:mod:`pyteomics.mgf` allows to extract headers with general
parameters from MGF files with :py:func:`pyteomics.mgf.read_header` function. It
also returns a :py:class:`dict`.

.. code-block:: python

    >>> header = mgf.read_header('tests/test.mgf')
    >>> auxiliary.print_tree(header)
    itolu
    itol
    username
    com
    useremail
    it_mods
    charge
    mods
    mass

Class-based interface
.....................

Since version 3.4.3, MGF parsing functionality is encapsulated in a class:
:py:class:`pyteomics.mgf.MGF`. This class can be used for:

 - sequential parsing of the file (the same as :py:func:`read`)::

 .. code-block:: python

    >>> with mgf.MGF('tests/test.mgf') as reader:
    ..:     for spectrum in reader:
    ..:         ...

 - accessing the file header (the same as :py:func:`read_header`)::

 .. code-block:: python

    >>> f = mgf.MGF('tests/test.mgf')
    >>> f.header
    {'charge': [2, 3],
     'com': 'Based on http://www.matrixscience.com/help/data_file_help.html',
     'it_mods': 'Oxidation (M)',
     'itol': '1',
     'itolu': 'Da',
     'mass': 'Monoisotopic',
     'mods': 'Carbamidomethyl (C)',
     'useremail': 'leu@altered-state.edu',
     'username': 'Lou Scene'}

 - direct access to spectra by title (the same as :py:func:`get_spectrum`)::

 .. code-block:: python

    >>> f = mgf.MGF('tests/test.mgf')
    >>> f['Spectrum 2']
    {'charge array': masked_array(data = [3 2 1 1 1 1],
                  mask = False,
            fill_value = 0),
     'intensity array': array([  237.,   128.,   108.,  1007.,   974.,    79.]),
     'm/z array': array([  345.1,   370.2,   460.2,  1673.3,  1674. ,  1675.3]),
     'params': {'charge': [2, 3],
      'com': 'Based on http://www.matrixscience.com/help/data_file_help.html',
      'it_mods': 'Oxidation (M)',
      'itol': '1',
      'itolu': 'Da',
      'mass': 'Monoisotopic',
      'mods': 'Carbamidomethyl (C)',
      'pepmass': (1084.9, 1234.0),
      'rtinseconds': '25',
      'scans': '3',
      'title': 'Spectrum 2',
      'useremail': 'leu@altered-state.edu',
      'username': 'Lou Scene'}}

.. note ::
    :py:class:`MGF`'s support for direct indexing is rudimentary, because it does not in fact keep an index and has
    to search through the file line-wise on every call. :py:class:`pyteomics.mgf.IndexedMGF` iis designed for
    random access and more (see `Indexed Parsers`_ for details).

Writing
.......

Creation of MGF files is implemented in :py:func:`pyteomics.mgf.write` function.
The user can specify the header, an iterable of spectra in the same format as
returned by :py:func:`read`, and the output path.

.. code-block:: python

    >>> spectra = mgf.read('tests/test.mgf')
    >>> mgf.write(spectra=spectra, header=header)
    USERNAME=Lou Scene
    ITOL=1
    USEREMAIL=leu@altered-state.edu
    MODS=Carbamidomethyl (C)
    IT_MODS=Oxidation (M)
    CHARGE=2+ and 3+
    MASS=Monoisotopic
    ITOLU=Da
    COM=Taken from http://www.matrixscience.com/help/data_file_help.html

    BEGIN IONS
    TITLE=Spectrum 1
    PEPMASS=983.6
    846.6 73.0
    846.8 44.0
    847.6 67.0
    1640.1 291.0
    1640.6 54.0
    1895.5 49.0
    END IONS

    BEGIN IONS
    TITLE=Spectrum 2
    RTINSECONDS=25
    PEPMASS=1084.9
    SCANS=3
    345.1 237.0
    370.2 128.0
    460.2 108.0
    1673.3 1007.0
    1674.0 974.0
    1675.3 79.0
    END IONS

MS1
---

`MS1 <http://dx.doi.org/10.1002/rcm.1603>`_ is a simple
human-readable format for MS1 data. It allows storing MS1 peak lists and
exprimental parameters. Just like MS1 format is quite similar to MGF,
the corresponding module (:py:mod:`pyteomics.ms1`) provides the same functions
and classes with very similar signatures for reading headers and spectra from files.

Writing is not supported at this time.

FASTA
-----

FASTA is a common format for protein sequence databases.

Reading
.......

To extract data from FASTA databases, use the :py:func:`pyteomics.fasta.read`
function.

.. code-block:: python

    >>> from pyteomics import fasta
    >>> with fasta.read('/path/to/file/my.fasta') as db:
    >>>     for entry in db:
    >>>         ...

Just like other parsers in **Pyteomics**, :py:func:`pyteomics.fasta.read`
returns a *generator object* instead of a
:py:class:`list` to prevent excessive memory use. The generator yields
(description, sequence) tuples, so it's natural to use it as follows:

.. code-block:: python

    >>> with fasta.read('/path/to/file/my.fasta') as db:
    >>>     for descr, seq in db:
    >>>         ...

You can also use attributes to access description and sequence:

.. code-block:: python

    >>> with fasta.read('my.fasta') as reader:
    >>>     descriptions = [item.description for item in reader]

Description parsing
...................

You can specify a function that will be applied to the FASTA headers for
your convenience. :py:data:`pyteomics.fasta.std_parsers` has some pre-defined
parsers that can be used for this purpose.

.. code-block:: python

    >>> with fasta.read('HUMAN.fasta', parser=fasta.std_parsers['uniprotkb']) as r:
    >>>    print(next(r).description)
    {'PE': 2, 'gene_id': 'LCE6A', 'GN': 'LCE6A', 'id': 'A0A183', 'taxon': 'HUMAN',
     'SV': 1, 'OS': 'Homo sapiens', 'entry': 'LCE6A_HUMAN',
     'name': 'Late cornified envelope protein 6A', 'db': 'sp'}

or try guessing the header format:

.. code-block:: python

    >>> with fasta.read('HUMAN.fasta', parser=fasta.parse) as r:
    >>>    print(next(r).description)
    {'PE': 2, 'gene_id': 'LCE6A', 'GN': 'LCE6A', 'id': 'A0A183', 'taxon': 'HUMAN',
     'SV': 1, 'OS': 'Homo sapiens', 'entry': 'LCE6A_HUMAN',
     'name': 'Late cornified envelope protein 6A', 'db': 'sp'}

Class-based interface
.....................

The :py:class:`pyteomics.fasta.FASTA` class is available for text-based (old style) parsing
(the same as shown with :py:func:`read` above). Also, the new binary-mode, indexed parser,
:py:class:`pyteomics.fasta.IndexedFASTA` implements all the perks of the `Indexed Parsers`_.
Both classes also have a number of flavor-specific subclasses that implement header parsing.

Additionally, flavored indexed parsers allow accessing the protein entries by the extracted ID field,
while the regular :py:class:`pyteomics.fasta.IndexedFASTA` uses full description string for identification::

    In [1]: from pyteomics import fasta

    In [2]: db = fasta.IndexedUniProt('sprot_human.fasta') # A SwissProt database

    In [3]: len(db['Q8IYH5'].sequence)
    Out[3]: 903

    In [4]: db['Q8IYH5'] == db['sp|Q8IYH5|ZZZ3_HUMAN ZZ-type zinc finger-containing protein 3 OS=Homo sapiens GN=ZZZ3 PE=1 SV=1']
    Out[4]: True


Writing
.......

You can also create a FASTA file using a sequence of `(description, sequence)`
:py:class:`tuples`.

.. code-block:: python

    >>> entries = [('Protein 1', 'PEPTIDE'*1000), ('Protein 2', 'PEPTIDE'*2000)]
    >>> fasta.write(entries, 'target-file.fasta')

Decoy databases
...............

Another common task is to generate a *decoy database*. **Pyteomics** allows
that by means of the :py:func:`pyteomics.fasta.decoy_db` and
:py:func:`pyteomics.fasta.write_decoy_db` functions.

.. code-block:: python

    >>> fasta.write_decoy_db('mydb.fasta', 'mydb-with-decoy.fasta')

The only required argument is the first one, indicating the source database. The
second argument is the target file and defaults to system standard output.

If you need to modify a single sequence, use the
:py:func:`pyteomics.fasta.decoy_sequence` function. It supports three modes:
``'reverse'``, ``'shuffle'``, and ``'fused'`` (see :py:func:`pyteomics.fasta.reverse`,
:py:func:`pyteomics.fasta.shuffle` and :py:func:`pyteomics.fasta.fused_decoy` for documentation).

.. code-block:: python

    >>> fasta.decoy_sequence('PEPTIDE', 'reverse')
    'EDITPEP'
    >>> fasta.decoy_sequence('PEPTIDE', 'shuffle')
    'TPPIDEE'
    >>> fasta.decoy_sequence('PEPTIDE', 'shuffle')
    'PTIDEPE'


mzTab
-----

mzTab is a HUPO-PSI standardized text-based format for describing identification
and quantification of peptides and small molecules. You can read an mzTab file into
a set of :py:class:`pandas.DataFrame` objects with the :py:class:`pyteomics.mztab.MzTab`
class.

.. code-block:: python

  >>> from pyteomics import mztab
  >>> tables = mztab.MzTab("path/to/file.mzTab")
  >>> psms = tables.spectrum_match_table
  >>> # do something with DataFrame
