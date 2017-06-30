===========
Data Access
===========

The following section is dedicated to data manipulation. **Pyteomics** aims to
support the most common formats of (LC-)MS/MS data, peptide identification
results and protein databases.

General Notes
=============

- Each module mentioned below corresponds to a file format. In each module, the
  top-level function :py:func:`read` allows iteration over entries in a file.
  It works like the built-in :py:func:`open`, allowing direct iteration and
  supporting the ``with`` syntax, which we recommend using. So you can do:

  .. code-block :: python

       >>> from pyteomics import mgf
       >>> reader = mgf.read('tests/test.mgf')
       >>> for spectrum in reader:
       >>>    ...
       >>> reader.close()

  ... but it is recommended to do:

  .. code-block :: python

       >>> from pyteomics import mgf
       >>> with mgf.read('tests/test.mgf') as reader:
       >>>     for spectrum in reader:
       >>>        ...

- Apart from :py:func:`read`, which reads just one file, all modules described
  here have functions for reading multiple files: :py:func:`chain` and
  :py:func:`chain.from_iterable`.
  ``chain('f1', 'f2')`` is equivalent to ``chain.from_iterable(['f1', 'f2'])``.
  :py:func:`chain` and :py:func:`chain.from_iterable` only support the
  ``with`` syntax. If you don't want to use the ``with`` syntax, you can just
  use the :py:mod:`itertools` functions :py:func:`chain` and
  :py:func:`chain.from_iterable`.

- Throughout this section we use
  :py:func:`pyteomics.auxiliary.print_tree` to display the structure of the
  data returned by various parsers.

.. contents:: Document contents
    :backlinks: top

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
as :py:class:`numpy.ndarray` objects `mass array` and `intensity array`.
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
with very similar signatures for reading headers and spectra from files.

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
:py:func:`pyteomics.fasta.decoy_sequence` function. It supports two modes:
``'reverse'`` and ``'random'``.

.. code-block:: python

    >>> fasta.decoy_sequence('PEPTIDE', 'reverse')
    'EDITPEP'
    >>> fasta.decoy_sequence('PEPTIDE', 'random')
    'TPPIDEE'
    >>> fasta.decoy_sequence('PEPTIDE', 'random')
    'PTIDEPE'


XML formats
===========

XML parsers are implemented as classes and provide an
object-oriented interface. The functional interface is preserved for backward
compatibility and wraps the actual class-based machinery.
That means that reader objects returned
by :py:func:`read` functions have additional methods.

One of the most important methods is :py:meth:`iterfind`. It allows reading
additional information from XML files.

mzML and mzXML
--------------

**mzML** and **mzXML** are XML-based formats for experimental data obtained on MS/MS or LC-MS
setups. **Pyteomics** offers you the functionality of :py:mod:`pyteomics.mzml` and
:py:mod:`pyteomics.mzxml`
modules to gain access to the information contained in those files from Python.
The interfaces of the two modules are very similar, this section will use **mzML**
for demonstration.

The user can iterate through MS/MS spectra contained in a file via the
:py:func:`pyteomics.mzml.read` function or :py:class:`pyteomics.mzml.MzML` class.
Here is an example of the output:

.. code-block:: python

    >>> from pyteomics import mzml, auxiliary
    >>> with mzml.read('tests/test.mzML') as reader:
    >>>     auxiliary.print_tree(next(reader))
    count
    index
    highest observed m/z
    ms level
    total ion current
    intensity array
    lowest observed m/z
    defaultArrayLength
    profile spectrum
    MSn spectrum
    positive scan
    base peak intensity
    m/z array
    base peak m/z
    id
    scanList
     -> count
     -> scan [list]
     ->  -> scan start time
     ->  -> preset scan configuration
     ->  -> filter string
     ->  -> instrumentConfigurationRef
     ->  -> scanWindowList
     ->  ->  -> count
     ->  ->  -> scanWindow [list]
     ->  ->  ->  -> scan window lower limit
     ->  ->  ->  -> scan window upper limit
     ->  -> [Thermo Trailer Extra]Monoisotopic M/Z:
     -> no combination

Additionally, :py:class:`pyteomics.mzml.MzML` objects support direct indexing
with spectrum IDs:

.. code-block:: python

    >>> from pyteomics import mzml
    >>> with mzml.MzML('tests/test.mzML') as reader:
    >>>    spectrum = reader["controllerType=0 controllerNumber=1 scan=1"]
    >>>    mz = spectrum['m/z array'] # do something with "spectrum"
    >>>    ...

:py:class:`pyteomics.mzml.PreIndexedMzML` offers the same functionality,
but it uses byte offset information found at the end of the file.
Unlike the rest of the functions and classes, :py:class:`pyteomics.mzml.PreIndexedMzML`
does not have a counterpart in :py:mod:`pyteomics.mzxml`.

pepXML
------

`pepXML <http://tools.proteomecenter.org/wiki/index.php?title=Formats:pepXML>`_
is a widely used XML-based format for peptide identifications.
It contains information about the MS data, the parameters of the search engine
used and the assigned sequences. To access these data, use
:py:mod:`pyteomics.pepxml` module.

The function :py:func:`pyteomics.pepxml.read` iterates through Peptide-Spectrum
matches in a pepXML file and returns them as a custom dict. Alternatively, you
can use the :py:class:`pyteomics.pepxml.PepXML` interface.

.. code-block:: python

    >>> from pyteomics import pepxml, auxiliary
    >>> with pepxml.read('tests/test.pep.xml') as reader:
    >>>     auxiliary.print_tree(next(reader))
    end_scan
    search_hit [list]
     -> hit_rank
     -> calc_neutral_pep_mass
     -> modifications
     -> modified_peptide
     -> peptide
     -> num_matched_ions
     -> search_score
     ->  -> deltacn
     ->  -> spscore
     ->  -> sprank
     ->  -> deltacnstar
     ->  -> xcorr
     -> num_missed_cleavages
     -> analysis_result [list]
     ->  -> peptideprophet_result
     ->  ->  -> all_ntt_prob
     ->  ->  -> parameter
     ->  ->  ->  -> massd
     ->  ->  ->  -> fval
     ->  ->  ->  -> nmc
     ->  ->  ->  -> ntt
     ->  ->  -> probability
     ->  -> analysis
     -> tot_num_ions
     -> num_tot_proteins
     -> is_rejected
     -> proteins [list]
     ->  -> num_tol_term
     ->  -> protein
     ->  -> peptide_next_aa
     ->  -> protein_descr
     ->  -> peptide_prev_aa
     -> massdiff
    index
    assumed_charge
    spectrum
    precursor_neutral_mass
    start_scan

Reading into a pandas.DataFrame
...............................

If you like working with tabular data using :py:mod:`pandas`, you can load data from pepXML files
directly into :py:class:`pandas.DataFrames`
using the :py:func:`pyteomics.pepxml.DataFrame` function. It can read multiple files
at once (using :py:func:`pyteomics.pepxml.chain`) and return a combined table with
essential information about search results. This function requires :py:mod:`pandas`.

X!Tandem
--------

`X!Tandem search engine <http://www.thegpm.org/tandem/>`_ has its own output
format that contains more info than pepXML. **Pyteomics** has a reader for it
in the :py:mod:`pyteomics.tandem` module.

.. code-block:: python

    >>> from pyteomics import tandem, auxiliary
    >>> with tandem.read('tests/test.t.xml') as reader:
    ...     auxiliary.print_tree(next(reader))
    ...
    rt
    support
     -> fragment ion mass spectrum
     ->  -> M+H
     ->  -> note
     ->  -> charge
     ->  -> Ydata
     ->  ->  -> units
     ->  ->  -> values
     ->  -> Xdata
     ->  ->  -> units
     ->  ->  -> values
     ->  -> label
     ->  -> id
     -> supporting data
     ->  -> convolution survival function
     ->  ->  -> Ydata
     ->  ->  ->  -> units
     ->  ->  ->  -> values
     ->  ->  -> Xdata
     ->  ->  ->  -> units
     ->  ->  ->  -> values
     ->  ->  -> label
     ->  -> b ion histogram
     ->  ->  -> Ydata
     ->  ->  ->  -> units
     ->  ->  ->  -> values
     ->  ->  -> Xdata
     ->  ->  ->  -> units
     ->  ->  ->  -> values
     ->  ->  -> label
     ->  -> y ion histogram
     ->  ->  -> Ydata
     ->  ->  ->  -> units
     ->  ->  ->  -> values
     ->  ->  -> Xdata
     ->  ->  ->  -> units
     ->  ->  ->  -> values
     ->  ->  -> label
     ->  -> hyperscore expectation function
     ->  ->  -> a1
     ->  ->  -> a0
     ->  ->  -> Ydata
     ->  ->  ->  -> units
     ->  ->  ->  -> values
     ->  ->  -> Xdata
     ->  ->  ->  -> units
     ->  ->  ->  -> values
     ->  ->  -> label
    mh
    maxI
    expect
    sumI
    act
    fI
    z
    id
    protein [list]
     -> peptide
     ->  -> pre
     ->  -> end
     ->  -> seq
     ->  -> b_ions
     ->  -> nextscore
     ->  -> mh
     ->  -> y_ions
     ->  -> start
     ->  -> hyperscore
     ->  -> expect
     ->  -> delta
     ->  -> id
     ->  -> post
     ->  -> missed_cleavages
     ->  -> b_score
     ->  -> y_score
     -> uid
     -> sumI
     -> label
     -> note
     -> expect
     -> file
     ->  -> URL
     ->  -> type
     -> id

:py:func:`pyteomics.tandem.read` returns a
:py:class:`pyteomics.tandem.TandemXML` instance, which can also be
created directly.

Reading into a pandas.DataFrame
...............................

You can also load data from X!Tandem files directly into :py:class:`pandas.DataFrames`
using the :py:func:`pyteomics.tandem.DataFrame` function. It can read multiple files
at once (using :py:func:`pyteomics.tandem.chain`) and return a combined table with
essential information about search results. Of course, this function requires :py:mod:`pandas`.

mzIdentML
---------

`mzIdentML <http://www.psidev.info/mzidentml>`_  is one of the standards
developed by the Proteomics Informatics working group of the HUPO Proteomics
Standard Initiative.

The module interface is similar to that of the other reader modules.
The :py:func:`pyteomics.mzid.read` function returns a
:py:class:`pyteomics.mzid.MzIdentML` instance, which you can just as easily
use directly.

.. code-block:: python

    >>> from pyteomics import mzid, auxiliary
    >>> with mzid.read('tests/test.mzid') as reader:
    >>>     auxiliary.print_tree(next(reader))
    SpectrumIdentificationItem [list]
     -> PeptideEvidenceRef [list]
     ->  -> peptideEvidence_ref
     -> ProteinScape:SequestMetaScore
     -> chargeState
     -> rank
     -> ProteinScape:IntensityCoverage
     -> calculatedMassToCharge
     -> peptide_ref
     -> passThreshold
     -> experimentalMassToCharge
     -> id
    spectrumID
    id
    spectraData_ref


Element IDs and references
..........................

In *mzIdentML*, some elements contain references to other elements in the same
file. The references are simply XML attributes whose name ends with ``_ref`` and
the value is an ID, identical to the value of the ``id`` attribute of a certain
element.

The parser can retrieve information from these references on the fly, which can
be enabled by passing ``retrieve_refs=True`` to the
:py:meth:`pyteomics.mzid.MzIdentML.iterfind` method, to
:py:class:`pyteomics.mzid.MzIdentML` constructor, or to
:py:func:`pyteomics.mzid.read`. Retrieval of data by ID is implemented in
the :py:meth:`pyteomics.mzid.MzIdentML.get_by_id` method. Alternatively, the
:py:class:`MzIdentML` object itself can be indexed with element IDs::

    >>> from pyteomics import mzid
    >>> m = mzid.MzIdentML('tests/test.mzid')
    >>> m['ipi.HUMAN_decoy']
    {'DatabaseName': 'database IPI_human',
     'decoy DB accession regexp': '^SHD',
     'decoy DB generation algorithm': 'PeakQuant.DecoyDatabaseBuilder',
     'id': 'ipi.HUMAN_decoy',
     'location': 'file://www.medizinisches-proteom-center.de/DBServer/ipi.HUMAN/3.15/ipi.HUMAN_decoy.fasta',
     'name': ['decoy DB from IPI_human',
      'DB composition target+decoy',
      'decoy DB type shuffle'],
     'numDatabaseSequences': 58099,
     'releaseDate': '2006-02-22T09:30:47Z',
     'version': '3.15'}
    >>> m.close()


.. note:: Since version 3.3, :py:class:`pyteomics.mzid.MzIdentML` objects keep an index of byte
          offsets for some of the elements. It helps achieve acceptable performance
          when using ``retrieve_refs=True``, or when accessing individual elements by their ID.

          This behavior can be disabled by passing
          ``use_index=False`` to the object constructor.
          An alternative, older mechanism is caching of element IDs. To build
          a cache for a file, you can pass ``build_id_cache=True`` and ``use_index=False``
          to the :py:class:`MzIdentML` constructor, or to :py:func:`pyteomics.mzid.read`,
          or call the :py:meth:`pyteomics.mzid.MzIdentML.build_id_cache` method
          prior to reading the data.


FeatureXML
----------

:py:mod:`pyteomics.openms.featurexml` implements a simple parser for **.featureXML** files
used in the `OpenMS <http://open-ms.sourceforge.net/about/>`_ framework. The usage
is identical to other XML parsing modules. Since **featureXML** has feature IDs,
:py:class:`FeatureXML` objects also support direct indexing as well as iteration::

    >>> from pyteomics.openms import featurexml

    >>> # function style, iteration
    ... with featurexml.read('tests/test.featureXML') as f:
    ...     qual = [feat['overallquality'] for feat in f]
    ...

    >>> qual # qualities of the two features in the test file
    [0.791454, 0.945634]

    >>> # object-oriented style, direct indexing
    >>> f = featurexml.FeatureXML('tests/test.featureXML')
    >>> f['f_189396504510444007']['overallquality']
    0.945634
    >>> f.close()

As always, :py:func:`pyteomics.openms.featurexml.read`
and :py:class:`pyteomics.openms.featurexml.FeatureXML` are interchangeable.

TrafoXML
--------

**.trafoXML** is another OpenMS format based on XML. It describes a
tranformation produced by an RT alignment algorithm. The file basically contains a series
of `(from; to)` pairs corresponding to original and transformed retention times:

   >>> from pyteomics.openms import trafoxml
   >>> from_rt, to_rt = [], []
   >>> with trafoxml.read('test/test.trafoXML') as f:
   ...    for pair in f:
   ...        from_rt.append(pair['from'])
   ...        to_rt.append(pair['to'])

   >>> # plot the transformation
   >>> import pylab
   >>> pylab.plot(from_rt, to_rt)

As always, :py:func:`pyteomics.openms.trafoxml.read`
and :py:class:`pyteomics.openms.trafoxml.TrafoXML` are interchangeable.

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