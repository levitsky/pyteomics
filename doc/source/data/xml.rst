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
with spectrum IDs and all other features of `Indexed Parsers`_.

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

If you like working with tabular data using :py:mod:`pandas`, you can load pepXML files
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
          offsets for some of the elements (see `Indexed Parsers`_).
          Indexing helps achieve acceptable performance
          when using ``retrieve_refs=True``, or when accessing individual elements by their ID.

          This behavior can be disabled by passing
          ``use_index=False`` to the object constructor.
          An alternative, older mechanism is caching of element IDs. To build
          a cache for a file, you can pass ``build_id_cache=True`` and ``use_index=False``
          to the :py:class:`MzIdentML` constructor, or to :py:func:`pyteomics.mzid.read`,
          or call the :py:meth:`pyteomics.mzid.MzIdentML.build_id_cache` method
          prior to reading the data.

Reading into a pandas.DataFrame
...............................

:py:mod:`pyteomics.mzid` also provides a :py:func:`pyteomics.mzid.DataFrame` function
that reads one or several files into a single Pandas :py:class:`DataFrame`.
This function requires :py:mod:`pandas`.

TraML
-----

`TraML <http://www.psidev.info/traml>`_ is also a PSI format. It stores a lot of information on SRM experiments.
The parser, :py:class:`pyteomics.traml.TraML`, iterates over `<Transition>` elements by default.
Like `MzIdentML`_, it has a `retrieve_refs` parameter that helps pull in the information from other parts of the file.
:py:class:`TraML` is one of the `Indexed Parsers`_.

FeatureXML
----------

:py:mod:`pyteomics.openms.featurexml` implements a simple parser for **.featureXML** files
used in the `OpenMS <http://open-ms.sourceforge.net/about/>`_ framework. The usage
is identical to other XML parsing modules. Since **featureXML** has feature IDs,
:py:class:`FeatureXML` objects also support direct indexing as well as iteration, among
the many features of `Indexed Parsers`_::

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
of `(from; to)` pairs corresponding to original and transformed retention times::

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
TrafoXML parsers do not support indexing because there are no IDs for specific data points in this format.

Controlled Vocabularies
=======================

`Controlled Vocabularies <http://www.psidev.info/controlled-vocabularies>`_
are the universal annotation system used in the PSI formats, including
**mzML** and **mzIdentML**. :py:class:`pyteomics.mzml.MzML`, :py:class:`pyteomics.traml.TraML` and :py:class:`pyteomics.mzid.MzIdentML`
retain the annotation information. It can be accessed using the helper function, :py:func:`pyteomics.auxiliary.cvquery`::

    >>> from pyteomics import auxiliary as aux, mzid, mzml
    >>> f = mzid.MzIdentML('tests/test.mzid')
    >>> s = next(f)
    >>> s
    {'SpectrumIdentificationItem': [{'ProteinScape:SequestMetaScore': 7.59488518903425, 'calculatedMassToCharge': 1507.695, 'PeptideEvidenceRef': [{'peptideEvidence_ref': 'PE1_SEQ_spec1_pep1'}], 'chargeState': 1, 'passThreshold': True, 'peptide_ref': 'prot1_pep1', 'rank': 1, 'id': 'SEQ_spec1_pep1', 'ProteinScape:IntensityCoverage': 0.3919545603809718, 'experimentalMassToCharge': 1507.696}], 'spectrumID': 'databasekey=1', 'id': 'SEQ_spec1', 'spectraData_ref': 'LCMALDI_spectra'}
    >>> aux.cvquery(s)
    {'MS:1001506': 7.59488518903425, 'MS:1001505': 0.3919545603809718}
    >>> f.close()
