Data access
===========

The following section is dedicated to data manipulation. **Pyteomics** aims to
support the most common formats of (LC-)MS/MS data, peptide identification
files and protein databases.

mzML
----

**mzML** is an XML-based format for experimental data obtained on MS/MS or LC-MS
setups. **Pyteomics** offers you the functionality of :py:mod:`pyteomics.mzml`
module to gain access to the information contained in mzML files from Python.

The main function in this module is :py:func:`pyteomics.mzml.read`. It allows
the user to iterate through MS/MS spectra
contained in an mzML file. Here is an example of its output:

.. code-block:: python

    >>> from pyteomics import mzml
    >>> with mzml.read('tests/test.mzML') as reader:
    >>>     print next(reader) # Retrieve the first spectrum from the file and print it.
    {'MSn spectrum': '',
     'base peak intensity': 1471973.875,
     'base peak m/z': 810.415283203125,
     'defaultArrayLength': '19914',
     'highest observed m/z': 2000.0099466203771,
     'id': 'controllerType=0 controllerNumber=1 scan=1',
     'index': '0',
     'intensity array': array([ 0.,  0.,  0., ...,  0.,  0.,  0.], dtype=float32),
     'lowest observed m/z': 200.00018816645022,
     'm/z array': array([  200.00018817,   200.00043034,   200.00067252, ...,  1999.96151259,
                          1999.98572931,  2000.00994662]),
     'ms level': 1,
     'no combination': '',
     'positive scan': '',
     'precursorList': [],
     'profile spectrum': '',
     'scanList': [{'[Thermo Trailer Extra]Monoisotopic M/Z:': '810.41522216796875',
                   'filter string': 'FTMS + p ESI Full ms [200.00-2000.00]',
                   'preset scan configuration': '1',
                   'scan start time': '0.0049350000000000002',
                   'scanWindowList': [{'scan window lower limit': '200',
                                       'scan window upper limit': '2000'}]}],
     'total ion current': 15245068.0}
   
At the moment, the interface of :py:func:`pyteomics.mzml.read` is
relatively low-level. It iterates through the spectra in the file and returns
each one as a dict with selected fields stored. The interface is relatively
raw and can be modified in the subsequent releases of Pyteomics.

MGF
---

Mascot Generic Format
(`MGF <http://www.matrixscience.com/help/data_file_help.html>`_) is a simple
human-readable format for MS/MS data. It allows storing MS/MS peak lists and
exprimental parameters. :py:mod:`pyteomics.mgf` is a module that implements
reading and writing MGF files.

Similar to :py:mod:`pyteomics.mzml`, :py:mod:`pyteomics.mgf` has a
:py:func:`read` function. It allows iterating through spectrum entries.
Spectra are represented as :py:class:`dictionaries`. MS/MS peak lists are stored
as :py:class:`numpy.array` objects :py:obj:`masses` and :py:obj:`intensities`.
Parameters are stored as a :py:class:`dict` under 'params' key.

Here is an example of use:

.. code-block:: python

    >>> from pyteomics import mgf
    >>> with mgf.read('tests/test.mgf') as reader:
    >>>     print next(reader) # Retrieve the first spectrum from the file and print it.
    {'m/z array': array([  345.1,   370.2,   460.2,  1673.3,  1674. ,  1675.3]),
    'charge array': array([ 3,  2,  1,  1,  1,  1]),
    'params': {'username': 'Lou Scene', 'useremail': 'leu@altered-state.edu',
    'mods': 'Carbamidomethyl (C)', 'itolu': 'Da', 'title': 'Spectrum 2',
    'rtinseconds': '25', 'itol': '1', 'charge': '2+ and 3+',
    'mass': 'Monoisotopic', 'it_mods': 'Oxidation (M)',
    'pepmass': (1084.9, 1234.0),
    'com': 'Based on http://www.matrixscience.com/help/data_file_help.html',
    'scans': '3'},
    'intensity array': array([  237.,   128.,   108.,  1007.,   974.,    79.])}
    
Also, :py:mod:`pyteomics.mgf` allows to extract headers with general
parameters from MGF files with :py:func:`read_header` function. It also returns
a :py:class:`dict`.

.. code-block:: python

    >>> header = mgf.read_header('tests/test.mgf')
    >>> print header
    {'username': 'Lou Scene', 'itol': '1', 'useremail': 'leu@altered-state.edu',
    'mods': 'Carbamidomethyl (C)', 'it_mods': 'Oxidation (M)',
    'charge': '2+ and 3+', 'mass': 'Monoisotopic', 'itolu': 'Da',
    'com': 'Taken from http://www.matrixscience.com/help/data_file_help.html'}


Creation of MGF files is implemented in :py:func:`write` function. The user
can specify the header, list of spectra in the same format as returned by
:py:func:`read` and the output path.

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


pepXML
------

`pepXML <http://tools.proteomecenter.org/wiki/index.php?title=Formats:pepXML>`_
is a widely used XML-based format for peptide identifications.
It contains information about the MS data, the parameters of the search engine
used and the assigned sequences. To access these data, use
:py:mod:`pyteomics.pepxml` module.

:py:mod:`pyteomics.pepxml` has the same structure as :py:mod:`pyteomics.mzml`.
The function :py:func:`pyteomics.pepxml.read` iterates through Peptide-Spectrum
matches in a pepXML file and returns them as a custom dict.

.. code-block:: python

    >>> from pyteomics import pepxml
    >>> reader = pepxml.read('tests/test.pep.xml')
    >>> print next(reader)
    {'end_scan': 100,
    'index': 1,
    'assumed_charge': 1,
    'spectrum': 'pps_sl20060731_18mix_25ul_r1_1154456409.0100.0100.1',
    'search_hit': [
        {'hit_rank': 1,
        'calc_neutral_pep_mass': 860.892,
        'modifications': [],
        'modified_peptide': 'SLNGEWR',
        'peptide': 'SLNGEWR',
        'num_matched_ions': 11,
        'search_score': {
            'spscore': 894.0,
            'sprank': 1.0,
            'deltacnstar': 0.0,
            'deltacn': 0.081,
            'xcorr': 1.553},
        'proteins': [
            {'num_tol_term': 2,
            'protein': 'sp|P00722|BGAL_ECOLI',
            'peptide_next_aa': 'F',
            'protein_descr': 'BETA-GALACTOSIDASE (EC 3.2.1.23) (LACTASE) - Escherichia coli.',
            'peptide_prev_aa': 'R'}],
        'num_missed_cleavages': 0,
        'analysis_result': [
            {'peptideprophet_result':
                {'parameter': {'massd': -0.5, 'nmc': 0.0, 'ntt': 2.0, 'fval': 1.4723},
                'all_ntt_prob': [0.0422, 0.509, 0.96],
                'probability': 0.96},
                'analysis': 'peptideprophet'}],
        'tot_num_ions': 12,
        'num_tot_proteins': 1,
        'is_rejected': False,
        'massdiff': -0.5}],
    'precursor_neutral_mass': 860.392,
    'start_scan': 100}

X!Tandem
--------

`X!Tandem search engine <http://www.thegpm.org/tandem/>`_ has its own output
format that contains more info than pepXML. **Pyteomics** has a reader for it
in the :py:mod:`pyteomics.tandem` module.

.. code-block:: python

    >>> from pyteomics import tandem
    >>> with tandem.read('tests/test.t.xml') as reader:
    ...     print next(reader)
    ...
    {'support':
        {'fragment ion mass spectrum':
            {'note': 'scan=9161 cs=2', 'charge': 2, 'Ydata':
                {'units': 'UNKNOWN',
                'values':
                array([   9.,   23.,   13.,   12.,   11.,   24.,    9.,  100.,    7.,
                 18.,   18.,    9.,   12.,    8.,    7.,   10.,   21.,   15.,
                 24.,   14.,   20.,   10.,   14.,    7.,   13.,    9.,    7.,
                 13.,    7.,    8.,    7.,   21.,   13.,   10.,   12.,   84.,
                 21.,   34.,   16.,   11.,   13.,   43.,   12.,    7.,   21.,
                 11.,   12.,   28.,    9.,    8.])}, 'Xdata':
                {'units': 'MASSTOCHARGERATIO',
                'values':
                array([  424.241,   460.272,   553.307,   575.344,   597.315,   665.184,
                 725.412,   743.416,   745.566,   794.265,   812.213,   813.284,
                 872.11 ,   873.41 ,   883.199,   907.471,   911.529,   925.485,
                 928.72 ,   937.587,   994.221,  1012.29 ,  1029.29 ,  1043.54 ,
                1057.58 ,  1084.37 ,  1093.22 ,  1116.39 ,  1243.35 ,  1261.63 ,
                1285.44 ,  1302.6  ,  1369.64 ,  1373.52 ,  1386.46 ,  1389.52 ,
                1391.66 ,  1503.76 ,  1504.75 ,  1572.66 ,  1613.72 ,  1631.67 ,
                1633.76 ,  1728.85 ,  1745.87 ,  1746.87 ,  1837.87 ,  1855.79 ,
                1857.   ,  1969.52 ])},
             'M+H': 2314.84, 'id': '10745', 'label': '10745.spectrum'},
         'supporting data':
            {'b ion histogram':
               {'Ydata':
                {'units': 'counts',
                 'values': array([2736, 3890, 1074,  201,   25,    0,    1,    0])},
                'Xdata':
                {'units': 'number of ions',
                'values': array([0, 1, 2, 3, 4, 5, 6, 7])},
                'label': '10745.b'},
             'y ion histogram':
                {'Ydata':
                 {'units': 'counts',
                  'values': array([2890, 3802, 1013,  197,   16,    5,    1,    0,    3,    0])},
                 'Xdata':
                 {'units': 'number of ions',
                 'values': array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])},
                 'label': '10745.y'},
             'convolution survival function':
                {'Ydata':
                 {'units': 'counts',
                  'values': array([7927, 7927, 7927, 7927, 7562, 6422, 5231, 3521, 2005,  833,  336, 72,   11,    2,    0])},
                 'Xdata':
                 {'units': 'score',
                 'values': array([ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14])},
                 'label': '10745.convolute'},
             'hyperscore expectation function':
                {'label': '10745.hyper', 'a1': -0.266291, 'a0': 5.60936,
                'Ydata':
                 {'units': 'counts',
                  'values': array([7922, 7922, 7922, 7922, 7557, 6417, 5382, 3923, 2715, 1693, 1057,
                                    572,  322,  180,   89,   35,   15,    7,    0,    0,    0,    0,
                                      5,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
                                      0,    0,    0,    3,    0])},
                'Xdata':
                {'units': 'score',
                 'values': array([ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16,
                                  17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33,
                                  34, 35, 36, 37])}}}},
     'mh': 2314.835, 'maxI': 22514.9, 'expect': 8.4e-05, 'sumI': 5.29,
     'fI': 225.149,
     'protein':
        {'file':
            {'URL': 'C:\\DMS_Temp_Org\\ID_001140_4BD5AF39.fasta',
             'type': 'peptide'},
         'peptide':
            {'pre': 'VVAR', 'seq': 'EQALQIEISMNEGKPADIAEK',
            'b_ions': 4, 'y_ions': 8, 'hyperscore': 36.4, 'b_score': 8.7,
            'expect': 8.4e-05, 'delta': 0.676, 'post': 'MVVG', 'id': '10745.1.1',
            'end': 211, 'mh': 2314.159, 'missed_cleavages': 0, 'start': 191,
            'y_score': 12.2, 'nextscore': 21.8}, 'uid': '1370',
            'label': 'SO_1630 translation elongation factor Ts (Tsf)',
            'note': 'SO_1630 translation elongation factor Ts (Tsf)',
            'expect': -137.0, 'sumI': 6.79, 'id': '10745.1'},
            'z': 2, 'id': '10745'}

mzIdentML
---------

`mzIdentML <http://www.psidev.info/mzidentml>`_  is one of the standards
developed by the Proteomics Informatics working group of the HUPO Proteomics
Standard Initiative.

The module interface is similar to that of the other reader modules.

.. code-block:: python

    >>> from pyteomics import mzid
    >>> with mzid.read('tests/test.mzid') as reader:
    >>>     print next(reader)
    {'SpectrumIdentificationItem': [
        {'ProteinScape:IntensityCoverage': 0.3919545603809718,
        'PeptideEvidenceRef': [
            {'peptideEvidence_ref': 'PE1_SEQ_spec1_pep1'}],
        'passThreshold': True,
        'rank': 1,
        'chargeState': 1,
        'calculatedMassToCharge': 1507.695,
        'peptide_ref': 'prot1_pep1',
        'experimentalMassToCharge': 1507.696,
        'id': 'SEQ_spec1_pep1',
        'ProteinScape:SequestMetaScore': 7.59488518903425}],
    'spectrumID': 'databasekey=1',
    'id': 'SEQ_spec1',
    'spectraData_ref': 'LCMALDI_spectra'}

You can tune the amount of information you get from the file. The available
options to the :py:func:`read` function are `recursive` (:py:const:`True` by
default) and `retrieve_refs` (:py:const:`False` by default). The latter pulls
additional info from the file that is present only as references in the example
above.

Additional function :py:func:`get_by_id` allows to extract info from any element
using its unique ID.

FASTA
-----

To extract data from FASTA databases, use the :py:func:`pyteomics.fasta.read`
function.

.. code-block:: python

    >>> from pyteomics import fasta
    >>> proteins = list(fasta.read('/path/to/file/my.fasta'))

:py:func:`pyteomics.fasta.read` returns a *generator object* instead of a
:py:class:`list` to prevent excessive memory use. The generator yields
(description, sequence) tuples, so it's natural to use it as follows:

.. code-block:: python

    >>> from pyteomics import fasta
    >>> for descr, seq in fasta.read('my.fasta'):
    >>>    ...

You can also use attributes to access description and sequence:

.. code-block:: python

    >>> from pyteomics import fasta
    >>> for protein in fasta.read('my.fasta'):
    >>>    print protein.description
    >>>    print protein.sequence
    
Note the new recommended `with` syntax:

.. code-block:: python

    >>> from pyteomics import fasta
    >>> with fasta.read('my.fasta') as reader:
    >>>    for descr, seq in reader:
    >>>       ...

You can specify a function that will be applied to the FASTA headers for
your convenience. :py:data:`pyteomics.fasta.std_parsers` has some pre-defined
parsers that can be used for this purpose.

You can also create a FASTA file using a sequence of (description, sequence)
:py:class:`tuples`.

.. code-block:: python

    >>> from pyteomics import fasta
    >>> entries = [('Protein 1', 'PEPTIDE'*1000), ('Protein 2', 'PEPTIDE'*2000)]
    >>> fasta.write(entries, 'target-file.fasta')

Another common task is to generate a *decoy database*. **Pyteomics** allows
that by means of the :py:func:`pyteomics.fasta.decoy_db` function.

.. code-block:: python

    >>> from pyteomics import fasta
    >>> fasta.decoy_db('mydb.fasta', 'mydb-with-decoy.fasta')

The only required argument is the first one, indicating the source database. The
second argument is the target file and defaults to system standard output.

If you need to modify a single sequence, use the :py:func:`pyteomics.fasta.decoy_sequence`
method. It currently supports two modes: *‘reverse’* and *‘random’*.

.. code-block:: python

    >>> from pyteomics import fasta
    >>> fasta.decoy_sequence('PEPTIDE', 'reverse')
    'EDITPEP'
    >>> fasta.decoy_sequence('PEPTIDE', 'random')
    ‘TPPIDEE'
    >>> fasta.decoy_sequence('PEPTIDE', 'random')
    'PTIDEPE'

