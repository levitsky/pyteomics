Data access
===========

The following section is dedicated to data manipulation. **Pyteomics** aims to 
support the most common formats of (LC-)MS/MS data, peptide identification
files and protein databases 

mzML
----

**mzML** is an XML-based format for experimental data obtained on MS/MS or LC-MS
setups. **Pyteomics** offers you the functionality of :py:mod:`pyteomics.mzml` module
to gain access to the information contained in .mzML files from Python.

The main function in this module is
:py:func:`pyteomics.mzml.iter_spectrum`. It allows the user to iterate through MS/MS spectra
contained in an mzML file. Here is an example of its output:

.. code-block:: python

    >>> from pyteomics import mzml
    >>> reader = mzml.iter_spectrum('tests/test.mzML')
    >>> print reader.next() # Retrieve the first spectrum from the file and print it.
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
     'ms level': 1.0,
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
   
At the moment, the interface of :py:func:`pyteomics.mzml.iter_spectrum` is
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
:py:func:`iter_spectrum` function. It allows iterating through spectrum entries.
Spectra are represented as :py:class:`dictionaries`. MS/MS peak lists are stored
as :py:class:`numpy.array` objects :py:obj:`masses` and :py:obj:`intensities`.
Parameters are stored as a :py:class:`dict` under 'params' key.

Here is an example of use:

.. code-block:: python

    >>> from pyteomics import mgf
    >>> reader = mgf.iter_spectrum('tests/test.mgf')
    >>> print reader.next() # Retrieve the first spectrum from the file and print it.
    {'intensities': array([  73.,   44.,   67.,  291.,   54.,   49.]), 
    'masses': array([  846.6,   846.8,   847.6,  1640.1,  1640.6,  1895.5]), 
    'params': {'username': 'Lou Scene', 'useremail': 'leu@altered-state.edu',
    'mods': 'Carbamidomethyl (C)', 'itolu': 'Da', 'title': 'Spectrum 1',
    'itol': '1', 'charge': '2+ and 3+', 'mass': 'Monoisotopic',
    'it_mods': 'Oxidation (M)', 'pepmass': '983.6',
    'com': 'Taken from http://www.matrixscience.com/help/data_file_help.html'}}

Also, :py:mod:`pyteomics.mgf` allows to extract headers with general search 
parameters from MGF files with :py:func:`read_header` function. It also returns
a :py:class:`dict`.

.. code-block:: python

    >>> header = mgf.read_header('tests/test.mgf')
    >>> print header
    {'username': 'Lou Scene', 'itol': '1', 'useremail': 'leu@altered-state.edu',
    'mods': 'Carbamidomethyl (C)', 'it_mods': 'Oxidation (M)',
    'charge': '2+ and 3+', 'mass': 'Monoisotopic', 'itolu': 'Da',
    'com': 'Taken from http://www.matrixscience.com/help/data_file_help.html'}


Creation of MGF files is implemented in :py:func:`write_mgf` function. The user
can specify the header, list of spectra in the same format as returned by
:py:func:`iter_spectrum` and the output path.

.. code-block:: python

    >>> spectra = [s for s in mgf.iter_spectrum('tests/test.mgf')]
    >>> mgf.write_mgf(spectra=spectra, header=header)
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
used and the assigned sequences. To access these data, use :py:mod:`pyteomics.pepxml`
module.

:py:mod:`pyteomics.pepxml` has the same structure as :py:mod:`pyteomics.mzml`. The function
:py:func:`pyteomics.pepxml.iter_psm` iterates through Peptide-Spectrum matches in a .pepXML file 
and returns them as a custom dict.

.. code-block:: python

    >>> from pyteomics import pepxml
    >>> reader = pepxml.iter_psm('tests/test.pep.xml')
    >>> print reader.next()
    {'end_scan': 100.0,
    'search_hits': [{'hit_rank': 1.0,
    'calc_neutral_pep_mass': 860.892,
    'deltacn': 0.081,
    'modified_peptide': 'SLNGEWR',
    'peptide': 'SLNGEWR',
    'num_matched_ions': 11.0,
    'modifications': [],
    'peptideprophet': 0.96,
    'proteins': [{'num_tol_term': 2.0,
    'protein': 'sp|P00722|BGAL_ECOLI',
    'peptide_prev_aa': 'R',
    'protein_descr': 'BETA-GALACTOSIDASE (EC 3.2.1.23) (LACTASE) - Escherichia coli.',
    'peptide_next_aa': 'F'}],
    'sprank': 1.0,
    'num_missed_cleavages': 0.0,
    'tot_num_ions': 12.0,
    'num_tot_proteins': 1.0,
    'deltacnstar': 0.0,
    'is_rejected': '0',
    'spscore': 894.0,
    'xcorr': 1.553,
    'massdiff': -0.5}],
    'index': 1.0,
    'assumed_charge': 1.0,
    'spectrum': 'pps_sl20060731_18mix_25ul_r1_1154456409.0100.0100.1',
    'precursor_neutral_mass': 860.392,
    'start_scan': 100.0}
                                                                                       
FASTA
-----

To extract data from FASTA databases, use the :py:func:`pyteomics.fasta.read_fasta` 
function.

.. code-block:: python

    >>> from pyteomics import fasta
    >>> proteins = list(fasta.read_fasta('/path/to/file/my.fasta'))

:py:func:`pyteomics.fasta.read_fasta` returns a *generator object* instead of a
:py:class:`list` to prevent excessive memory use. 

You can also create a FASTA file using a list of (description, sequence)
:py:class:`tuples`.

.. code-block:: python

    >>> from pyteomics import fasta
    >>> entries = [('Protein 1', 'PEPTIDE'*1000), ('Protein 2', 'PEPTIDE'*2000)]
    >>> fasta.write_fasta(entries, 'target-file.fasta')

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

