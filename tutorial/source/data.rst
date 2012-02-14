Data access
===========

The following section is dedicated to data manipulation. **Pyteomics** aims to 
support the most common formats of (LC-)MS/MS data, peptide identification
files and protein databases 

mzML
----

**mzML** is an XML-based format for experimental data obtained on MS/MS or LC-MS
setups. **Pyteomics** offers you the functionality of :py:mod:`mzml` module
to gain access to the information contained in .mzML files from Python.

Examples
........

The main function in this module is
:py:func:`iter_spectrum`. It allows the user to iterate through MS/MS spectra
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
   
At the moment, the interface of :py:func:`iter_spectrum` is relatively 
low-level. It iterates through the spectra in the file and returns each one as 
a dict with selected fields stored. The interface is relatively raw and can be
modified in the subsequent releases of Pyteomics.

pepXML
------

**.pepXML** is a widely used XML-based format for peptide identifications.
It contains information about the MS data, the parameters of the search engine 
used and the assigned sequences. To access these data, use :py:mod:`pepxml`
module.

:py:mod:`pepxml` has the same structure as :py:mod:`mzml`. The function
:py:func:`iter_psm` iterates through Peptide-Spectrum matches in a .pepXML file 
and returns them as a custom dict.

.. code-block:: python

    >>> from pyteomics import pepxml
    >>> reader = pepxml.iter_psm('tests/test.pep.xml')
    >>> print reader.next()
    {'assumed_charge': 1.0,
     'calc_neutral_pep_mass': 860.892,
     'deltacn': 0.081,
     'deltacnstar': 0.0,
     'end_scan': 100.0,
     'hit_rank': 1.0,
     'index': 1.0,
     'is_rejected': '0',
     'massdiff': -0.5,
     'modifications': [],
     'modified_peptide': 'SLNGEWR',
     'num_matched_ions': 11.0,
     'num_missed_cleavages': 0.0,
     'num_tot_proteins': 1.0,
     'peptide': 'SLNGEWR',
     'peptideprophet': 0.96,
     'precursor_neutral_mass': 860.392,
     'proteins': [{'num_tol_term': 2.0,
                   'peptide_next_aa': 'F',
                   'peptide_prev_aa': 'R',
                   'protein': 'sp|P00722|BGAL_ECOLI',
                   'protein_descr': 'BETA-GALACTOSIDASE (EC 3.2.1.23) (LACTASE) - Escherichia coli.'}],
     'spectrum': 'pps_sl20060731_18mix_25ul_r1_1154456409.0100.0100.1',
     'sprank': 1.0,
     'spscore': 894.0,
     'start_scan': 100.0,
     'tot_num_ions': 12.0,
     'xcorr': 1.553}
                                                                                       
FASTA
-----

To extract data from FASTA databases, use the :py:func:`read_fasta` function
from :py:mod:`fasta.py`.

.. code-block:: python

    >>> from pyteomics import fasta
    >>> proteins = list(fasta.read_fasta('/path/to/file/my.fasta'))

:py:func:`read_fasta` returns a *generator object* instead of a *list*
to prevent excessive memory use. 

You can also create a FASTA file using a list of (description, sequence) *tuples*.

.. code-block:: python

    >>> from pyteomics import fasta
    >>> entries = [('Protein 1', 'PEPTIDE'*1000), ('Protein 2', 'PEPTIDE'*2000)]
    >>> fasta.write_fasta(entries, 'target-file.fasta')

Another common task is to generate a *decoy database*. **Pyteomics** allows
that by means of the :py:func:`decoy_db` function. 

.. code-block:: python

    >>> from pyteomics import fasta
    >>> fasta.decoy_db('mydb.fasta', 'mydb-with-decoy.fasta')

The only required argument is the first one, indicating the source database. The
second argument is the target file and defaults to system standard output. 

If you need to modify a single sequence, use the :py:func:`decoy_sequence`
method. It currently supports two modes: *‘reverse’* and *‘random’*.

.. code-block:: python

    >>> from pyteomics import fasta
    >>> fasta.decoy_sequence('PEPTIDE', 'reverse')
    'EDITPEP'
    >>> fasta.decoy_sequence('PEPTIDE', 'random')
    ‘TPPIDEE'
    >>> fasta.decoy_sequence('PEPTIDE', 'random')
    'PTIDEPE'

