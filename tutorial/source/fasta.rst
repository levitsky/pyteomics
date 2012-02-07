FASTA
=====

To extract data from FASTA databases, use the :py:func:`read_fasta` function
from :py:mod:`fasta.py`. It’s as easy as it gets: you just need to specify
the file::

    from pyteomics.fasta import read_fasta
    proteins = read_fasta(‘/path/to/file/my.fasta’)

:py:func:`read_fasta` returns a *generator object* instead of a *list*
to prevent excessive memory use. To get a list if you actually need it,
type::

    protein_list = [p for p in read_fasta(‘/path/to/file/my.fasta’)]


You can also create a FASTA file using a list of (description, sequence) *tuples*.

::

    >>> from pyteomics.fasta import write_fasta
    >>> entries = [(“Protein 1”, ‘PEPTIDE’*1000), (“Protein 2”, ‘PEPTIDE’*2000)]
    >>> write_fasta(entries, ‘target-file.fasta’)

A common task is to generate a so-called *decoy database*. **Pyteomics** allows
that by means of the :py:func:`decoy_db` function.  ::

    >>> from pyteomics.fasta import decoy_db
    >>> decoy_db(‘mydb.fasta’, ‘mydb-with-decoy.fasta’)

The only required argument is the first one, indicating the source database. The
second argument is the target file and defaults to system standard output. 

.. seealso::

    See the documentation for more info about possible parameters.

If you need to modify a single sequence, use the :py:func:`decoy_sequence`
method. It currently supports to modes: *‘reverse’* and *‘random’*.

    >>> from pyteomics.fasta import decoy_sequence
    >>> decoy_sequence('PEPTIDE', 'reverse')
    'EDITPEP'
    >>> decoy_sequence('PEPTIDE', 'random')
    ‘TPPIDEE'
    >>> decoy_sequence('PEPTIDE', 'random')
    'PTIDEPE'

