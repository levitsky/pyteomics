Peptide sequence formats. *Parser* module
=========================================

Everywhere in **Pyteomics** we use a widely established way of peptide
sequence notation. We call it the **modX notation**, meaning that each amino
acid residue is represented with a capital letter representing a standard amino
acid, preceded by an arbitrary number of small letters to designate a
modification. Terminal groups can be specified using a hyphen (‘-’).
*“H-HoxMMdaN”* is an example of valid sequence notation. See 
`parser documentation page` for additional information. Also, *parser* module
contains some conversion tools for sequence information (see below).

Examples
--------

:py:func:`peptide_length` is the easy way to calculate the peptide length. It
can process sequences with and without explicit terminal groups correctly::

    >>> from pyteomics.parser import peptide_length
    >>> peptide_length(‘H-PEPTIDE’) == peptide_length(‘PEPTIDE’)
    True

A *modX* sequence can be translated to a list of amino acid residues with
:py:func:`parse_sequence` function::

    >>> from pyteomics.parser import parse_sequence
    >>> for label in parse_sequence('PEPTIDE'):
       print label
    P
    E
    P
    T
    I
    D
    E
    >>> from pyteomics.parser import std_labels # std_labels is a pre-defined list of standard amino acid labels
    >>> parse_sequence('aVRILLaVIGNE', labels=std_labels+['aV'])
    ['aV', 'R', 'I', 'L', 'L', 'aV', 'I', 'G', 'N', 'E']

See more examples and information in the documentation.

.. note::

    You don’t have to import each function separately. See, for example,
    http://effbot.org/zone/import-confusion.htm.

The :py:func:`amino_acid_composition` function accepts a sequence and returns
a *dictionary* with amino acid labels as *keys* and integer numbers as *values*,
corresponding to the number of times each residue occurs in the sequence::

    >>> from pyteomics.parser import amino_acid_composition
    >>> amino_acid_composition('PEPTIDE')
    {'I': 1.0, 'P': 2.0, 'E': 2.0, 'T': 1.0, 'D': 1.0}

Lastly, :py:func:`cleave` is a method to perform *in silico* cleavage. You need
to specify the sequence, the cleavage rule and the number of missed cleavages
allowed::

    >>> from pyteomics.parser import cleave, expasy_rules
    >>> cleave('AKAKBK', expasy_rules['trypsin'], 0)
    ['AK', 'BK']

:py:data:`expasy_rules` is a *dictionary* with the most common rules of cleavage.
It’s contained in the :py:mod:`pyteomics.parser` module, too. As you can see, the
returned *list* contains unique peptides.
The following example function finds the common fragments occurring in
proteolytic digests of all polypeptides in a given list::

    def digest_overlap(peptides, rule, missed_cleavages):
        sets = [set(cleave(peptide, rule, missed_cleavages)) for peptide in peptides]
        return set.intersection(*sets)

Then one can invoke it like this::

    >>> digest_overlap(['PEKPKTIKDKE', 'TIKDEKPEKP'], expasy_rules['trypsin'], 1)
    set(['TIK'])

This means, there is only one common peptide in the two digests in the specified conditions.
