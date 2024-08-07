Peptide sequence formats. *Parser* module
=========================================

modX
----

**Pyteomics** uses a custom IUPAC-derived peptide sequence notation named **modX**.
As in the IUPAC notation, each amino acid residue is represented by a capital
letter, but it may preceded by an arbitrary number of small letters to show
modification. Terminal groups are separated from the backbone sequence by
a hyphen (‘-’). By default, both termini are assumed to be unmodified, which can be
shown explicitly by 'H-' for N-terminal hydrogen and '-OH' for C-terminal hydroxyl.

``“H-HoxMMdaN-OH”`` is an example of a valid sequence in *modX*. See
:doc:`api/parser` for additional information. Note that it is recommended to include
either 0 or 2 terminal groups in a *modX* sequence.

Sequence operations
-------------------

Parsing
.......

There are two helper functions to check if a label is in modX format or represents
a terminal modification: :py:func:`pyteomics.parser.is_modX` and
:py:func:`pyteomics.parser.is_term_group`:

.. code-block:: python

    >>> parser.is_modX('A')
    True
    >>> parser.is_modX('pT')
    True
    >>> parser.is_modX('pTx')
    False
    >>> parser.is_term_group('pT')
    False
    >>> parser.is_term_group('Ac-')
    True


A *modX* sequence can be translated to a list of amino acid residues with
:py:func:`pyteomics.parser.parse` function:

.. code-block:: python

    >>> from pyteomics import parser
    >>> parser.parse('PEPTIDE')
    ['P', 'E', 'P', 'T', 'I', 'D', 'E']
    >>> parser.parse('PEPTIDE', show_unmodified_termini=True)
    ['H-', 'P', 'E', 'P', 'T', 'I', 'D', 'E', '-OH']
    >>> parser.parse('Ac-PEpTIDE', labels=parser.std_labels+['Ac-', 'pT'])
    ['Ac-', 'P', 'E', 'pT', 'I', 'D', 'E']

In the last example we supplied two arguments, the sequence itself
and 'labels'. The latter is used to specify what labels are allowed for amino
acid residues and terminal modifications. :py:data:`std_labels` is a predefined
set of labels for the twenty standard amino acids, 'H-' for N-terminal hydrogen
and '-OH' for C-terminal hydroxyl. In this example we specified the codes for
phosphorylated threonine and N-terminal acetylation.

Since version 2.5, specifying ``labels`` is *never mandatory*. If this argument
is not supplied, no checks will be made. However, the last example won't work
without ``labels``, because it has **only one terminal group shown**, which is
**discouraged**.

:py:func:`parse` has another mode, in which it returns tuples:

.. code-block:: python

    >>> parser.parse('Ac-PEpTIDE-OH', split=True)
    [('Ac-', 'P'), ('E',), ('p', 'T'), ('I',), ('D',), ('E',)]

or:

.. code-block:: python

    >>> parser.parse('Ac-PEpTIDE-OH', split=True, labels=parser.std_labels+['Ac-', 'p'])
    [('Ac-', 'P'), ('E',), ('p', 'T'), ('I',), ('D',), ('E',)]

Also, note what we supply as `labels` here: 'p' instead of 'pT'. That means that
'p' is a modification applicable to any residue.

In modX, standard :py:func:`len` function cannot be used to determine the length
of a peptide because of the modifications.
Use :py:func:`pyteomics.parser.length` instead:

.. code-block:: python

    >>> from pyteomics import parser
    >>> parser.length('aVRILLaVIGNE')
    10

The :py:func:`pyteomics.parser.amino_acid_composition` function accepts a sequence
and returns a *dictionary* with amino acid labels as *keys* and integer numbers as
*values*, corresponding to the number of times each residue occurs in the sequence:

.. code-block:: python

    >>> from pyteomics import parser
    >>> parser.amino_acid_composition('PEPTIDE')
    {'I': 1.0, 'P': 2.0, 'E': 2.0, 'T': 1.0, 'D': 1.0}

*In silico* digestion
.....................

:py:func:`pyteomics.parser.cleave` performs *in silico* cleavage.
The required arguments are the sequence, the rule for enzyme specificity and the
number of missed cleavages allowed (optional). :py:func:`cleave` returns a
:py:class:`set` of product peptides; you can get original indices of peptides with :py:func:`xcleave`.

.. code-block:: python

    >>> from pyteomics import parser
    >>> parser.cleave('AKAKBK', parser.expasy_rules['trypsin'], 0)
    {'AK', 'BK'}
    >>> parser.xcleave('AKAKBK', 'trypsin', 0)
    [(0, 'AK'), (2, 'AK'), (4, 'BK')]

:py:data:`pyteomics.parser.expasy_rules` and :py:data:`pyteomics.parser.psims_rules` are predefined :py:class:`dicts`
with the clevage rules for the most common proteases. Their keys are recognized by :py:func:`cleave`.

Variable modifications
......................

All possible modified sequences of a peptide can be obtained with
:py:func:`pyteomics.parser.isoforms`:

.. code-block:: python

    >>> from pyteomics import parser
    >>> forms = parser.isoforms('PEPTIDE', variable_mods={'p': ['T'], 'ox': ['P']})
    >>> for seq in forms: print seq
    ...
    oxPEPpTIDE
    oxPEPTIDE
    oxPEoxPpTIDE
    oxPEoxPTIDE
    PEPpTIDE
    PEPTIDE
    PEoxPpTIDE
    PEoxPTIDE
