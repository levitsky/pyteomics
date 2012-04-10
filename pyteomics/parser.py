"""
parser - operations on modX peptide sequences
=============================================

modX is a simple extension of the `IUPAC one-letter peptide sequence
representation <http://www.chem.qmul.ac.uk/iupac/AminoAcid/A2021.html>`_.

The labels (or codes) for the 20 standard amino acids in modX are the same as
in IUPAC nomeclature. A label for a modified amino acid has a general
form of 'modX', i.e.:

- it starts with an arbitrary number of lower-case symbols or numbers
  (a modification);

- it ends with a single upper-case symbol (an amino acid residue).

The valid examples of modX amino acid labels are: 'G', 'pS', 'oxM'. This rule
allows to combine read- and parseability.

Besides the sequence of amino acid residues, modX has a rule to specif
terminal modifications of a polypeptide. Such a label should start or
end with a hyphen. The default N-terminal amine group and C-terminal
carboxyl group may not be shown explicitly.

Therefore, the valid examples of peptide sequences in modX are: "GAGA",
"H-PEPTIDE", "TEST-NH2".

Operations on polypeptide sequences
-----------------------------------

  :py:func:`parse_sequence` - convert a sequence string into a list of
  amino acid residues.
  
  :py:func:`amino_acid_composition` - get numbers of each amino acid
  residue in a peptide.
  
  :py:func:`cleave` - cleave a polypeptide using a given rule of
  enzymatic digestion.

  :py:func:`isoforms` - find the set of unique modified peptides
  given the initial sequence and modifications.

Auxiliary commands
------------------

  :py:func:`peptide_length` - calculate the number of amino acid
  residues in a polypeptide.

  :py:func:`is_term_mod` - check if supplied code corresponds to a
  terminal modification.

Data
----

  :py:data:`std_amino_acids` - a list of the 20 standard amino acid IUPAC codes.
  
  :py:data:`std_nterm` - the standard N-terminal modification (the
  unmodified group is a single atom of hydrogen).
  
  :py:data:`std_cterm` - the standard C-terminal modification (the
  unmodified group is hydroxyl).

  :py:data:`std_labels` - a list of all standard sequence
  elements, amino acid residues and terminal modifications.

  :py:const:`expasy_rules` - a dict with the regular expressions of
  cleavage rules for the most popular proteolytic enzymes.

.. ipython::
   :suppress:

   In [1]: import pyteomics.parser; from pprint import pprint

-------------------------------------------------------------------------------

"""

# Licensed under the MIT license:
# http://www.opensource.org/licenses/mit-license.php

import re
from collections import deque
from itertools import chain
from .auxiliary import PyteomicsError

std_amino_acids = ['Q','W','E','R','T','Y','I','P','A','S',
                   'D','F','G','H','K','L','C','V','N','M']
"""modX labels for the 20 standard amino acids."""

std_nterm = 'H-'
"""modX label for the unmodified N-terminus."""

std_cterm = '-OH'
"""modX label for the unmodified C-terminus."""

std_labels = std_amino_acids + [std_nterm, std_cterm]
"""modX labels for the standard amino acids and unmodified termini."""

def is_term_mod(label):
    """Check if `label` corresponds to a terminal modification.

    Parameters
    ----------
    label : str

    Returns
    -------
    out : bool
    """
    return label.startswith('-') or label.endswith('-')

def is_modX(label):
    """Check if `label` is a valid 'modX' label.

    Parameters
    ----------
    label : str

    Returns
    -------
    out : bool
    """
    return label and ((len(label) == 1 and label[0].isupper()) or
            (label[:-1].islower() and label[-1].isupper()))

def peptide_length(sequence, **kwargs):
    """Calculate the number of amino acid residues in a polypeptide
    written in modX notation.

    Parameters
    ----------
    sequence : str or list or dict
        A string with a polypeptide sequence, a list with a parsed sequence or
        a dict of amino acid composition.        
    labels : list, optional    
        A list of allowed labels for amino acids and terminal modifications
        (default is `std_labels`, the 20 standard amino acids, N-terminal H-
        and C-terminal -OH).

    Examples
    --------
    >>> peptide_length('PEPTIDE')
    7
    >>> peptide_length('H-PEPTIDE-OH')
    7
    """    
    labels = kwargs.get('labels', std_labels)

    if isinstance(sequence, str) or isinstance(sequence, list):
        if isinstance(sequence, str):
            parsed_sequence = parse_sequence(sequence, labels=labels)
        else:
            parsed_sequence = sequence            
        num_term_groups = 0
        if is_term_mod(parsed_sequence[0]):
            num_term_groups += 1
        if is_term_mod(parsed_sequence[-1]):
            num_term_groups += 1
        return len(parsed_sequence) - num_term_groups
    elif isinstance(sequence, dict):
        return sum([amount for aa, amount in list(sequence.items()) 
                    if not is_term_mod(aa)])

    raise PyteomicsError('Unsupported type of a sequence.')
    return None

def parse_sequence(sequence,               
                   show_unmodified_termini=False, split=False,
                   **kwargs):
    """Parse a sequence string written in modX notation into a list of
    tuples representing amino acid residues and their modifications. 

    Parameters
    ----------
    sequence : str
        The sequence of a polypeptide.
    show_unmodified_termini : bool, optional   
        If :py:const:`True` then the unmodified N- and C-termini are explicitly
        shown in the returned list. Default value is :py:const:`False`.
    split : bool, optional
        If :py:const:`True` then the result will be a list of tuples with 1 to 3
        elements: terminal modification, modification, residue. Default value is
        :py:const:`False`.
    labels : list, optional
        A list of allowed labels for amino acids and terminal modifications
        (default is the 20 standard amino acids, N-terminal H- and C-terminal
        -OH).

    Returns
    -------
    out : list
        List of tuples with labels of modifications and amino acid residues.

    Examples
    --------
    >>> parse_sequence('PEPTIDE', split=True)
    [('P',), ('E',), ('P',), ('T',), ('I',), ('D',), ('E',)]
    >>> parse_sequence('H-PEPTIDE')
    ['P', 'E', 'P', 'T', 'I', 'D', 'E']
    >>> parse_sequence('PEPTIDE', show_unmodified_termini=True)
    ['H-', 'P', 'E', 'P', 'T', 'I', 'D', 'E', '-OH']
    >>> parse_sequence('TEpSToxM', labels=std_labels + ['pS', 'oxM'])
    ['T', 'E', 'pS', 'T', 'oxM']
    """
    def split_label(label):
        if not is_modX(label):
            raise PyteomicsError('Cannot split a non-modX label: %s' % label)
        if len(label) == 1:
            return (label, )
        else:
            return (label[:-1], label[-1])

    labels = kwargs.get('labels', std_labels)
    backbone_sequence = str(sequence)

    # Determine the N-terminal modification.
    nterm = std_nterm
    for aa in labels:
        if is_term_mod(aa) and backbone_sequence.startswith(aa):
            nterm = aa
            backbone_sequence = backbone_sequence[len(nterm):]
            break

    # Determine the C-terminal modification.
    cterm = std_cterm
    j = backbone_sequence.find('-')
    if j != -1:
        for aa in labels:
            if is_term_mod(aa) and backbone_sequence[j:] == aa:
                cterm = aa
                backbone_sequence = backbone_sequence[:j]
                break

    # Parse the polypeptide backbone.
    i = 0
    parsed_sequence = []
    while i < len(backbone_sequence):
        amino_acid_found = False
        for aa in labels:
            if backbone_sequence.startswith(aa, i) and is_modX(aa):
                parsed_sequence.append(aa)
                amino_acid_found = True
                break

        j = i+2
        while j <= len(backbone_sequence):
            try:
                mod, res = split_label(backbone_sequence[i:j])
            except PyteomicsError:
                pass
            else:
                if mod in labels and res in labels:
                    parsed_sequence.append(backbone_sequence[i:j])
                    amino_acid_found = True
                    break
            finally:
                j += 1

        if not amino_acid_found:
            raise PyteomicsError(
                'Unknown amino acid in sequence %s at position %d: %s' % (
                    backbone_sequence, i+1, backbone_sequence[i:]))
        i += len(parsed_sequence[-1])

    # Append labels of unmodified termini.
    if show_unmodified_termini:
        parsed_sequence.insert(0, nterm)
        parsed_sequence.append(cterm)

    # Make a list of tuples instead of list of labels
    if split:
        if not parsed_sequence or all(map(is_term_mod, parsed_sequence)):
            return map(tuple, parsed_sequence)

        
        tuples = []
        start = 0
        if is_term_mod(parsed_sequence[0]):
            tuples.append((parsed_sequence[0],) + split_label(parsed_sequence[1]))
            start = 2
        tuples.extend(split_label(x) for x in parsed_sequence[start:-1])
        if is_term_mod(parsed_sequence[-1]):
            tuples.append(split_label(parsed_sequence[-2]) + (parsed_sequence[-1],))
        else:
            tuples.append(split_label(parsed_sequence[-1]))
        
        return tuples

    return parsed_sequence

def amino_acid_composition(sequence,
                           show_unmodified_termini=False,
                           term_aa=False,
                           **kwargs):
    """Calculate amino acid composition of a polypeptide.

    Parameters
    ----------
    sequence : str or list
        The sequence of a polypeptide or a list with a parsed sequence.
    show_unmodified_termini : bool, optional
        If :py:const:`True` then the unmodified N- and C-terminus are explicitly
        shown in the returned list. Default value is :py:const:`False`.
    term_aa : bool, optional
        If :py:const:`True` then the terminal amino acid residues are
        artificially modified with `nterm` or `cterm` modification.
        Default value is :py:const:`False`.
    labels : list, optional
        A list of allowed labels for amino acids and terminal modifications
        (default is the 20 standard amino acids, N-terminal 'H-' and C-terminal
        '-OH').

    Returns
    -------
    out : dict
        A dictionary of amino acid composition.

    Examples
    --------
    >>> amino_acid_composition('PEPTIDE')
    {'I': 1.0, 'P': 2.0, 'E': 2.0, 'T': 1.0, 'D': 1.0}
    >>> amino_acid_composition('PEPTDE', term_aa=True)
    {'ctermE': 1.0, 'E': 1.0, 'D': 1.0, 'P': 1.0, 'T': 1.0, 'ntermP': 1.0}
    >>> amino_acid_composition('PEPpTIDE',\
    labels=std_labels+['pT'])
    {'I': 1.0, 'P': 2.0, 'E': 2.0, 'D': 1.0, 'pT': 1.0}
    """
    labels = kwargs.get('labels', std_labels)

    if isinstance(sequence, str):
        parsed_sequence = parse_sequence(sequence, show_unmodified_termini,
                                         labels=labels)
    elif isinstance(sequence, list):
        parsed_sequence = sequence
    else:
        raise PyteomicsError('Unsupported type of a sequence.'
                'Must be str or list, not %s' % type(sequence))

    aa_dict = {}

    # Process terminal amino acids.    
    if term_aa:
        nterm_aa_position = 1 if is_term_mod(parsed_sequence[0]) else 0
        cterm_aa_position = (
            len(parsed_sequence) - 2 if is_term_mod(parsed_sequence[-1])
            else len(parsed_sequence) - 1)
        if len(parsed_sequence) > 1:
            aa_dict['cterm' + parsed_sequence.pop(cterm_aa_position)] = 1
        aa_dict['nterm' + parsed_sequence.pop(nterm_aa_position)] = 1

    # Process core amino acids.
    for aa in parsed_sequence:
        aa_dict[aa] = aa_dict.get(aa, 0) + 1
        
    return aa_dict

def cleave(sequence, rule, missed_cleavages=0):
    """Cleaves a polypeptide sequence using a given rule.
    
    Parameters
    ----------
    sequence : str
        The sequence of a polypeptide.
    rule : str    
        A string with a regular expression describing the C-terminal site of
        cleavage.    
    missed_cleavages : int, optional
        The maximal number of allowed missed cleavages. Defaults to 0.

    Returns
    -------
    out : list
        A list of unique (!) peptides.

    Examples
    --------
    >>> cleave('AKAKBK', expasy_rules['trypsin'], 0)
    ['AK', 'BK']
    >>> cleave('AKAKBKCK', expasy_rules['trypsin'], 2)
    ['CK', 'AKBK', 'BKCK', 'AKAK', 'AKBKCK', 'AK', 'AKAKBK', 'BK']

    """
    peptides = set()
    cleavage_sites = deque([0], maxlen=missed_cleavages+2)
    for i in chain(map(lambda x: x.end(), re.finditer(rule, sequence)),
                   [None]):
        cleavage_sites.append(i)
        for j in range(0, len(cleavage_sites)-1):
            peptides.add(sequence[cleavage_sites[j]:cleavage_sites[-1]])
    return list(peptides)

expasy_rules = {
    'arg-c':         'R',
    'asp-n':         '\w(?=D)',
    'bnps-skatole' : 'W',
    'caspase 1':     '(?<=[FWYL]\w[HAT])D(?=[^PEDQKR])',
    'caspase 2':     '(?<=DVA)D(?=[^PEDQKR])',
    'caspase 3':     '(?<=DMQ)D(?=[^PEDQKR])',
    'caspase 4':     '(?<=LEV)D(?=[^PEDQKR])',
    'caspase 5':     '(?<=[LW]EH)D',
    'caspase 6':     '(?<=VE[HI])D(?=[^PEDQKR])',
    'caspase 7':     '(?<=DEV)D(?=[^PEDQKR])',
    'caspase 8':     '(?<=[IL]ET)D(?=[^PEDQKR])',
    'caspase 9':     '(?<=LEH)D',
    'caspase 10':    '(?<=IEA)D',
    'chymotrypsin low specificity' : '([FY](?=[^P]))|(W(?=[^MP]))',
    'chymotrypsin high specificity':
        '([FLY](?=[^P]))|(W(?=[^MP]))|(M(?=[^PY]))|(H(?=[^DMPW]))',
    'clostripain':   'R',
    'cnbr':          'M',
    'enterokinase':  '(?<=[DN][DN][DN])K',
    'factor xa':     '(?<=[AFGILTVM][DE]G)R',
    'formic acid':   'D',
    'glutamyl endopeptidase': 'E',
    'granzyme b':    '(?<=IEP)D',
    'hydroxylamine': 'N(?=G)',
    'iodosobezoic acid': 'W',
    'lysc':          'K',
    'ntcb':          '\w(?=C)',
    'pepsin ph1.3':  '((?<=[^HKR][^P])[^R](?=[FLWY][^P]))|'
                     '((?<=[^HKR][^P])[FLWY](?=\w[^P]))',
    'pepsin ph2.0':  '((?<=[^HKR][^P])[^R](?=[FL][^P]))|'
                     '((?<=[^HKR][^P])[FL](?=\w[^P]))',
    'proline endopeptidase': '(?<=[HKR])P(?=[^P])',
    'proteinase k':  '[AEFILTVWY]',
    'staphylococcal peptidase i': '(?<=[^E])E',
    'thermolysin':   '[^DE](?=[AFILMV])',
    'thrombin':      '((?<=G)R(?=G))|'
                     '((?<=[AFGILTVM][AFGILTVWA]P)R(?=[^DE][^DE]))',
    'trypsin':       '([KR](?=[^P]))|((?<=W)K(?=P))|((?<=M)R(?=P))'
    }
"""
This dict contains regular expressions for cleavage rules of the most
popular proteolytic enzymes. The rules were taken from the
`PeptideCutter tool
<http://ca.expasy.org/tools/peptidecutter/peptidecutter_enzymes.html>`_
at Expasy.

.. ipython::
   
   In [2]: pprint(pyteomics.parser.expasy_rules)
"""

def isoforms(sequence, **kwargs):
    """
    Apply variable and fixed modifications to the polypeptide and return a
    :py:class:`set` of unique modified sequences.

    Parameters
    ----------
    sequence : str
        Peptide sequence to modify.
    variable_mods : dict, optional
        A dict of variable modifications in the following format:
        :py:const:`{'label1': ['X', 'Y', ...], 'label2': ['X', 'A', 'B', ...]}`

        **Note**: several variable modifications can occur on amino acids of the
        same type, but in the output each amino acid residue will be modified
        at most once.
    fixed_mods : dict, optional
        A dict of fixed modifications in the same format.
        **Note**: if a residue is affected by a fixed modification, no variable
        modifications will be applied to it.
    labels : list, optional
        A list of amino acid labels containing all the labels present in
        `sequence`. Modified entries will be added automatically.
        Defaults to :py:data:`std_labels`.

    Returns
    -------
    isoforms : set
        A set of all possible unique polypeptide sequences resulting from the
        specified modifications.
    """
    variable_mods = kwargs.get('variable_mods', {})
    fixed_mods = kwargs.get('fixed_mods', {})
    labels = kwargs.get('labels', std_labels)
    
    mods = {}
    mods.update(fixed_mods)
    mods.update(variable_mods)

    parsed = parse_sequence(sequence, 
            labels=labels+[m+aa for m in mods for aa in mods[m]])
    for cmod in fixed_mods:
        for group in parsed:
            if group in fixed_mods[cmod]:
                i = parsed.index(group)
                parsed[i] = cmod+group

    for aa in range(len(parsed)):
        if any([parsed[aa] in x for x in variable_mods.values()]):
            variable = aa
            break
    if not 'variable' in locals():
        return [''.join(parsed)]

    mod_peptides = []
    for mod in variable_mods:
        if parsed[variable] in mods[mod]:
            mod_parsed = parsed[:]
            mod_parsed[variable] = mod + parsed[variable]
            mod_peptides += [''.join(parsed[:variable+1]) + x
                    for x in isoforms(
                        ''.join(parsed[variable+1:]), **kwargs)]
            mod_peptides += [''.join(mod_parsed[:variable+1]) + x
                    for x in isoforms(
                        ''.join(mod_parsed[variable+1:]), **kwargs)]
    return set(mod_peptides)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
