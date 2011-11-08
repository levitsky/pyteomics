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

Operations on polypeptide sequences:
------------------------------------

  :py:func:`parse_sequence` - convert a sequence string into a list of
  amino acid residues.
  
  :py:func:`amino_acid_composition` - get numbers of each amino acid
  residue in a peptide.
  
  :py:func:`cleave` - cleave a polypeptide using a given rule of
  enzymatic digestion.

Auxiliary commands:
-------------------

  :py:func:`peptide_length` - calculate the number of amino acid
  residues in a polypeptide.

  :py:func:`is_term_mod` - check if supplied code corresponds to a
  terminal modification.

Data:
-----

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
from itertools import imap, chain
from auxiliary import PyteomicsError

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
    """Check if a label corresponds to a terminal modification.
    Return Bool.
    """
    return label.startswith('-') or label.endswith('-')

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

    if isinstance(sequence, basestring) or isinstance(sequence, list):
        if isinstance(sequence, basestring):
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
        return sum([amount for aa, amount in sequence.items() 
                    if not is_term_mod(aa)])

    raise PyteomicsError('Unsupported type of a sequence.')
    return None

def parse_sequence(sequence,               
                   show_unmodified_termini=False,
                   **kwargs):
    """Parse a sequence string written in modX notation into a list of
    amino acid residues and terminal modifications. 

    Parameters
    ----------
    sequence : str
        The sequence of a polypeptide.
    show_unmodified_termini : bool    
        If True then the unmodified N- and C-termini are explicitly shown in
        the returned list.
    labels : list, optional
        A list of allowed labels for amino acids and terminal modifications
        (default is the 20 standard amino acids, N-terminal H- and C-terminal
        -OH).

    Returns
    -------
    out : list
        List of labels of terminal modifications and amino acid residues.

    Examples
    --------
    >>> parse_sequence('PEPTIDE')
    ['P', 'E', 'P', 'T', 'I', 'D', 'E']
    >>> parse_sequence('H-PEPTIDE')
    ['P', 'E', 'P', 'T', 'I', 'D', 'E']
    >>> parse_sequence('PEPTIDE', show_unmodified_termini=True)
    ['H-', 'P', 'E', 'P', 'T', 'I', 'D', 'E', '-OH']
    >>> parse_sequence('TEpSToxM', labels=std_labels + ['pS', 'oxM'])
    ['T', 'E', 'pS', 'T', 'oxM']
    """

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
            if backbone_sequence.startswith(aa, i):
                parsed_sequence.append(aa)
                amino_acid_found = True
                break
        if not amino_acid_found:
            raise PyteomicsError(
                'Unknown amino acid in sequence %s at position %d: %s' % (
                    backbone_sequence, i+1, backbone_sequence[i:]))
            return []
        i = i + len(parsed_sequence[-1])

    # Append labels of unmodified termini.
    if show_unmodified_termini:
        parsed_sequence.insert(0, nterm)
        parsed_sequence.append(cterm)

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
    show_unmodified_termini : bool    
        If True then the unmodified N- and C-terminus are explicitly shown in
        the returned list.
    term_aa : bool
        If True then the terminal amino acid residues are artificially
        modified with "nterm" or "cterm" modification.     
    labels : list, optional
        A list of allowed labels for amino acids and terminal modifications
        (default is the 20 standard amino acids, N-terminal H- and C-terminal
        -OH).

    Returns
    -------
    out : dict or None
        a dictionary of amino acid content. Returns None if `sequence` is not
        of supported type.

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

    if isinstance(sequence, basestring):
        parsed_sequence = parse_sequence(sequence, show_unmodified_termini,
                                         labels=labels)
    elif isinstance(sequence, list):
        parsed_sequence = sequence
    else:
        raise PyteomicsError('Unsupported type of a sequence.')

    aa_dict = {}

    # Process terminal amino acids.    
    if term_aa:
        nterm_aa_position = 1 if is_term_mod(parsed_sequence[0]) else 0
        cterm_aa_position = (
            len(parsed_sequence) - 2 if is_term_mod(parsed_sequence[-1])
            else len(parsed_sequence) - 1)
        if len(parsed_sequence) > 1:
            aa_dict['cterm' + parsed_sequence.pop(cterm_aa_position)] = 1.0
        aa_dict['nterm' + parsed_sequence.pop(nterm_aa_position)] = 1.0

    # Process core amino acids.
    for aa in parsed_sequence:
        aa_dict[aa] = aa_dict.get(aa, 0.0) + 1.0
        
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
    missed_cleavages : int
        The maximal number of allowed missed cleavages.

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
    for i in chain(imap(lambda x: x.end(), re.finditer(rule, sequence)),
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

if __name__ == "__main__":
    import doctest
    doctest.testmod()
