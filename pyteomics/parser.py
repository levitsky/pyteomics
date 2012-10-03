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

  :py:func:`parse` - convert a sequence string into a list of
  amino acid residues.

  :py:func:`tostring` - convert a parsed sequence to a string.
  
  :py:func:`amino_acid_composition` - get numbers of each amino acid
  residue in a peptide.
  
  :py:func:`cleave` - cleave a polypeptide using a given rule of
  enzymatic digestion.

  :py:func:`isoforms` - generate all unique modified peptide sequences
  given the initial sequence and modifications.

Auxiliary commands
------------------

  :py:func:`peptide_length` - calculate the number of amino acid
  residues in a polypeptide.

  :py:func:`is_modX` - check if supplied code corresponds to a modX label.

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

-------------------------------------------------------------------------------

"""

#   Copyright 2012 Anton Goloborodko, Lev Levitsky
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

import re
from collections import deque
from itertools import chain, product
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
            parsed_sequence = parse(sequence, labels=labels)
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

def _split_label(label):
        if not is_modX(label):
            raise PyteomicsError('Cannot split a non-modX label: %s' % label)
        if len(label) == 1:
            return (label, )
        else:
            return (label[:-1], label[-1])

def parse(sequence,               
                   show_unmodified_termini=False, split=False,
                   allow_unknown_modifications=False,
                   **kwargs):
    """Parse a sequence string written in modX notation into a list of
    labels or (if `split` argument is :py:const:`True`) into a list of
    tuples representing amino acid residues and their modifications. 

    Parameters
    ----------
    sequence : str
        The sequence of a polypeptide.
    show_unmodified_termini : bool, optional   
        If :py:const:`True` then the unmodified N- and C-termini are explicitly
        shown in the returned list. Default value is :py:const:`False`.
    split : bool, optional
        If :py:const:`True` then the result will be a list of tuples with 1 to 4
        elements: terminal modification, modification, residue. Default value is
        :py:const:`False`.
    allow_unknown_modifications : bool, optional
        If :py:const:`True` then do not raise an exception when an unknown 
        modification of a known amino acid residue is found in the sequence.
        Default value is :py:const:`False`.
    labels : list, optional
        A list of allowed labels for amino acids, modifications and terminal
        modifications (default is the 20 standard amino acids, N-terminal H- 
        and C-terminal -OH). 

        New in ver. 1.2.2: separate labels for modifications (such as 'p' or 'ox')
        are now allowed.

    Returns
    -------
    out : list
        List of tuples with labels of modifications and amino acid residues.

    Examples
    --------
    >>> parse('PEPTIDE', split=True)
    [('P',), ('E',), ('P',), ('T',), ('I',), ('D',), ('E',)]
    >>> parse('H-PEPTIDE')
    ['P', 'E', 'P', 'T', 'I', 'D', 'E']
    >>> parse('PEPTIDE', show_unmodified_termini=True)
    ['H-', 'P', 'E', 'P', 'T', 'I', 'D', 'E', '-OH']
    >>> parse('TEpSToxM', labels=std_labels + ['pS', 'oxM'])
    ['T', 'E', 'pS', 'T', 'oxM']
    >>> parse('zPEPzTIDzE', True, True, labels=std_labels+['z'])
    [('H-', 'z', 'P'), ('E',), ('P',), ('z', 'T'), ('I',), ('D',), ('z', 'E', '-OH')]
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
            if backbone_sequence.startswith(aa, i) and is_modX(aa):
                parsed_sequence.append(aa)
                amino_acid_found = True
                break

        if not amino_acid_found:
            j = i+2
            while j <= len(backbone_sequence):
                try:
                    mod, res = _split_label(backbone_sequence[i:j])
                except PyteomicsError:
                    pass
                else:
                    if ((mod in labels and res in labels) or 
                        (allow_unknown_modifications and res in labels)):

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

    # Append terminal labels.
    if show_unmodified_termini or nterm != std_nterm:
        parsed_sequence.insert(0, nterm)
    if show_unmodified_termini or cterm != std_cterm:
        parsed_sequence.append(cterm)

    # Make a list of tuples instead of list of labels
    if split:
        if sum(1 for aa in parsed_sequence if not is_term_mod(aa)) == 1:
            result = []
            for x in parsed_sequence:
                if is_term_mod(x): result.append(x)
                else: result.extend(_split_label(x))
            return [tuple(result)]

        tuples = []
        start = 0
        if is_term_mod(parsed_sequence[0]):
            tuples.append((parsed_sequence[0],) + _split_label(parsed_sequence[1]))
            start = 2
        tuples.extend(_split_label(x) for x in parsed_sequence[start:-2])
        if is_term_mod(parsed_sequence[-1]):
            tuples.append(_split_label(parsed_sequence[-2]) + (parsed_sequence[-1],))
        else:
            tuples.extend(_split_label(x) for x in parsed_sequence[-2:])
        
        return tuples

    return parsed_sequence

def tostring(parsed_sequence, show_unmodified_termini=True):
    """Create a string from a parsed sequence.
    
    Parameters
    ----------
    parsed_sequence : iterable
        Expected to be in one of the formats returned by
        :py:func:`parse`, i.e. list of labels or list of tuples.
    show_unmodified_termini : bool, optional
        Defines the behavior towards standard terminal groups in the input.
        :py:const:`True` means that they will be preserved if present (default).
        :py:const:`False` means that they will be removed. Standard terminal
        modifications will not be added if not shown in `parsed_sequence`,
        regardless of this setting.

    Returns
    -------
    sequence : str
    """
    labels = []
    for group in parsed_sequence:
        if type(group) == str:
            if (group not in (std_cterm, std_nterm)) or show_unmodified_termini:
                labels.append(group)
        else: #treat `group` as a tuple
            group_l = list(group)
            if not show_unmodified_termini:
                if std_cterm in group_l: group_l.remove(std_cterm)
                if std_nterm in group_l: group_l.remove(std_nterm)
            labels.append(''.join(group_l))
    return ''.join(labels)
            
def amino_acid_composition(sequence,
                           show_unmodified_termini=False,
                           term_aa=False,
                           allow_unknown_modifications=False,
                           **kwargs):
    """Calculate amino acid composition of a polypeptide.

    Parameters
    ----------
    sequence : str or list
        The sequence of a polypeptide or a list with a parsed sequence.
    show_unmodified_termini : bool, optional
        If :py:const:`True` then the unmodified N- and C-terminus are explicitly
        shown in the returned dict. Default value is :py:const:`False`.
    term_aa : bool, optional
        If :py:const:`True` then the terminal amino acid residues are
        artificially modified with `nterm` or `cterm` modification.
        Default value is :py:const:`False`.
    allow_unknown_modifications : bool, optional
        If :py:const:`True` then do not raise an exception when an unknown 
        modification of a known amino acid residue is found in the sequence.
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
    {'I': 1, 'P': 2, 'E': 2, 'T': 1, 'D': 1}
    >>> amino_acid_composition('PEPTDE', term_aa=True)
    {'ctermE': 1, 'E': 1, 'D': 1, 'P': 1, 'T': 1, 'ntermP': 1}
    >>> amino_acid_composition('PEPpTIDE',\
    labels=std_labels+['pT'])
    {'I': 1, 'P': 2, 'E': 2, 'D': 1, 'pT': 1}
    """
    labels = kwargs.get('labels', std_labels)

    if isinstance(sequence, str):
        parsed_sequence = parse(sequence, show_unmodified_termini,
            allow_unknown_modifications=allow_unknown_modifications,
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
    out : set
        A set of unique (!) peptides.

    Examples
    --------
    >>> cleave('AKAKBK', expasy_rules['trypsin'], 0)
    set(['AK', 'BK'])
    >>> cleave('AKAKBKCK', expasy_rules['trypsin'], 2)
    set(['CK', 'AKBK', 'BKCK', 'AKAK', 'AKBKCK', 'AK', 'AKAKBK', 'BK'])

    """
    peptides = set()
    cleavage_sites = deque([0], maxlen=missed_cleavages+2)
    for i in chain(map(lambda x: x.end(), re.finditer(rule, sequence)),
                   [None]):
        cleavage_sites.append(i)
        for j in range(0, len(cleavage_sites)-1):
            peptides.add(sequence[cleavage_sites[j]:cleavage_sites[-1]])
    if '' in peptides:
        peptides.remove('')
    return peptides

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
"""

def isoforms(sequence, **kwargs):
    """
    Apply variable and fixed modifications to the polypeptide and yield 
    the unique modified sequences.

    Parameters
    ----------
    
    sequence : str
        Peptide sequence to modify.
    
    variable_mods : dict, optional
        A dict of variable modifications in the following format:
        :py:const:`{'label1': ['X', 'Y', ...], 'label2': ['X', 'A', 'B', ...]}`

        **Note**: several variable modifications can occur on amino acids of the
        same type, but in the output each amino acid residue will be modified
        at most once (apart from terminal modifications).
    
    fixed_mods : dict, optional
        A dict of fixed modifications in the same format.
        
        **Note**: if a residue is affected by a fixed modification, no variable
        modifications will be applied to it (apart from terminal modifications).
    
    labels : list, optional
        A list of amino acid labels containing all the labels present in
        `sequence`. Modified entries will be added automatically.
        Defaults to :py:data:`std_labels`.

    override : bool, optional
        Defines how to handle the residues that are modified in the input.
        :py:const:`False` means that they will be preserved (default).
        :py:const:`True` means they will be treated as unmodified.

        **Note**: If :py:const:`True`, then supplying fixed mods is pointless.
        
    show_unmodified_termini : bool, optional
        If :py:const:`True` then the unmodified N- and C-termini are explicitly
        shown in the returned sequences. Default value is :py:const:`False`.

    Yields
    -------
    isoform : str
        All possible unique polypeptide sequences resulting from
        the specified modifications.
    """
    def main(group): # index of the residue (capital letter) in `group`
        temp = [i for i, x in enumerate(group) if 
                (not is_term_mod(x)) and x.isupper()]
        if len(temp) != 1: raise PyteomicsError('Invalid group: %s' % group)
        return temp[0], group[temp[0]]

    def apply_mod(label, mod):
    # `label` is assumed to be a tuple (see split option of parse)
    # unmodified termini are assumed shown
    # if the modification is not applicable, `label` is returned
        group = list(label)
        m = main(group)[0]
        if m == 0 and not is_term_mod(mod): group.insert(0, mod)
        elif mod.startswith('-') and (group[-1] == std_cterm or (
            group[-1].startswith('-') and override)):
            group[-1] = mod
        elif mod.endswith('-') and (group[0] == std_nterm or (
            group[0].endswith('-') and override)):
            group[0] = mod
        elif not is_term_mod(mod):
            if not group[m-1].endswith('-'):
                if override: group[m-1] = mod
            else: group.insert(m, mod)
        return tuple(group)            

    variable_mods = kwargs.get('variable_mods', {})
    fixed_mods = kwargs.get('fixed_mods', {})
    labels = kwargs.get('labels', std_labels)
    length = peptide_length(sequence, labels=labels)
    parsed = parse(sequence, True, True,
            labels=labels+list(fixed_mods.keys()))
    override = kwargs.get('override', False)
    show_unmodified_termini = kwargs.get('show_unmodified_termini', False)

    # Apply fixed modifications
    for cmod in fixed_mods:
        for i, group in enumerate(parsed):
            if main(group)[1] in fixed_mods[cmod]:
                parsed[i] = apply_mod(group, cmod)

    # Create a list of possible states for each group
    # Start with N-terminal mods and regular mods on the N-terminal residue
    second = set(apply_mod(parsed[0], m) for m, r in variable_mods.items()
            if main(parsed[0])[1] in r and not is_term_mod(m)).union([parsed[0]])
    first = chain((apply_mod(group, mod) for group in second 
            for mod, res in variable_mods.items()
        if (mod.endswith('-') or mod.startswith('-') and len(parsed) == 1)
        and main(group)[1] in res), second)
    states = [set(first).union([parsed[0]])]
    # Continue with regular mods
    states.extend([set([group]).union(set([apply_mod(group, mod) for mod in variable_mods if
        main(group)[1] in variable_mods[mod] and not is_term_mod(mod)])) for group in parsed[1:-1]])
    # Finally add C-terminal mods and regular mods on the C-terminal residue
    if len(parsed) > 1:
        second = set(apply_mod(parsed[-1], m) for m, r in variable_mods.items()
                if main(parsed[-1])[1] in r and not is_term_mod(m)).union([parsed[-1]])
        first = chain((apply_mod(group, mod) for group in second 
                for mod, res in variable_mods.items()
            if mod.startswith('-') and main(group)[1] in res), second)
        states.append(set(first).union([parsed[-1]]))

    return (tostring(x, show_unmodified_termini) for x in product(*states))
    

if __name__ == "__main__":
    import doctest
    doctest.testmod()
