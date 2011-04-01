# Licensed under the MIT license:
# http://www.opensource.org/licenses/mit-license.php 

import re
from collections import deque
from itertools import imap, chain

std_amino_acids = ['Q','W','E','R','T','Y','I','P','A','S',
                   'D','F','G','H','K','L','C','V','N','M']
std_nterm = 'H-'
std_cterm = '-OH'

std_chem_groups = std_amino_acids + [std_nterm, std_cterm]

def is_term_group(label):
    """Check if a label corresponds to a terminal group.
    Return Bool.
    """
    return label.startswith('-') or label.endswith('-')

def peptide_length(sequence, chem_groups=std_chem_groups):
    """Calculate the number of amino acid residues in a polypeptide
    written in modX notation.

    Keyword arguments:
    peptide -- a string with a polypeptide sequence, a list with a
               parsed sequence or a dict of amino acid composition;
    chem_groups -- a list of all possible amino acids and terminal
                   groups (default 20 standard amino acids, N-terminal
                   H- and C-terminal -OH).

    >>> peptide_length('PEPTIDE')
    7
    >>> peptide_length('H-PEPTIDE-OH')
    7
    """

    if isinstance(sequence, basestring) or isinstance(sequence, list):
        if isinstance(sequence, basestring):
            parsed_sequence = parse_sequence(sequence, chem_groups)
        else:
            parsed_sequence = sequence
        num_term_groups = 0
        if is_term_group(parsed_sequence[0]):
            num_term_groups += 1
        if is_term_group(parsed_sequence[-1]):
            num_term_groups += 1
        return len(parsed_sequence) - num_term_groups
    elif isinstance(sequence, dict):
        return sum([amount for aa, amount in sequence.items() 
                    if not is_term_group(aa)])

    raise Exception('Unsupported type of a sequence.')
    return None

def parse_sequence(sequence,               
                   chem_groups=std_chem_groups,
                   show_standard_term_groups=False):
    """Parse a sequence string written in modX notation into a list of
    amino acids. 

    Keyword arguments:
    sequence -- a polypeptide sequence;
    chem_groups -- a list of all possible amino acids and terminal
                   groups (default 20 standard amino acids, N-terminal
                   H- and C-terminal -OH);
    show_standard_term_groups -- if True then standard N-terminal and
                                 C-Terminal groups are always appended
                                 to the parsed sequence.

    >>> parse_sequence('PEPTIDE')
    ['P', 'E', 'P', 'T', 'I', 'D', 'E']
    >>> parse_sequence('H-PEPTIDE')
    ['P', 'E', 'P', 'T', 'I', 'D', 'E']
    >>> parse_sequence('PEPTIDE', show_standard_term_groups=True)
    ['H-', 'P', 'E', 'P', 'T', 'I', 'D', 'E', '-OH']
    >>> parse_sequence('TEpSToxM', std_chem_groups + ['pS', 'oxM'])
    ['T', 'E', 'pS', 'T', 'oxM']
    """

    stripped_sequence = str(sequence)

    # Determine the N-terminal group.
    nterm = std_nterm
    for aa in chem_groups:
        if is_term_group(aa) and stripped_sequence.startswith(aa):
            nterm = aa
            stripped_sequence = stripped_sequence[len(nterm):]
            break

    # Determine the C-terminal group.            
    cterm = std_cterm
    j = sequence.find('-')
    if j != -1:
        for aa in chem_groups:
            if is_term_group(aa) and stripped_sequence[j+1:] == aa:
                cterm = aa
                stripped_sequence = stripped_sequence[:j]
                break

    # Parse the polypeptide backbone.
    i = 0
    parsed_sequence = []
    while i<len(stripped_sequence):
        amino_acid_found = False
        for aa in chem_groups:
            if stripped_sequence.startswith(aa, i):
                parsed_sequence.append(aa)
                amino_acid_found = True
                break
        if not amino_acid_found:
            raise Exception(
                'Unknown amino acid in sequence %s at position %d: %s' % (
                    stripped_sequence, i+1, stripped_sequence[i:]))
            return []
        i = i + len(parsed_sequence[-1])

    # Append standard term groups.
    if show_standard_term_groups:
        parsed_sequence.insert(0, nterm)
        parsed_sequence.append(cterm)

    return parsed_sequence

def get_aminoacid_composition(sequence,
                              chem_groups=std_chem_groups,
                              show_standard_term_groups=False,
                              term_aa=False):
    """Calculate amino acid composition of a polypeptide.

    Keyword arguments:
    sequence -- a peptide sequence;
    chem_groups -- a list of all possible amino acids and terminal
                   groups (default 20 standard amino acids, N-terminal
                   NH2- and C-terminal -OH).
    show_standard_term_groups -- if True then standard N-terminal and
                                 C-Terminal groups are always appended
                                 to the parsed sequence.
    term_aa -- if True than terminal amino acids are treated as being
               modified with 'ntermX'/'ctermX' modifications. False by
               default.

    Returns a dictionary of amino acid content.
    >>> get_aminoacid_composition('PEPTIDE')
    {'I': 1.0, 'P': 2.0, 'E': 2.0, 'T': 1.0, 'D': 1.0}
    >>> get_aminoacid_composition('PEPTDE', term_aa=True)
    {'ctermE': 1.0, 'E': 1.0, 'D': 1.0, 'P': 1.0, 'T': 1.0, 'ntermP': 1.0}
    >>> get_aminoacid_composition('PEPpTIDE',\
                                  chem_groups=std_chem_groups+['pT'])
    {'I': 1.0, 'P': 2.0, 'E': 2.0, 'D': 1.0, 'pT': 1.0}
    """

    parsed_sequence = parse_sequence(sequence, chem_groups, 
                                   show_standard_term_groups)
    if not parsed_sequence:
        return {}

    aa_dict = {}

    # Process terminal amino acids.    
    if term_aa:
        nterm_aa_position = 1 if is_term_group(parsed_sequence[0]) else 0
        cterm_aa_position = (
            len(parsed_sequence) - 2 if is_term_group(parsed_sequence[-1])
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

    Keyword arguments:
    sequence -- a polypeptide sequence;
    rule -- a string with regular expression describing the C-terminal
            site of cleavage.
    missed_cleavages -- the maximal number of allowed missed cleavages.

    Return a list of unique (!) peptides.
    """
    peptides = set()
    cleavage_sites = deque([0], maxlen=missed_cleavages+2)
    for i in chain(imap(lambda x: x.end(), re.finditer(rule, sequence)),
                   [None]):
        cleavage_sites.append(i)
        for i in range(1, len(cleavage_sites)):
            peptides.add(sequence[cleavage_sites[0]:cleavage_sites[i]])
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

if __name__ == "__main__":
    import doctest
    doctest.testmod()
