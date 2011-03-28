import re
from collections import deque
from itertools import imap, chain

def read_fasta(fasta_file):
    """Reads fasta file and iterates over proteins.

    Keyword arguments:
    fasta_file -- a file object with fasta database
    """
    accumulated_strings = []
    for string in fasta_file:
        if string.startswith('>') and accumulated_strings:
            description = accumulated_strings[0][1:].strip()
            sequence = "".join([i.strip() for i in accumulated_strings[1:]])
            yield (description, sequence)
            accumulated_strings = []
        accumulated_strings.append(string)

def cleave(sequence, rule, missed_cleavages=0):
    """Cleaves a polypeptide sequence using given rule.

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
    
