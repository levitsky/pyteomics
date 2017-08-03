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

Besides the sequence of amino acid residues, modX has a rule to specify
terminal modifications of a polypeptide. Such a label should start or
end with a hyphen. The default N-terminal amine group and C-terminal
carboxyl group may not be shown explicitly.

Therefore, valid examples of peptide sequences in modX are: "GAGA",
"H-PEPTIDE-OH", "H-TEST-NH2". It is not recommmended to specify only one
terminal group.

Operations on polypeptide sequences
-----------------------------------

  :py:func:`parse` - convert a sequence string into a list of
  amino acid residues.

  :py:func:`tostring` - convert a parsed sequence to a string.

  :py:func:`amino_acid_composition` - get numbers of each amino acid
  residue in a peptide.

  :py:func:`cleave` - cleave a polypeptide using a given rule of
  enzymatic digestion.

  :py:func:`num_sites` - count the number of cleavage sites in a sequence.

  :py:func:`isoforms` - generate all unique modified peptide sequences
  given the initial sequence and modifications.

Auxiliary commands
------------------

  :py:func:`coverage` - calculate the sequence coverage of a protein by peptides.

  :py:func:`length` - calculate the number of amino acid
  residues in a polypeptide.

  :py:func:`valid` - check if a sequence can be parsed successfully.

  :py:func:`fast_valid` - check if a sequence contains of known one-letter
  codes.

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

  :py:data:`expasy_rules` - a dict with the regular expressions of
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
import itertools as it
from .auxiliary import PyteomicsError, memoize, BasicComposition

std_amino_acids = ['Q', 'W', 'E', 'R', 'T', 'Y', 'I', 'P', 'A', 'S',
                   'D', 'F', 'G', 'H', 'K', 'L', 'C', 'V', 'N', 'M']
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
    return label[0] == '-' or label[-1] == '-'

def match_modX(label):
    """Check if `label` is a valid 'modX' label.

    Parameters
    ----------
    label : str

    Returns
    -------
    out : re.match or None
    """
    return re.match(_modX_split, label)

def is_modX(label):
    """Check if `label` is a valid 'modX' label.

    Parameters
    ----------
    label : str

    Returns
    -------
    out : bool
    """
    return bool(match_modX(label))

def length(sequence, **kwargs):
    """Calculate the number of amino acid residues in a polypeptide
    written in modX notation.

    Parameters
    ----------
    sequence : str or list or dict
        A string with a polypeptide sequence, a list with a parsed sequence or
        a dict of amino acid composition.
    labels : list, optional
        A list of allowed labels for amino acids and terminal modifications.

    Examples
    --------
    >>> length('PEPTIDE')
    7
    >>> length('H-PEPTIDE-OH')
    7
    """
    if not sequence: return 0

    if isinstance(sequence, str) or isinstance(sequence, list):
        if isinstance(sequence, str):
            parsed_sequence = parse(sequence, **kwargs)
        else:
            parsed_sequence = sequence
        num_term_groups = 0
        if is_term_mod(parsed_sequence[0]):
            num_term_groups += 1
        if is_term_mod(parsed_sequence[-1]):
            num_term_groups += 1
        return len(parsed_sequence) - num_term_groups
    elif isinstance(sequence, dict):
        return sum(amount for aa, amount in sequence.items()
                    if not is_term_mod(aa))

    raise PyteomicsError('Unsupported type of sequence.')

def _split_label(label):
    try:
        mod, X = match_modX(label).groups()
    except AttributeError:
        raise PyteomicsError('Cannot split a non-modX label: %s' % label)
    if not mod:
        return (X,)
    else:
        return mod, X

_modX_sequence = re.compile(r'^([^-]+-)?((?:[a-z]*[A-Z])+)(-[^-]+)?$')
_modX_group = re.compile(r'[a-z]*[A-Z]')
_modX_split = re.compile(r'([a-z]*)([A-Z])')

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
        This also includes terminal groups.
        Default value is :py:const:`False`.

        .. note::
            Since version 2.5, this parameter has effect only if `labels`
            are provided.
    labels : container, optional
        A container of allowed labels for amino acids,
        modifications and terminal modifications.
        If not provided, no checks will be done.
        Separate labels for modifications (such as 'p' or 'ox')
        can be supplied, which means they are applicable to all residues.

        .. warning::
            If `show_unmodified_termini` is set to :py:const:`True`, standard
            terminal groups need to be present in `labels`.

        .. warning::
            Avoid using sequences with only one terminal group, as they are
            ambiguous. If you provide one, `labels` (or :py:const:`std_labels`)
            will be used to resolve the ambiguity.

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
    sequence = str(sequence)

    try:
        n, body, c = re.match(_modX_sequence, sequence).groups()
    except AttributeError:
        raise PyteomicsError('Not a valid modX sequence: ' + sequence)

    # Check for allowed labels, if they were explicitly given
    labels = kwargs.get('labels')
    # labels help save the day when only one terminal group is given
    if c is None and n is not None:
        if labels is None:
            labels = std_labels
        # we can try to resolve the ambiguity
        if n != std_nterm and n not in labels:
            # n is the body then
            c = '-' + body
            body = n[:-1]
            n = None

    # Actual parsing
    if split:
        parsed_sequence = [g if g[0] else (g[1],) for g in re.findall(
            _modX_split, body)]
    else:
        parsed_sequence = re.findall(_modX_group, body)
    nterm, cterm = (n or std_nterm), (c or std_cterm)

    # Check against `labels` if given
    if labels is not None:
        labels = set(labels)
        for term, std_term in zip([n, c], [std_nterm, std_cterm]):
            if term and term not in labels and not allow_unknown_modifications:
                raise PyteomicsError(
                            'Unknown label: {}'.format(term))
        for group in parsed_sequence:
            if split:
                mod, X = group if len(group) == 2 else ('', group[0])
            else:
                mod, X = re.match(_modX_split, group).groups()
            if ((not mod) and X not in labels) or not ((mod+X in labels) or (
                X in labels and (
                    mod in labels or allow_unknown_modifications))):
                raise PyteomicsError(
                        'Unknown label: {}'.format(group))

    # Append terminal labels
    if show_unmodified_termini or nterm != std_nterm:
        if split:
            parsed_sequence[0] = (nterm,) + parsed_sequence[0]
        else:
            parsed_sequence.insert(0, nterm)
    if show_unmodified_termini or cterm != std_cterm:
        if split:
            parsed_sequence[-1] = parsed_sequence[-1] + (cterm,)
        else:
            parsed_sequence.append(cterm)


    return parsed_sequence

def valid(*args, **kwargs):
    """Try to parse sequence and catch the exceptions.
    All parameters are passed to :py:func:`parse`.

    Returns
    -------
    out : bool
        :py:const:`True` if the sequence was parsed successfully, and
        :py:const:`False` otherwise.
    """
    try:
        parse(*args, **kwargs)
    except PyteomicsError:
        return False
    return True

def fast_valid(sequence, labels=set(std_labels)):
    """Iterate over `sequence` and check if all items are in `labels`.
    With strings, this only works as expected on sequences without
    modifications or terminal groups.

    Parameters
    ----------
    sequence : iterable (expectedly, str)
        The sequence to check. A valid sequence would be a string of
        labels, all present in `labels`.
    labels : iterable, optional
        An iterable of known labels.

    Returns
    -------
    out : bool
    """
    return set(sequence).issubset(labels)

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
        groups will not be added if not shown in `parsed_sequence`,
        regardless of this setting.

    Returns
    -------
    sequence : str
    """
    parsed_sequence = list(parsed_sequence)
    labels = []
    nterm = parsed_sequence[0]
    cterm = parsed_sequence[-1]

    if isinstance(nterm, str):
        if nterm != std_nterm or show_unmodified_termini:
            labels.append(nterm)
        labels.extend(parsed_sequence[1:-1])
        if len(parsed_sequence) > 1 and (cterm != std_cterm or show_unmodified_termini):
            labels.append(cterm)
    else:
        if len(parsed_sequence) == 1:
            g = nterm
            if nterm[0] == std_nterm and not show_unmodified_termini:
                g = g[1:]
            if nterm[-1] == std_cterm and not show_unmodified_termini:
                g = g[:-1]
            return ''.join(g)
        if nterm[0] != std_nterm or show_unmodified_termini:
            labels.append(''.join(nterm))
        else:
            labels.append(''.join(nterm[1:]))
        labels.extend(''.join(g) for g in parsed_sequence[1:-1])
        if len(parsed_sequence) > 1:
            if cterm[-1] != std_cterm or show_unmodified_termini:
                labels.append(''.join(cterm))
            else:
                labels.append(''.join(cterm[:-1]))
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
        A list of allowed labels for amino acids and terminal modifications.

    Returns
    -------
    out : dict
        A dictionary of amino acid composition.

    Examples
    --------
    >>> amino_acid_composition('PEPTIDE') == \
    {'I': 1, 'P': 2, 'E': 2, 'T': 1, 'D': 1}
    True
    >>> amino_acid_composition('PEPTDE', term_aa=True) == \
    {'ctermE': 1, 'E': 1, 'D': 1, 'P': 1, 'T': 1, 'ntermP': 1}
    True
    >>> amino_acid_composition('PEPpTIDE', labels=std_labels+['pT']) == \
    {'I': 1, 'P': 2, 'E': 2, 'D': 1, 'pT': 1}
    True
    """
    labels = kwargs.get('labels')

    if isinstance(sequence, str):
        parsed_sequence = parse(sequence, show_unmodified_termini,
            allow_unknown_modifications=allow_unknown_modifications,
            labels=labels)
    elif isinstance(sequence, list):
        if sequence and isinstance(sequence[0], tuple):
            parsed_sequence = parse(tostring(sequence, True),
                show_unmodified_termini,
                allow_unknown_modifications=allow_unknown_modifications,
                labels=labels)
        else:
            parsed_sequence = sequence
    else:
        raise PyteomicsError('Unsupported type of a sequence.'
                'Must be str or list, not %s' % type(sequence))

    aa_dict = BasicComposition()

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
        aa_dict[aa] += 1

    return aa_dict

@memoize()
def cleave(sequence, rule, missed_cleavages=0, min_length=None):
    """Cleaves a polypeptide sequence using a given rule.

    Parameters
    ----------
    sequence : str
        The sequence of a polypeptide.

        .. note::
            The sequence is expected to be in one-letter uppercase notation.
            Otherwise, some of the cleavage rules in :py:data:`expasy_rules`
            will not work as expected.

    rule : str or compiled regex
        A `regular expression <https://docs.python.org/library/re.html#regular-expression-syntax>`_
        describing the site of cleavage. It is recommended
        to design the regex so that it matches only the residue whose C-terminal
        bond is to be cleaved. All additional requirements should be specified
        using `lookaround assertions
        <http://www.regular-expressions.info/lookaround.html>`_.
        :py:data:`expasy_rules` contains cleavage rules for popular cleavage agents.
    missed_cleavages : int, optional
        Maximum number of allowed missed cleavages. Defaults to 0.
    min_length : int or None, optional
        Minimum peptide length. Defaults to :py:const:`None`.

        .. note ::
            This checks for string length, which is only correct for one-letter
            notation and not for full *modX*. Use :py:func:`length` manually if
            you know what you are doing and apply :py:func:`cleave` to *modX*
            sequences.

    Returns
    -------
    out : set
        A set of unique (!) peptides.

    Examples
    --------
    >>> cleave('AKAKBK', expasy_rules['trypsin'], 0) == {'AK', 'BK'}
    True
    >>> cleave('GKGKYKCK', expasy_rules['trypsin'], 2) == \
    {'CK', 'GKYK', 'YKCK', 'GKGK', 'GKYKCK', 'GK', 'GKGKYK', 'YK'}
    True

    """
    return set(_cleave(sequence, rule, missed_cleavages, min_length))

def _cleave(sequence, rule, missed_cleavages=0, min_length=None):
    """Like :py:func:`cleave`, but the result is a list. Refer to
    :py:func:`cleave` for explanation of parameters.
    """
    peptides = []
    ml = missed_cleavages+2
    trange = range(ml)
    cleavage_sites = deque([0], maxlen=ml)
    cl = 1
    for i in it.chain([x.end() for x in re.finditer(rule, sequence)],
                   [None]):
        cleavage_sites.append(i)
        if cl < ml:
            cl += 1
        for j in trange[:cl-1]:
            seq = sequence[cleavage_sites[j]:cleavage_sites[-1]]
            if seq:
                if min_length is None or len(seq) >= min_length:
                    peptides.append(seq)
    return peptides

def num_sites(sequence, rule, **kwargs):
    """Count the number of sites where `sequence` can be cleaved using
    the given `rule` (e.g. number of miscleavages for a peptide).

    Parameters
    ----------
    sequence : str
        The sequence of a polypeptide.
    rule : str or compiled regex
        A regular expression describing the site of cleavage. It is recommended
        to design the regex so that it matches only the residue whose C-terminal
        bond is to be cleaved. All additional requirements should be specified
        using `lookaround assertions
        <http://www.regular-expressions.info/lookaround.html>`_.
    labels : list, optional
        A list of allowed labels for amino acids and terminal modifications.

    Returns
    -------
    out : int
        Number of cleavage sites.
    """
    return len(_cleave(sequence, rule, **kwargs)) - 1

expasy_rules = {
    'arg-c':         r'R',
    'asp-n':         r'\w(?=D)',
    'bnps-skatole' : r'W',
    'caspase 1':     r'(?<=[FWYL]\w[HAT])D(?=[^PEDQKR])',
    'caspase 2':     r'(?<=DVA)D(?=[^PEDQKR])',
    'caspase 3':     r'(?<=DMQ)D(?=[^PEDQKR])',
    'caspase 4':     r'(?<=LEV)D(?=[^PEDQKR])',
    'caspase 5':     r'(?<=[LW]EH)D',
    'caspase 6':     r'(?<=VE[HI])D(?=[^PEDQKR])',
    'caspase 7':     r'(?<=DEV)D(?=[^PEDQKR])',
    'caspase 8':     r'(?<=[IL]ET)D(?=[^PEDQKR])',
    'caspase 9':     r'(?<=LEH)D',
    'caspase 10':    r'(?<=IEA)D',
    'chymotrypsin high specificity' : r'([FY](?=[^P]))|(W(?=[^MP]))',
    'chymotrypsin low specificity':
        r'([FLY](?=[^P]))|(W(?=[^MP]))|(M(?=[^PY]))|(H(?=[^DMPW]))',
    'clostripain':   r'R',
    'cnbr':          r'M',
    'enterokinase':  r'(?<=[DE]{3})K',
    'factor xa':     r'(?<=[AFGILTVM][DE]G)R',
    'formic acid':   r'D',
    'glutamyl endopeptidase': r'E',
    'granzyme b':    r'(?<=IEP)D',
    'hydroxylamine': r'N(?=G)',
    'iodosobenzoic acid': r'W',
    'lysc':          r'K',
    'ntcb':          r'\w(?=C)',
    'pepsin ph1.3':  r'((?<=[^HKR][^P])[^R](?=[FLWY][^P]))|'
                     r'((?<=[^HKR][^P])[FLWY](?=\w[^P]))',
    'pepsin ph2.0':  r'((?<=[^HKR][^P])[^R](?=[FL][^P]))|'
                     r'((?<=[^HKR][^P])[FL](?=\w[^P]))',
    'proline endopeptidase': r'(?<=[HKR])P(?=[^P])',
    'proteinase k':  r'[AEFILTVWY]',
    'staphylococcal peptidase i': r'(?<=[^E])E',
    'thermolysin':   r'[^DE](?=[AFILMV])',
    'thrombin':      r'((?<=G)R(?=G))|'
                     r'((?<=[AFGILTVM][AFGILTVWA]P)R(?=[^DE][^DE]))',
    'trypsin':       r'([KR](?=[^P]))|((?<=W)K(?=P))|((?<=M)R(?=P))'
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

        Keys in the dict are modification labels (terminal modifications allowed).
        Values are iterables of residue labels (one letter each) or
        :py:const:`True`. If a value for a modification is :py:const:`True`,
        it is applicable to any residue (useful for terminal modifications).
        You can use values such as 'ntermX' or 'ctermY' to specify that a
        mdofication only occurs when the residue is in the terminal position.
        This is *not needed* for terminal modifications.

        .. note:: Several variable modifications can occur on amino acids of the
                  same type, but in the output each amino acid residue will be
                  modified at most once (apart from terminal modifications).

    fixed_mods : dict, optional
        A dict of fixed modifications in the same format.

        **Note**: if a residue is affected by a fixed modification, no variable
        modifications will be applied to it (apart from terminal modifications).

    labels : list, optional
        A list of amino acid labels containing all the labels present in
        `sequence`. Modified entries will be added automatically.
        Defaults to :py:data:`std_labels`.
        Not required since version 2.5.

    max_mods : int or None, optional
        Number of modifications that can occur simultaneously on a peptide,
        excluding fixed modifications. If :py:const:`None` or if ``max_mods``
        is greater than the number of modification sites, all possible
        isoforms are generated. Default is :py:const:`None`.

    override : bool, optional
        Defines how to handle the residues that are modified in the input.
        :py:const:`False` means that they will be preserved (default).
        :py:const:`True` means they will be treated as unmodified.

    show_unmodified_termini : bool, optional
        If :py:const:`True` then the unmodified N- and C-termini are explicitly
        shown in the returned sequences. Default value is :py:const:`False`.

    format : str, optional
        If :py:const:`'str'` (default), an iterator over sequences is returned.
        If :py:const:`'split'`, the iterator will yield results in the same
        format as :py:func:`parse` with the 'split' option, with unmodified
        terminal groups shown.

    Returns
    -------

    out : iterator over strings or lists
        All possible unique polypeptide sequences resulting from
        the specified modifications are yielded obe by one.
    """
    def main(group): # index of the residue (capital letter) in `group`
        if group[-1][0] == '-':
            i = -2
        else:
            i = -1
        return len(group) + i, group[i]

    def apply_mod(label, mod):
    # `label` is assumed to be a tuple (see split option of parse)
    # unmodified termini are assumed shown
    # if the modification is not applicable, `label` is returned
        group = list(label)
        m = main(group)[0]
        if m == 0 and not is_term_mod(mod):
            group.insert(0, mod)
        elif mod[0] == '-' and (group[-1] == std_cterm or (
            group[-1][0] == '-' and override)):
            group[-1] = mod
        elif mod[-1] == '-' and (group[0] == std_nterm or (
            group[0][-1] == '-' and override)):
            group[0] = mod
        elif not is_term_mod(mod):
            if m and not group[m-1][-1] == '-':
                if override:
                    group[m-1] = mod
            else:
                group.insert(m, mod)
        return tuple(group)

    variable_mods = kwargs.get('variable_mods', {})
    fixed_mods = kwargs.get('fixed_mods', {})
    parse_kw = {}
    if 'labels' in kwargs:
        parse_kw['labels'] = list(kwargs['labels']) + list(fixed_mods)
    parsed = parse(sequence, True, True, **parse_kw)
    override = kwargs.get('override', False)
    show_unmodified_termini = kwargs.get('show_unmodified_termini', False)
    max_mods = kwargs.get('max_mods')
    format_ = kwargs.get('format', 'str')

    # Apply fixed modifications
    for cmod in fixed_mods:
        for i, group in enumerate(parsed):
            if fixed_mods[cmod] == True or main(group)[1] in fixed_mods[cmod]:
                parsed[i] = apply_mod(group, cmod)

    # Create a list of possible states for each group
    # Start with N-terminal mods and regular mods on the N-terminal residue
    second = set(apply_mod(parsed[0], m) for m, r in variable_mods.items()
                if (r == True or
                    main(parsed[0])[1] in r or
                    'nterm' + main(parsed[0])[1] in r or
                    (len(parsed) == 1 and 'cterm' + main(parsed[0])[1] in r))
                and not is_term_mod(m)
                ).union([parsed[0]])
    first = it.chain((apply_mod(group, mod) for group in second
            for mod, res in variable_mods.items()
        if (mod.endswith('-') or (mod.startswith('-') and len(parsed) == 1))
        and (res == True or main(group)[1] in res)), second)
    states = [[parsed[0]] + list(set(first).difference({parsed[0]}))]
    # Continue with regular mods
    states.extend([group] + list(set(apply_mod(group, mod)
            for mod in variable_mods if (
                variable_mods[mod] == True or
                group[-1] in variable_mods[mod]) and not is_term_mod(mod)
                ).difference({group}))
        for group in parsed[1:-1])
    # Finally add C-terminal mods and regular mods on the C-terminal residue
    if len(parsed) > 1:
        second = set(apply_mod(parsed[-1], m) for m, r in variable_mods.items()
                    if (r == True or
                        main(parsed[-1])[1] in r or
                        'cterm' + main(parsed[-1])[1] in r)
                    and not is_term_mod(m)
                ).union((parsed[-1],))
        first = it.chain((apply_mod(group, mod) for group in second
                for mod, res in variable_mods.items()
            if mod.startswith('-') and (
                res == True or main(group)[1] in res)), second)
        states.append([parsed[-1]] + list(set(first).difference({parsed[-1]})))

    sites = [s for s in enumerate(states) if len(s[1]) > 1]
    if max_mods is None or max_mods > len(sites):
        possible_states = it.product(*states)
    else:
        def state_lists():
            for m in range(max_mods+1):
                for comb in it.combinations(sites, m):
                    skel = [[s[0]] for s in states]
                    for i, e in comb:
                        skel[i] = e[1:]
                    yield skel
        possible_states = it.chain.from_iterable(
                it.product(*skel) for skel in state_lists())

    if format_ == 'split':
        def strip_std_terms():
            for ps in possible_states:
                ps = list(ps)
                if not show_unmodified_termini:
                    if ps[0][0] == std_nterm:
                        ps[0] = ps[0][1:]
                    if ps[-1][-1] == std_cterm:
                        ps[-1] = ps[-1][:-1]
                yield ps
        return strip_std_terms()
    elif format_ == 'str':
        return (tostring(form, show_unmodified_termini)
            for form in possible_states)
    else:
        raise PyteomicsError('Unsupported value of "format": {}'.format(format_))

def coverage(protein, peptides):
    """Calculate how much of `protein` is covered by `peptides`.
    Peptides can overlap. If a peptide is found multiple times in `protein`,
    it contributes more to the overall coverage.

    Requires :py:mod:`numpy`.

    .. note::
        Modifications and terminal groups are discarded.

    Parameters
    ----------
    protein : str
        A protein sequence.
    peptides : iterable
        An iterable of peptide sequences.

    Returns
    -------
    out : float
        The sequence coverage, between 0 and 1.

    Examples
    --------
    >>> coverage('PEPTIDES'*100, ['PEP', 'EPT'])
    0.5
    """
    import numpy as np
    protein = re.sub(r'[^A-Z]', '', protein)
    mask = np.zeros(len(protein), dtype=np.int8)
    for peptide in peptides:
        indices = [m.start() for m in re.finditer(
            '(?={})'.format(re.sub(r'[^A-Z]', '', peptide)), protein)]
        for i in indices:
            mask[i:i+len(peptide)] = 1
    return mask.sum(dtype=float) / mask.size


if __name__ == "__main__":
    import doctest
    doctest.testmod()
