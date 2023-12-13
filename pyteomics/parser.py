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

  :py:func:`to_string` - convert a parsed sequence to a string.

  :py:func:`to_proforma` - convert a (parsed) *modX* sequence to ProForma.

  :py:func:`amino_acid_composition` - get numbers of each amino acid
  residue in a peptide.

  :py:func:`cleave`, :py:func:`icleave`, :py:func:`xcleave` - cleave a polypeptide using a given rule of
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

  :py:data:`expasy_rules` and :py:data:`psims_rules` - two dicts with the regular expressions of
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
import warnings
from .auxiliary import PyteomicsError, memoize, BasicComposition, cvstr, cvquery


std_amino_acids = ['Q', 'W', 'E', 'R', 'T', 'Y', 'I', 'P', 'A', 'S',
                   'D', 'F', 'G', 'H', 'K', 'L', 'C', 'V', 'N', 'M']
"""modX labels for the 20 standard amino acids."""

std_nterm = 'H-'
"""modX label for the unmodified N-terminus."""

std_cterm = '-OH'
"""modX label for the unmodified C-terminus."""

std_labels = std_amino_acids + [std_nterm, std_cterm]
"""modX labels for the standard amino acids and unmodified termini."""

_nterm_mod = r'[^-]+-$'
_cterm_mod = r'-[^-]+$'


def is_term_mod(label):
    """Check if `label` corresponds to a terminal modification.

    Parameters
    ----------
    label : str

    Returns
    -------
    out : bool

    Examples
    --------
    >>> is_term_mod('A')
    False
    >>> is_term_mod('Ac-')
    True
    >>> is_term_mod('-customGroup')
    True
    >>> is_term_mod('this-group-')
    False
    >>> is_term_mod('-')
    False
    """
    return (re.match(_nterm_mod, label) or re.match(_cterm_mod, label)) is not None


def match_modX(label):
    """Check if `label` is a valid 'modX' label.

    Parameters
    ----------
    label : str

    Returns
    -------
    out : re.match or None
    """
    return re.match(_modX_single, label)


def is_modX(label):
    """Check if `label` is a valid 'modX' label.

    Parameters
    ----------
    label : str

    Returns
    -------
    out : bool

    Examples
    --------
    >>> is_modX('M')
    True
    >>> is_modX('oxM')
    True
    >>> is_modX('oxMet')
    False
    >>> is_modX('160C')
    True
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

    Returns
    -------
    out : int

    Examples
    --------
    >>> length('PEPTIDE')
    7
    >>> length('H-PEPTIDE-OH')
    7
    """
    if not sequence:
        return 0

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
        return sum(amount for aa, amount in sequence.items() if not is_term_mod(aa))

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


_modX_sequence = re.compile(r'^([^-]+-)?((?:[^A-Z-]*[A-Z])+)(-[^-]+)?$')
_modX_group = re.compile(r'[^A-Z-]*[A-Z]')
_modX_split = re.compile(r'([^A-Z-]*)([A-Z])')
_modX_single = re.compile(r'^([^A-Z-]*)([A-Z])$')


def parse(sequence, show_unmodified_termini=False, split=False, allow_unknown_modifications=False, **kwargs):
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
    >>> parse('Pmod1EPTIDE')
    ['P', 'mod1E', 'P', 'T', 'I', 'D', 'E']
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
                raise PyteomicsError('Unknown label: {}'.format(term))
        for group in parsed_sequence:
            if split:
                mod, X = group if len(group) == 2 else ('', group[0])
            else:
                mod, X = re.match(_modX_split, group).groups()
            if ((not mod) and X not in labels) or not ((mod + X in labels) or (
                X in labels and (
                    mod in labels or allow_unknown_modifications))):
                raise PyteomicsError('Unknown label: {}'.format(group))

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


def to_string(parsed_sequence, show_unmodified_termini=True):
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


tostring = to_string


def to_proforma(sequence, **kwargs):
    """Converts a (parsed) *modX* sequence to a basic ProForma string.
    Modifications are represented as masses, if those are given in :arg:`aa_mass`,
    as chemical formulas (via :arg:`aa_comp`) or as names (using :arg:`mod_names`).

    Parameters
    ----------
    sequence : str or list
        A *modX* sequence, possibly in the parsed form.
    aa_mass : dict, keyword only, optional
        Used to render modifications as mass shifts.
    aa_comp : dict, keyword only, optional
        Used to render modifications as chemical formulas.
    mod_names : dict or callable, keyword only, optional
        Used to get the rendered name of modification from the mod label.
    prefix : str, keyword only, optional
        Prepend all modification names with the given prefix.

    Returns
    -------
    out : str
        A ProForma sequence.

    Examples
    --------
    >>> to_proforma('PEPTIDE')
    'PEPTIDE'
    >>> to_proforma('Ac-oxMYPEPTIDE-OH', aa_mass={'Ac-': 42.010565}, mod_names={'ox': 'Oxidation'}, prefix='U:')
    '[+42.0106]-M[U:Oxidation]YPEPTIDE'
    >>> to_proforma('oxidationMYPEPTIDE')  # last fallback is to just capitalize the label
    'M[Oxidation]YPEPTIDE'
    """
    from . import proforma
    from .mass.mass import std_aa_mass, std_aa_comp

    if isinstance(sequence, str):
        return to_proforma(parse(sequence), **kwargs)

    aa_mass = kwargs.get('aa_mass', std_aa_mass)
    aa_comp = kwargs.get('aa_comp', std_aa_comp)
    mod_names = kwargs.get('mod_names', {})
    prefix = kwargs.get('prefix', '')

    if isinstance(mod_names, dict):
        get_name = mod_names.get
    else:
        get_name = mod_names

    def get_tag(label):
        if label in aa_mass:
            return [proforma.MassModification(aa_mass[label])]
        if label in aa_comp:
            return [proforma.FormulaModification(''.join('{}{}'.format(k, v if v not in {0, 1} else '') for k, v in aa_comp[label].items()))]
        name = get_name(label)
        if not name:
            warnings.warn("Unable to resolve label `{}`. "
                "The ProForma string may be invalid. Specify `mod_names`, `aa_mass` or `aa_comp`.".format(label))
            name = label.capitalize()
        return [proforma.GenericModification(prefix + name)]

    i, j = 0, len(sequence)
    nterm = cterm = None
    pro_sequence = []
    if isinstance(sequence[0], str):  # regular parsed sequence
        if is_term_mod(sequence[0]) and sequence[0] != std_nterm:
            nterm = get_tag(sequence[0])
            i = 1
        if is_term_mod(sequence[-1]) and sequence[-1] != std_cterm:
            cterm = get_tag(sequence[-1])
            j -= 1
        for label in sequence[i:j]:
            if len(label) == 1:
                pro_sequence.append((label, None))
            else:
                mod, aa = _split_label(label)
                pro_sequence.append((aa, get_tag(mod)))
    else:  # split sequence
        if is_term_mod(sequence[0][0]) and sequence[0][0] != std_nterm:
            nterm = get_tag(sequence[0][0])
        if is_term_mod(sequence[-1][-1]) and sequence[-1][-1] != std_cterm:
            cterm = get_tag(sequence[-1][-1])
        if len(sequence) == 1:
            pro_sequence = [(sequence[0][-2] if cterm else sequence[0][-1], get_tag(sequence[0][1]) if len(sequence[0]) == 4 else None)]
        else:
            pro_sequence.append((sequence[0][-1], get_tag(sequence[0][-2]) if len(sequence[0]) == 3 else None))
            for group in sequence[1:-1]:
                pro_sequence.append((group[-1], get_tag(group[0]) if len(group) == 2 else None))
            if len(sequence[-1]) == 1 or (len(sequence[-1]) == 2 and cterm):
                pro_sequence.append((sequence[-1][0], None))
            else:
                pro_sequence.append((sequence[-1][1], get_tag(sequence[-1][0])))

    return proforma.to_proforma(pro_sequence, n_term=nterm, c_term=cterm)


def amino_acid_composition(sequence, show_unmodified_termini=False, term_aa=False, allow_unknown_modifications=False, **kwargs):
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
def cleave(*args, **kwargs):
    """Cleaves a polypeptide sequence using a given rule.

    .. seealso::
        :func:`icleave` and :func:`xcleave`, which produce both peptides and their indices.

    Parameters
    ----------
    sequence : str
        The sequence of a polypeptide.

        .. note::
            The sequence is expected to be in one-letter uppercase notation.
            Otherwise, some of the cleavage rules in :py:data:`expasy_rules`
            will not work as expected.

    rule : str or compiled regex
        A key present in :py:data:`expasy_rules`, :py:data:`psims_rules` (or an MS ontology accession) or a
        `regular expression <https://docs.python.org/library/re.html#regular-expression-syntax>`_
        describing the site of cleavage. It is recommended
        to design the regex so that it matches only the residue whose C-terminal
        bond is to be cleaved. All additional requirements should be specified
        using `lookaround assertions
        <http://www.regular-expressions.info/lookaround.html>`_.
        :py:data:`expasy_rules` contains cleavage rules for popular cleavage agents.

        .. seealso:: The `regex` argument.

    missed_cleavages : int, optional
        Maximum number of allowed missed cleavages. Defaults to 0.
    min_length : int or None, optional
        Minimum peptide length. Defaults to :py:const:`None`.

        .. note ::
            This checks for string length, which is only correct for one-letter
            notation and not for full *modX*. Use :py:func:`length` manually if
            you know what you are doing and apply :py:func:`cleave` to *modX*
            sequences.

    max_length : int or None, optional
        Maximum peptide length. Defaults to :py:const:`None`. See note above.

    semi : bool, optional
        Include products of semi-specific cleavage. Default is :py:const:`False`.
        This effectively cuts every peptide at every position and adds results to the output.

    exception : str or compiled RE or None, optional
        Exceptions to the cleavage rule. If specified, should be a key present in :py:const:`expasy_rules`
        or regular expression. Cleavage sites matching `rule` will be checked against `exception` and omitted
        if they match.

    regex : bool, optional
        If :py:const:`True`, the cleavage rule is always interpreted as a regex. Otherwise, a matching value
        is looked up in :py:data:`expasy_rules` and :py:data:`psims_rules`.

    Returns
    -------
    out : set
        A set of unique (!) peptides.

    Examples
    --------
    >>> cleave('AKAKBK', expasy_rules['trypsin'], 0) == {'AK', 'BK'}
    True
    >>> cleave('AKAKBK', 'trypsin', 0) == {'AK', 'BK'}
    True
    >>> cleave('AKAKBK', 'MS:1001251', 0) == {'AK', 'BK'}
    True
    >>> cleave('GKGKYKCK', 'Trypsin/P', 2) == \
    {'CK', 'GKYK', 'YKCK', 'GKGK', 'GKYKCK', 'GK', 'GKGKYK', 'YK'}
    True

    """
    return set(p for i, p in icleave(*args, **kwargs))


def icleave(sequence, rule, missed_cleavages=0, min_length=None, max_length=None, semi=False, exception=None, regex=False):
    """Like :py:func:`cleave`, but the result is an iterator and includes peptide indices.
    Refer to :py:func:`cleave` for explanation of parameters.

    Returns
    -------
    out : iterator
        An iterator over (index, sequence) pairs.

    """
    if not regex:
        if rule in expasy_rules:
            rule = expasy_rules[rule]
        elif rule in psims_rules:
            rule = psims_rules[rule]
        elif rule in _psims_index:
            rule = _psims_index[rule]
        elif re.search(r'[a-z]', rule):
            warnings.warn('Interpreting the rule as a regular expression: {}. Did you mistype the rule? '
                'Specify `regex=True` to silence this warning.'.format(rule))
    exception = expasy_rules.get(exception, exception)
    ml = missed_cleavages + 2
    trange = range(ml)
    cleavage_sites = deque([0], maxlen=ml)
    if min_length is None:
        min_length = 1
    if max_length is None:
        max_length = len(sequence)
    cl = 1
    if exception is not None:
        exceptions = {x.end() for x in re.finditer(exception, sequence)}
    for end in it.chain([x.end() for x in re.finditer(rule, sequence)], [None]):
        if exception is not None and end in exceptions:
            continue
        cleavage_sites.append(end)
        if cl < ml:
            cl += 1
        for j in trange[:cl - 1]:
            seq = sequence[cleavage_sites[j]:cleavage_sites[-1]]
            lenseq = len(seq)
            if end is not None:
                start = end - lenseq
            else:
                start = len(sequence) - lenseq
            if seq and min_length <= lenseq <= max_length:
                yield (start, seq)
                if semi:
                    for k in range(min_length, min(lenseq, max_length)):
                        yield (start, seq[:k])
                    for k in range(max(1, lenseq - max_length), lenseq - min_length + 1):
                        yield (start + k, seq[k:])


def xcleave(*args, **kwargs):
    """Like :py:func:`icleave`, but returns a list.

    Returns
    -------
    out : list
        A list of (index, sequence) pairs.

    Examples
    --------
    >>> xcleave('AKAKBK', 'trypsin', 1)
    [(0, 'AK'), (0, 'AKAK'), (2, 'AK'), (2, 'AKBK'), (4, 'BK')]
    """
    return list(icleave(*args, **kwargs))


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
    exception : str or compiled RE or None, optional
        Exceptions to the cleavage rule. If specified, should be a regular expression.
        Cleavage sites matching `rule` will be checked against `exception` and omitted
        if they match.

    Returns
    -------
    out : int
        Number of cleavage sites.
    """
    return sum(1 for _ in icleave(sequence, rule, **kwargs)) - 1


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
    'pepsin ph1.3':  r'((?<=[^HKR][^P])[^R](?=[FL][^P]))|'
                     r'((?<=[^HKR][^P])[FL](?=\w[^P]))',
    'pepsin ph2.0':  r'((?<=[^HKR][^P])[^R](?=[FLWY][^P]))|'
                     r'((?<=[^HKR][^P])[FLWY](?=\w[^P]))',
    'proline endopeptidase': r'(?<=[HKR])P(?=[^P])',
    'proteinase k':  r'[AEFILTVWY]',
    'staphylococcal peptidase i': r'(?<=[^E])E',
    'thermolysin':   r'[^DE](?=[AFILMV][^P])',
    'thrombin':      r'((?<=G)R(?=G))|'
                     r'((?<=[AFGILTVM][AFGILTVWA]P)R(?=[^DE][^DE]))',
    'trypsin':       r'([KR](?=[^P]))|((?<=W)K(?=P))|((?<=M)R(?=P))',
    'trypsin_exception': r'((?<=[CD])K(?=D))|((?<=C)K(?=[HY]))|((?<=C)R(?=K))|((?<=R)R(?=[HR]))',
}
"""
This dict contains regular expressions for cleavage rules of the most
popular proteolytic enzymes. The rules were taken from the
`PeptideCutter tool
<http://ca.expasy.org/tools/peptidecutter/peptidecutter_enzymes.html>`_
at Expasy.

.. note::
    'trypsin_exception' can be used as `exception` argument when calling
    :py:func:`cleave` with 'trypsin' `rule`::

        >>> parser.cleave('PEPTIDKDE', parser.expasy_rules['trypsin'])
        {'DE', 'PEPTIDK'}
        >>> parser.cleave('PEPTIDKDE', parser.expasy_rules['trypsin'], \
exception=parser.expasy_rules['trypsin_exception'])
        {'PEPTIDKDE'}
"""


psims_rules = {
    cvstr('2-iodobenzoate', 'MS:1001918'): r'(?<=W)',
    cvstr('Arg-C', 'MS:1001303'): r'(?<=R)(?!P)',
    cvstr('Asp-N', 'MS:1001304'): r'(?=[BD])',
    cvstr('Asp-N ambic', 'MS:1001305'): r'(?=[DE])',
    cvstr('CNBr', 'MS:1001307'): r'(?<=M)',
    cvstr('Chymotrypsin', 'MS:1001306'): r'(?<=[FYWL])(?!P)',
    cvstr('Formic acid', 'MS:1001308'): r'((?<=D))|((?=D))',
    cvstr('Lys-C', 'MS:1001309'): r'(?<=K)(?!P)',
    cvstr('Lys-C/P', 'MS:1001310'): r'(?<=K)',
    cvstr('PepsinA', 'MS:1001311'): r'(?<=[FL])',
    cvstr('TrypChymo', 'MS:1001312'): r'(?<=[FYWLKR])(?!P)',
    cvstr('Trypsin', 'MS:1001251'): r'(?<=[KR])(?!P)',
    cvstr('Trypsin/P', 'MS:1001313'): r'(?<=[KR])',
    cvstr('V8-DE', 'MS:1001314'): r'(?<=[BDEZ])(?!P)',
    cvstr('V8-E', 'MS:1001315'): r'(?<=[EZ])(?!P)',
    cvstr('glutamyl endopeptidase', 'MS:1001917'): r'(?<=[^E]E)',
    cvstr('leukocyte elastase', 'MS:1001915'): r'(?<=[ALIV])(?!P)',
    cvstr('proline endopeptidase', 'MS:1001916'): r'(?<=[HKR]P)(?!P)',
}
"""
This dict contains regular expressions for cleavage rules of the most
popular proteolytic enzymes. The rules were taken from the PSI `MS ontology
<http://purl.obolibrary.org/obo/MS_1001045>`_.

You can use names or accessions to access the rules.
Use :py:func:`pyteomics.auxiliary.cvquery` for accession access::

    >>> from pyteomics.auxiliary import cvquery
    >>> from pyteomics.parser import psims_rules
    >>> cvquery(psims_rules, 'MS:1001918')
    '(?<=W)'

"""

_psims_index = cvquery(psims_rules)

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
    def main(group):  # index of the residue (capital letter) in `group`
        if group[-1][0] == '-':
            i = -2
        else:
            i = -1
        return len(group) + i, group[i]

    def apply_mod(label, mod):
        # `label` is assumed to be a tuple (see split option of `parse`)
        # unmodified termini are assumed shown
        # if the modification is not applicable, `None` is returned
        group = list(label)
        m = main(group)[0]
        c = True  # whether the change is applied in the end
        if m == 0 and not is_term_mod(mod):
            group.insert(0, mod)
        elif mod[0] == '-' and (group[-1] == std_cterm or (group[-1][0] == '-' and override)):
            group[-1] = mod
        elif mod[-1] == '-' and (group[0] == std_nterm or (group[0][-1] == '-' and override)):
            group[0] = mod
        elif not is_term_mod(mod):
            if m and group[m - 1][-1] != '-':
                if override:
                    group[m - 1] = mod
                else:
                    c = False
            else:
                group.insert(m, mod)
        else:
            c = False
        if c:
            return tuple(group)

    variable_mods = kwargs.get('variable_mods', {})
    varmods_term, varmods_non_term = [], []
    for m, r in sorted(variable_mods.items()):
        if is_term_mod(m):
            varmods_term.append((m, r))
        else:
            varmods_non_term.append((m, r))
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
    for cmod, res in fixed_mods.items():
        for i, group in enumerate(parsed):
            if res is True or main(group)[1] in res:
                parsed[i] = apply_mod(group, cmod) or parsed[i]

    # Create a list of possible states for each group
    # Start with N-terminal mods and regular mods on the N-terminal residue
    states = [[parsed[0]]]
    m0 = main(parsed[0])[1]
    for m, r in varmods_non_term:
        if r is True or m0 in r or 'nterm' + m0 in r or len(parsed) == 1 and 'cterm' + m0 in r:
            applied = apply_mod(parsed[0], m)
            if applied is not None:
                states[0].append(applied)
    more_states = []
    for m, r in varmods_term:
        if r is True or m0 in r:
            if m[-1] == '-' or len(parsed) == 1:
                for group in states[0]:
                    applied = apply_mod(group, m)
                    if applied is not None:
                        more_states.append(applied)
    states[0].extend(more_states)

    # Continue with regular mods
    for group in parsed[1:-1]:
        gstates = [group]
        for m, r in varmods_non_term:
            if r is True or group[-1] in r:
                applied = apply_mod(group, m)
                if applied is not None:
                    gstates.append(applied)
        states.append(gstates)

    # Finally add C-terminal mods and regular mods on the C-terminal residue
    if len(parsed) > 1:
        states.append([parsed[-1]])
        m1 = main(parsed[-1])[1]
        for m, r in varmods_non_term:
            if r is True or m1 in r or 'cterm' + m1 in r or len(parsed) == 1 and 'nterm' + m1 in r:
                applied = apply_mod(parsed[-1], m)
                if applied is not None:
                    states[-1].append(applied)
        more_states = []
        for m, r in varmods_term:
            if r is True or m1 in r:
                if m[0] == '-' or len(parsed) == 1:
                    for group in states[-1]:
                        applied = apply_mod(group, m)
                        if applied is not None:
                            more_states.append(applied)
        states[-1].extend(more_states)

    sites = [s for s in enumerate(states) if len(s[1]) > 1]
    if max_mods is None or max_mods > len(sites):
        possible_states = it.product(*states)
    else:
        def state_lists():
            for m in range(max_mods + 1):
                for comb in it.combinations(sites, m):
                    skel = [[s[0]] for s in states]
                    for i, e in comb:
                        skel[i] = e[1:]
                    yield skel
        possible_states = it.chain.from_iterable(it.product(*skel) for skel in state_lists())

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
            mask[i:i + len(peptide)] = 1
    return mask.sum(dtype=float) / mask.size


if __name__ == "__main__":
    import doctest
    doctest.testmod()
