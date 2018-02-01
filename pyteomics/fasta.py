"""
fasta - manipulations with FASTA databases
==========================================

FASTA is a simple file format for protein sequence databases. Please refer to
`the NCBI website <http://www.ncbi.nlm.nih.gov/blast/fasta.shtml>`_
for the most detailed information on the format.

Data manipulation
-----------------

  :py:func:`read` - iterate through entries in a FASTA database.

  :py:func:`chain` - read multiple files at once.

  :py:func:`chain.from_iterable` - read multiple files at once, using an
  iterable of files.

  :py:func:`write` - write entries to a FASTA database.

  :py:func:`parse` - parse a FASTA header.

Decoy sequence generation
-------------------------

:py:func:`decoy_sequence` - generate a decoy sequence from a given sequence, using
one of the other functions listed in this section or any other callable.

:py:func:`reverse` - generate a reversed decoy sequence.

:py:func:`shuffle` - generate a shuffled decoy sequence.

:py:func:`fused_decoy` - generate a "fused" decoy sequence.


Decoy database generation
-------------------------
  
  :py:func:`write_decoy_db` - generate a decoy database and write it to a file.

  :py:func:`decoy_db` - generate entries for a decoy database from a given FASTA
  database.

  :py:func:`decoy_chain` - a version of :py:func:`decoy_db` for multiple files.

  :py:func:`decoy_chain.from_iterable` - like :py:func:`decoy_chain`, but with
  an iterable of files.

Auxiliary
---------

  :py:data:`std_parsers` - a dictionary with parsers for known FASTA header
  formats.

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

import itertools
import random
from collections import namedtuple
import re
from . import auxiliary as aux

Protein = namedtuple('Protein', ('description', 'sequence'))

@aux._file_reader()
def read(source=None, ignore_comments=False, parser=None):
    """Read a FASTA file and return entries iteratively.

    Parameters
    ----------
    source : str or file or None, optional
        A file object (or file name) with a FASTA database. Default is
        :py:const:`None`, which means read standard input.
    ignore_comments : bool, optional
        If True then ignore the second and subsequent lines of description.
        Default is :py:const:`False`.
    parser : function or None, optional
        Defines whether the fasta descriptions should be parsed. If it is a
        function, that function will be given the description string, and
        the returned value will be yielded together with the sequence.
        The :py:data:`std_parsers` dict has parsers for several formats.
        Hint: specify :py:func:`parse` as the parser to apply automatic
        format recognition.
        Default is :py:const:`None`, which means return the header "as is".

    Returns
    -------
    out : iterator of tuples
        A named 2-tuple with FASTA header (str or dict) and sequence (str).
        Attributes 'description' and 'sequence' are also provided.
    """
    f = parser or (lambda x: x)
    accumulated_strings = []

    # Iterate through '>' after the file is over to retrieve the last entry.
    for string in itertools.chain(source, '>'):
        stripped_string = string.strip()

        # Skip empty lines.
        if not stripped_string:
            continue

        is_comment = (stripped_string[0] in '>;')
        if is_comment:
            # If it is a continuing comment
            if len(accumulated_strings) == 1:
                if not ignore_comments:
                    accumulated_strings[0] += (' '+stripped_string[1:])
                else:
                    continue

            elif accumulated_strings:
                description = accumulated_strings[0]
                sequence = ''.join(accumulated_strings[1:])

                # Drop the translation stop sign.
                if sequence.endswith('*'):
                    sequence = sequence[:-1]
                yield Protein(f(description), sequence)
                accumulated_strings = [stripped_string[1:], ]
            else:
                # accumulated_strings is empty; we're probably reading
                # the very first line of the file
                accumulated_strings.append(stripped_string[1:])
        else:
            accumulated_strings.append(stripped_string)

@aux._file_writer()
def write(entries, output=None):
    """
    Create a FASTA file with `entries`.

    Parameters
    ----------
    entries : iterable of (str, str) tuples
        An iterable of 2-tuples in the form (description, sequence).
    output : file-like or str, optional
        A file open for writing or a path to write to. If the file exists,
        it will be opened for appending. Default is :py:const:`None`, which
        means write to standard output.
    file_mode : str, keyword only, optional
        If `output` is a file name, defines the mode the file will be opened in.
        Otherwise will be ignored. Default is 'a'.

    Returns
    -------
    output_file : file object
        The file where the FASTA is written.
    """
    for descr, seq in entries:
        output.write('>' + descr.replace('\n', '\n;') + '\n')
        output.write(''.join([('%s\n' % seq[i:i+70])
            for i in range(0, len(seq), 70)]) + '\n')

    return output.file

def reverse(sequence, keep_nterm=False, keep_cterm=False):
    """
    Create a decoy sequence by reversing the original one.

    Parameters
    ----------
    sequence : str
        The initial sequence string.
    keep_nterm : bool, optional
        If :py:const:`True`, then the N-terminal residue will be kept.
        Default is :py:const:`False`.
    keep_cterm : bool, optional
        If :py:const:`True`, then the C-terminal residue will be kept.
        Default is :py:const:`False`.

    Returns
    -------
    decoy_sequence : str
        The decoy sequence.
    """
    start = 1 if keep_nterm else 0
    end = len(sequence)-1 if keep_cterm else len(sequence)
    if start == end:
        return sequence
    return sequence[:start] + sequence[start:end][::-1] + sequence[end:]

def shuffle(sequence, keep_nterm=False, keep_cterm=False):
    """
    Create a decoy sequence by shuffling the original one.

    Parameters
    ----------
    sequence : str
        The initial sequence string.
    keep_nterm : bool, optional
        If :py:const:`True`, then the N-terminal residue will be kept.
        Default is :py:const:`False`.
    keep_cterm : bool, optional
        If :py:const:`True`, then the C-terminal residue will be kept.
        Default is :py:const:`False`.

    Returns
    -------
    decoy_sequence : str
        The decoy sequence.
    """
    start = 1 if keep_nterm else 0
    end = len(sequence)-1 if keep_cterm else len(sequence)
    if start == end:
        return sequence
    elif keep_cterm or keep_nterm:
        return sequence[:start] + shuffle(sequence[start:end]) + sequence[end:]

    modified_sequence = list(sequence)
    random.shuffle(modified_sequence)
    return ''.join(modified_sequence)

def fused_decoy(sequence, decoy_mode='reverse', sep='R', **kwargs):
    """
    Create a "fused" decoy sequence by concatenating a decoy sequence with the original one.
    The method and its use cases are described in:

    Ivanov, M. V., Levitsky, L. I., & Gorshkov, M. V. (2016).
    `Adaptation of Decoy Fusion Strategy for Existing Multi-Stage Search Workflows.
    <http://doi.org/10.1007/s13361-016-1436-7>`_
    Journal of The American Society for Mass Spectrometry, 27(9), 1579-1582.

    Parameters
    ----------
    sequence : str
        The initial sequence string.
    decoy_mode : str or callable, optional
        Type of decoy sequence to use. Should be one of the standard modes or any callable.
        Standard modes are:

        - 'reverse' for :py:func:`reverse`;
        - 'shuffle' for :py:func:`shuffle`;
        - 'fused' for :py:func:`fused_decoy` (if you love recursion).

        Default is 'reverse'.
    sep : str, optional
        Amino acid motif that separates the decoy sequence from the target one.
        This setting should reflect the enzyme specificity used in the search against the
        database being generated. Default is 'R', which is suitable for trypsin searches.
    **kwargs : given to the decoy generation function.

    Examples
    --------
    >>> fused_decoy('PEPT')
    'TPEPRPEPT'
    >>> fused_decoy('MPEPT', 'shuffle', 'K', keep_nterm=True)
    'MPPTEKMPEPT'
    """
    decoy = decoy_sequence(sequence, decoy_mode, **kwargs)
    return decoy + sep + sequence

_decoy_functions = {'reverse': reverse, 'shuffle': shuffle, 'fused': fused_decoy}

def decoy_sequence(sequence, mode='reverse', **kwargs):
    """
    Create a decoy sequence out of a given sequence string.

    Parameters
    ----------
    sequence : str
        The initial sequence string.
    mode : str or callable, optional
        Type of decoy sequence. Should be one of the standard modes or any callable.
        Standard modes are:

        - 'reverse' for :py:func:`reverse`;
        - 'shuffle' for :py:func:`shuffle`;
        - 'fused' for :py:func:`fused_decoy`.

        Default is 'reverse'.
    **kwargs : given to the decoy function.

    Returns
    -------
    decoy_sequence : str
        The decoy sequence.
    """
    fmode = mode
    if isinstance(mode, str):
        fmode = _decoy_functions.get(mode)
        if fmode is None:
            raise aux.PyteomicsError('Unsupported decoy mode: {}'.format(mode))
    return fmode(sequence, **kwargs)

@aux._file_reader()
def decoy_db(source=None, mode='reverse', prefix='DECOY_', decoy_only=False,
        ignore_comments=False, parser=None, **kwargs):
    """Iterate over sequences for a decoy database out of a given ``source``.

    Parameters
    ----------
    source : file-like object or str or None, optional
        A path to a FASTA database or a file object itself. Default is
        :py:const:`None`, which means read standard input.
    mode : str or callable, optional
        Algorithm of decoy sequence generation. 'reverse' by default.
        See :py:func:`decoy_sequence` for more information.
    prefix : str, optional
        A prefix to the protein descriptions of decoy entries. The default
        value is `'DECOY_'`.
    decoy_only : bool, optional
        If set to :py:const:`True`, only the decoy entries will be written to
        `output`. If :py:const:`False`, the entries from `source` will be
        written first.
        :py:const:`False` by default.
    ignore_comments : bool, optional
        If True then ignore the second and subsequent lines of description.
        Default is :py:const:`False`.
    parser : function or None, optional
        Defines whether the fasta descriptions should be parsed. If it is a
        function, that function will be given the description string, and
        the returned value will be yielded together with the sequence.
        The :py:data:`std_parsers` dict has parsers for several formats.
        Hint: specify :py:func:`parse` as the parser to apply automatic
        format guessing.
        Default is :py:const:`None`, which means return the header "as is".
    **kwargs : given to :py:func:`decoy_sequence`.

    Returns
    -------
    out : iterator
        An iterator over entries of the new database.
    """

    # store the initial position
    pos = source.tell()
    if not decoy_only:
        with read(source, ignore_comments, parser) as f:
            for x in f:
                yield x
        # return to the initial position in the source file to read again
        source.seek(pos)

    parser = parser or (lambda x: x)
    with read(source, ignore_comments) as f:
        for descr, seq in f:
            yield Protein(parser(prefix + descr), decoy_sequence(seq, mode, **kwargs))


@aux._file_writer()
def write_decoy_db(source=None, output=None, mode='reverse', prefix='DECOY_',
        decoy_only=False, **kwargs):
    """Generate a decoy database out of a given ``source`` and write to file.

    If `output` is a path, the file will be open for appending, so no information
    will be lost if the file exists. Although, the user should be careful when
    providing open file streams as `source` and `output`. The reading and writing
    will start from the current position in the files, which is where the last I/O
    operation finished. One can use the :py:func:`file.seek` method to change it.

    Parameters
    ----------
    source : file-like object or str or None, optional
        A path to a FASTA database or a file object itself. Default is
        :py:const:`None`, which means read standard input.
    output : file-like object or str, optional
        A path to the output database or a file open for writing.
        Defaults to :py:const:`None`, the results go to the standard output.
    mode : str or callable, optional
        Algorithm of decoy sequence generation. 'reverse' by default.
        See :py:func:`decoy_sequence` for more details.
    prefix : str, optional
        A prefix to the protein descriptions of decoy entries. The default
        value is `'DECOY_'`
    decoy_only : bool, optional
        If set to :py:const:`True`, only the decoy entries will be written to
        `output`. If :py:const:`False`, the entries from `source` will be
        written as well.
        :py:const:`False` by default.
    file_mode : str, keyword only, optional
        If `output` is a file name, defines the mode the file will be opened in.
        Otherwise will be ignored. Default is 'a'.
    **kwargs : given to :py:func:`decoy_sequence`.

    Returns
    -------
    output : file
        A (closed) file object for the created file.
    """
    with decoy_db(source, mode, prefix, decoy_only, **kwargs) as entries:
        write(entries, output)
        return output.file

# auxiliary functions for parsing of FASTA headers
def _split_pairs(s):
    return dict(map(lambda x: x.strip(), x.split('='))
            for x in re.split(' (?=\w+=)', s.strip()))

def _intify(d, keys):
    for k in keys:
        if k in d:
            d[k] = int(d[k])

# definitions for custom parsers
def _parse_uniprotkb(header):
    db, ID, entry, name, pairs, _ = re.match(
           r'^(\w+)\|([-\w]+)\|(\w+)\s+([^=]*\S)((\s+\w+=[^=]+(?!\w*=))+)\s*$',
           header).groups()
    gid, taxon = entry.split('_')
    info = {'db': db, 'id': ID, 'entry': entry,
            'name': name, 'gene_id': gid, 'taxon': taxon}
    info.update(_split_pairs(pairs))
    _intify(info, ('PE', 'SV'))
    return info

def _parse_uniref(header):
    assert 'Tax' in header
    ID, cluster, pairs, _ = re.match(
            r'^(\S+)\s+([^=]*\S)((\s+\w+=[^=]+(?!\w*=))+)\s*$',
            header).groups()
    info = {'id': ID, 'cluster': cluster}
    info.update(_split_pairs(pairs))
    gid, taxon = info['RepID'].split('_')
    type_, acc = ID.split('_')
    info.update({'taxon': taxon, 'gene_id': gid,
            'type': type_, 'accession': acc})
    _intify(info, ('n',))
    return info

def _parse_uniparc(header):
    ID, status = re.match(r'(\S+)\s+status=(\w+)\s*$', header).groups()
    return {'id': ID, 'status': status}

def _parse_unimes(header):
    assert 'OS=' in header and 'SV=' in header and 'PE=' not in header
    ID, name, pairs, _ = re.match(
            r'^(\S+)\s+([^=]*\S)((\s+\w+=[^=]+(?!\w*=))+)\s*$',
            header).groups()
    info = {'id': ID, 'name': name}
    info.update(_split_pairs(pairs))
    _intify(info, ('SV',))
    return info

def _parse_spd(header):
    assert '=' not in header
    ID, gene, d = map(lambda s: s.strip(), header.split('|'))
    gid, taxon = gene.split('_')
    return {'id': ID, 'gene': gene, 'description': d,
            'taxon': taxon, 'gene_id': gid}

std_parsers = {'uniprotkb': _parse_uniprotkb, 'uniref': _parse_uniref,
        'uniparc': _parse_uniparc, 'unimes': _parse_unimes, 'spd': _parse_spd}
"""A dictionary with parsers for known FASTA header formats. For now, supported
formats are those described at
`UniProt help page <http://www.uniprot.org/help/fasta-headers>`_."""

def parse(header, flavour='auto', parsers=None):
    """Parse the FASTA header and return a nice dictionary.

    Parameters
    ----------

    header : str
        FASTA header to parse
    flavour : str, optional
        Short name of the header format (case-insensitive). Valid values are
        :py:const:`'auto'` and keys of the `parsers` dict. Default is
        :py:const:`'auto'`, which means try all formats in turn and return the
        first result that can be obtained without an exception.
    parsers : dict, optional
        A dict where keys are format names (lowercased) and values are functions
        that take a header string and return the parsed header. Default is
        :py:const:`None`, which means use the default dictionary
        :py:data:`std_parsers`.

    Returns
    -------

    out : dict
        A dictionary with the info from the header. The format depends on the
        flavour."""

    # accept strings with and without leading '>'
    if header.startswith('>'):
        header = header[1:]

    # choose the format
    known = parsers or std_parsers
    if flavour.lower() == 'auto':
        for fl, parser in known.items():
            try:
                return parser(header)
            except:
                pass
        raise aux.PyteomicsError('Unknown FASTA header format: ' + header)
    elif flavour.lower() in known:
        try:
            return known[flavour.lower()](header)
        except Exception as e:
            raise aux.PyteomicsError('Could not parse header as "{}". '
                    'The error message was: {}: {}. Header: "{}"'.format(
                        flavour, type(e).__name__, e.args[0], header))

chain = aux._make_chain(read, 'read')
decoy_chain = aux._make_chain(decoy_db, 'decoy_db')
