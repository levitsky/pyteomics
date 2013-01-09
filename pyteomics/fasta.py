"""
fasta - manipulations with FASTA databases 
==========================================

FASTA is a simple file format for protein sequence databases. Please refer to
`the NCBI website <http://www.ncbi.nlm.nih.gov/blast/fasta.shtml>`_
for the most detailed information on the format.

Data manipulation
-----------------

  :py:func:`read` - iterate through entries in a FASTA database

  :py:func:`write` - write entries to a FASTA database

  :py:func:`parse` - parse a FASTA header

Decoy database generation
-------------------------

  :py:func:`decoy_sequence` - generate a decoy sequence from a given sequence

  :py:func:`decoy_db` - generate a decoy database from a given FASTA database

Auxiliary
----------

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
import re
from .auxiliary import PyteomicsError, _file_obj, _file_reader

@_file_reader()
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
        format guessing.
        Default is :py:const:`None`, which means return the header "as is".

    Yields
    ------
    
    description, sequence
        A 2-tuple with FASTA header (str) and sequence (str).
    """
    f = parser or (lambda x: x)
    accumulated_strings = []

    # Iterate through '>' after the file is over to retrieve the last entry.
    for string in itertools.chain(source, '>'):
        stripped_string = string.strip()

        # Skip empty lines.
        if not stripped_string:
            continue
        
        is_comment = (stripped_string.startswith('>')
                      or stripped_string.startswith(';'))
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
                yield f(description), sequence
                accumulated_strings = [stripped_string[1:], ]
            else:
                # accumulated_strings is empty; we're probably reading
                # the very first line of the file
                accumulated_strings.append(stripped_string[1:])
        else:
            accumulated_strings.append(stripped_string)

def write(entries, output=None):
    """
    Create a FASTA file with ``entries``.

    Parameters
    ----------
    entries : iterable of (str, str) tuples
        An iterable of 2-tuples in the form (description, sequence).
    output : file-like or str, optional
        A file open for writing or a path to write to. If the file exists,
        it will be opened for appending. Default is :py:const:`None`, which
        means write to standard output.

    Returns
    -------
    output_file : file object
        The file where the FASTA is written.
    """

    with _file_obj(output, 'a') as foutput:
        for descr, seq in entries:
            # write the description
            foutput.write('>' + descr.replace('\n', '\n;') + '\n')
            # write the sequence; it should be interrupted with \n every 70 characters
            foutput.write(''.join([('%s\n' % seq[i:i+70])
                for i in range(0, len(seq), 70)]) + '\n')

        return foutput.file

def decoy_sequence(sequence, mode):
    """
    Create a decoy sequence out of a given sequence string.

    Parameters
    ----------
    sequence : str
        The initial sequence string.
    mode : {'reverse', 'shuffle'}
        Type of decoy sequence. 

    Returns
    -------
    modified_sequence : str
        The modified sequence.
    """
    if mode == 'reverse':
        return sequence[::-1]
    if mode == 'shuffle':
        modified_sequence = list(sequence)
        random.shuffle(modified_sequence)
        return ''.join(modified_sequence)
    raise PyteomicsError(
            """`fasta.decoy_sequence`: `mode` must be 'reverse' or
            'shuffle', not {}""".format(mode))

def decoy_db(source=None, output=None, mode='reverse', prefix='DECOY_',
             decoy_only=False):
    """Generate a decoy database out of a given ``source``.
    
    If 'output' is a path, the file will be open for appending, so no information
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
    mode : {'reverse', 'shuffle'}, optional
        Algorithm of decoy sequence generation. 'reverse' by default.
    prefix : str, optional
        A prefix to the protein descriptions of decoy entries. The default 
        value is "DECOY\_"
    decoy_only : bool, optional
        If set to True, only the decoy entries will be written to 'output'.
        If False, the entries from 'source' will be written as well.
        False by default.
   
    Returns
    -------
    output : file
        A file object for the created file.
    """
    with _file_obj(source, 'r') as fsource, _file_obj(output, 'a') as foutput:

        # store the initial position
        pos = fsource.tell()
        if not decoy_only:
            write(read(fsource), foutput)

        # return to the initial position the source file to read again
        fsource.seek(pos)

        decoy_entries = ((prefix + descr, decoy_sequence(seq, mode))
            for descr, seq in read(fsource))

        write(decoy_entries, foutput)
        return foutput.file

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
    info = {'db': db, 'id': ID, 'entry': entry, 'name': name}
    info.update(_split_pairs(pairs))
    _intify(info, ('PE', 'SV'))
    return info

def _parse_uniref(header):
    assert 'Tax' in header and 'RepID=' in header
    ID, cluster, pairs, _ = re.match(
            r'^(\S+)\s+([^=]*\S)((\s+\w+=[^=]+(?!\w*=))+)\s*$',
            header).groups()
    info = {'id': ID, 'cluster': cluster}
    info.update(_split_pairs(pairs))
    _intify(info, ('n',))
    return info

def _parse_uniparc(header):
    ID, status = re.match(r'(\S+)\s+status=(\w+)\s*$', header).groups()
    return {'id': ID, 'status': status}

def _parse_unimes(header):
    assert 'OS=' in header and 'SV=' in header
    ID, name, pairs, _ = re.match(
            r'^(\S+)\s+([^=]*\S)((\s+\w+=[^=]+(?!\w*=))+)\s*$',
            header).groups()
    info = {'id': ID, 'name': name}
    info.update(_split_pairs(pairs))
    _intify(info, ('SV',))
    return info

std_parsers = {'uniprotkb': _parse_uniprotkb, 'uniref': _parse_uniref,
        'uniparc': _parse_uniparc, 'unimes': _parse_unimes}
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
        raise PyteomicsError('Unknown FASTA header format.')
    elif flavour.lower() in known:
        try:
            return known[flavour.lower()](header)
        except Exception as e:
            raise PyteomicsError('Could not parse as {}. '
                    'The error message was: {}'.format(
                        flavour, e.message))
