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

Decoy database generation
-------------------------

  :py:func:`decoy_sequence` - generate a decoy sequence from a given sequence

  :py:func:`decoy_db` - generate a decoy database from a given FASTA database

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

#   The FASTA format description can be found at
#   http://www.ncbi.nlm.nih.gov/blast/fasta.shtml 

import itertools
import sys
import random
from .auxiliary import PyteomicsError

def read(fasta_file, ignore_comments=False):
    """Read a FASTA file and return entries iteratively.

    Parameters
    ----------

    fasta_file : str or file
        A file object (or file name) with a FASTA database.
    ignore_comments : bool, optional
        If True then ignore the second and subsequent lines of description.
        Default is :py:const:`False`.

    Yield a tuple (description, sequence).
    """
    accumulated_strings = []
    if type(fasta_file) == str:
        source = open(fasta_file)
    else:
        source = fasta_file
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
                yield (description, sequence)
                accumulated_strings = [stripped_string[1:], ]
            else:
                # accumulated_strings is empty; we're probably reading
                # the very first line of the file
                accumulated_strings.append(stripped_string[1:])
        else:
            accumulated_strings.append(stripped_string)

def write(entries, output=None, close=True):
    """
    Create a FASTA file with ``entries``.

    Parameters
    ----------
    entries : iterable of (str, str) tuples
        An iterable of 2-tuples in the form (description, sequence).
    output : file-like or str, optional
        A file open for writing or a path to write to. If the file exists,
        it will be open for appending. Default is :py:const:`None`, which
        means write to standard output.
    close : bool, optional
        If True, the file will be closed after writing. Defaults to True.

    Returns
    -------
    output_file : file object
        The file where the FASTA is written.
    """
    if hasattr(output, 'write'):
        output_file = output
    elif isinstance(output, str):
        output_file = open(output, 'a')
    elif output is None:
        output_file = sys.stdout
    else:
        raise PyteomicsError("""Wrong argument type:
        `output` must be file-like or str or None, not %s (%s)""" % (
            type(output), output))
    
    for protein in entries:
        # write the description
        output_file.write('>' + protein[0].replace('\n', '\n;') + '\n')
        # write the sequence; it should be interrupted with \n every 70 characters
        output_file.write(''.join([('%s\n' % protein[1][i:i+70])
            for i in range(0, len(protein[1]), 70)]) + '\n')

    if close and output: output_file.close()
    return output_file

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
        modified_sequence = sequence[::-1]
    elif mode == 'shuffle':
        modified_sequence = list(sequence)
        random.shuffle(modified_sequence)
        modified_sequence = ''.join(modified_sequence)
    else:
        raise PyteomicsError(
                """`fasta.decoy_sequence`: `mode` must be 'reverse' or
                'shuffle', not '%s'""" % mode)
    return modified_sequence

def decoy_db(source, output=None, mode='reverse', prefix='DECOY_',
             decoy_only=False, close=True):
    """Generate a decoy database out of a given ``source``.
    
    If 'output' is a path, the file will be open for appending, so no information
    will be lost if the file exists. Although, the user should be careful when
    providing open file streams as 'source' and 'output'. The reading and writing
    will start from the current position in the files, which is where the last I/O
    operation finished. One can use the :py:func:`file.seek` method to change it.

    Parameters
    ----------
    source : file-like object or str
        A path to a FASTA database or a file object itself.
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
    close : bool, optional
        If :py:const:`True`, the target file will be closed in the end (default).
        Set to False, if you need to perform I/O operations with the file
        from your Python script after it's created.
    
    Returns
    -------
    output_file : file
        A file object with the created file.
    """
    if hasattr(source, 'seek'):
        source_file = source
    elif isinstance(source, str):
        source_file = open(source)
    else:
        raise PyteomicsError("""Wrong argument type:
        `source` must be file or str, not %s""" % type(source))

    if isinstance(output, str):
        output_file = open(output, 'a')
    elif output is None:
        output_file = sys.stdout
    else:
        output_file = output

    if not decoy_only:
        write(read(source_file, False), output_file, False)

    # return to the beginning of the source file to read again
    source_file.seek(0)

    decoy_entries = ((prefix + protein[0],
        decoy_sequence(protein[1], mode))
        for protein in read(source_file, False))

    write(decoy_entries, output_file, close=(close if output else False))
    return output_file
