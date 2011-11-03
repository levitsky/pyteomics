"""
fasta
======

FASTA is a simple file format for protein sequence databases. More information
can be found at http://www.ncbi.nlm.nih.gov/blast/fasta.shtml.
"""

# Licensed under the MIT license:
# http://www.opensource.org/licenses/mit-license.php

# The FASTA format description can be found at
# http://www.ncbi.nlm.nih.gov/blast/fasta.shtml 

import itertools
import sys
from auxiliary import PyteomicsError

def read_fasta(fasta_file, ignore_comments=True):
    """Read a FASTA file and return entries iteratively.

    Parameters
    ----------
    fasta_file : str or file
        A file object (or file name) with a FASTA database.
    ignore_comments : bool, optional
        If True then ignore the second and subsequent lines of description.

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

def decoy_sequence(sequence, mode):
    """
    Modify a sequence string (used for decoy database generation).

    Parameters
    ----------
    sequence : str
        The initial sequence string.
    mode : str
        Type of modification. Supported values are: 'reverse'.

    Returns
    -------
    modified_sequence : str
        The modified sequence.
    """
    if mode == 'reverse':
        modified_sequence = sequence[::-1]
    return modified_sequence

def write_fasta(entries, output=None, close=True):
    """
    Create a FASTA file with 'entries'.

    Parameters
    ----------
    entries : list of (str, str) tuples
        A list containing 2-tuples in the form (description, sequence).
    output : file or str, optional
        A file open for writing or a path to write to. If the file exists,
        it will be open for appending.
    close : bool, optional
        If True, the file will be closed after writing. Defaults to True.

    Returns
    -------
    output_file : file object
        The file where the FASTA is written.
    """
    if type(output) == file:
        output_file = output
    elif type(output) == str:
        output_file = open(output, 'a')
    elif output == None:
        output_file = sys.stdout
    else:
        raise PyteomicsError("""Wrong argument type:
        'output' must be file or str or None, not %s""" % type(source))
    
    for protein in entries:
        # write the description
        output_file.write('>' + protein[0] + '\n')
        # write the sequence; it should be interrupted with \n every 70 characters
        output_file.write(''.join([('%s\n' % protein[1][i:i+70])
            for i in range(0, len(protein[1]), 70)]))
    if close and output: output_file.close()
    return output_file

def decoy_db(source, output=sys.stdout, mode='reverse', prefix='DECOY_', decoy_only=False):
    """Read the FASTA entries from 'source'. The file written to 'output' will
    contain the entries from 'source' and the modified sequences.
    
    Currently, the only supported type of modification is simple reversion.
    The protein description will be prepended with 'prefix' for modified
    sequences.

    If 'output' is a path, a file will be open for appending, so no information
    will be lost if the file exists.

    Parameters
    ----------
    source : file-like object or str
        A file with (or the filename of) a FASTA database
    output : file-like object or str, optional
        A file open for writing or a path to write to.
        If not specified, the results will be printed to standard output.
    mode : str, optional
        Algorithm of decoy sequence generation. Only 'reverse' is supported
        for now.
    prefix : str, optional
        A prefix that will prepend the protein descriptions for decoy entries.
        The default value is 'DECOY_'
    decoy_only : bool, optional
        If set to True, only the decoy entries will be written to 'output'.
        If False, the entries from 'source' will be written as well.

    """
    if type(source) == file:
        source_file = source
    elif type(source) == str:
        source_file = open(source)
    else:
        raise PyteomicsError("""Wrong argument type:
        'source' must be file or str, not %s""" % type(source))

    if type(output) == file:
        output_file = output
    elif type(output) == str:
        output_file = open(output, 'a')
    else:
        raise PyteomicsError("""Wrong argument type:
        'output' must be file or str, not %s""" % type(source))

    if not decoy_only:
        write_fasta(read_fasta(source), output_file, False)

    # return to the beginning of the source file to read again
    source_file.seek(0)

    decoy_entries = [(prefix + protein[0],
        decoy_sequence(protein[1], mode))
        for protein in read_fasta(source_file)]

    write_fasta(decoy_entries, output_file)
