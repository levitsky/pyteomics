# Licensed under the MIT license:
# http://www.opensource.org/licenses/mit-license.php

# The FASTA format description can be found at
# http://www.ncbi.nlm.nih.gov/blast/fasta.shtml 

import itertools

def read_fasta(fasta_file, ignore_comments = True):
    """Read a FASTA file and return entries iteratively.

    Keyword arguments:
    fasta_file -- a file object with a fasta database.
    ignore_comments -- if True then ignore the second and following lines
                       of description.

    Yield a tuple (description, sequence).
    """
    accumulated_strings = []

    # Iterate through '>' after the file is over to retrieve the last entry.
    for string in itertools.chain(fasta_file, '>'):
        stripped_string = string.strip()

        # Skip empty lines.
        if not stripped_string:
            continue
        
        is_comment = (stripped_string.startswith('>')
                      or stripped_string.startswith(';'))
        
        if is_comment:
            # If it is a continuing comment
            if len(accumulated_strings) == 1 and not ignore_comments:
                accumulated_strings[0] += stripped_string[1:]
            else:
                description = accumulated_strings[0]
                sequence = ''.join(accumulated_strings[1:])

                # Drop the translation stop sign.
                if sequence.endswith('*'):
                    sequence = sequence[:-1]
                yield (description, sequence)
                accumulated_strings = [stripped_string[1:], ]
        else:
            accumulated_strings.append(stripped_string)

