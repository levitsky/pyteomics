# -*- coding: utf-8 -*-
"""
Calculate GRand AVerage of hYdropathicity (GRAVY) index for amino acid sequence

Indices obtained from Kyte J., Doolittle F. J. Mol. Biol. 157:105-132 (1982)
"""
import re
from collections import Counter

"""
Hydropathicity indices
"""
AA_hidropaticity = {
    "A": 1.800,
    "R": -4.500,
    "N": -3.500,
    "D": -3.500,
    "C": 2.500,
    "Q": -3.500,
    "E": -3.500,
    "G": -0.400,
    "H": -3.200,
    "I": 4.500,
    "L": 3.800,
    "K": -3.900,
    "M": 1.900,
    "F": 2.800,
    "P": -1.600,
    "S": -0.800,
    "T": -0.700,
    "W": -0.900,
    "Y": -1.300,
    "V": 4.200,
}

# Regex to select letters that are not in standard FASTA format
invalidFasta = re.compile("[^{}]".format("".join(AA_hidropaticity.keys())))

"""
Calculate GRAVY index

Parameters:
   sequence
      string
      Amino acid sequence

Returns:
   numeric
   GRAVY index
"""


def gravy(sequence):
    sequence = re.sub(invalidFasta, "", sequence)
    return sum([AA_hidropaticity[aa] for aa in sequence]) / len(sequence)


"""
Calculate amino acid content of a sequence

Parameters:
   sequence
      string
      Amino acid sequence

   frequencies
      bool
      return amino acid counts (`False`) or amino acid frequencies (`True`)
      Default: `True`

Returns:
   `dict`
   Dictionary with amino acids as keys and counts/frequencies as values      
"""


def AAContent(sequence, frequencies=True):
    AACont = Counter(sequence)

    if frequencies:
        for k, v in AACont.iteritems():
            AACont[k] = float(v) / len(sequence)

    return AACont
