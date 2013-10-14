"""
tandem - X!Tandem output file reader
====================================

Summary
-------

`X!Tandem <http://thegpm.org/tandem/>`_  is an open-source proteomic search
engine with a very simple, sophisticated application programming interface
(API): it simply takes an XML file of instructions on its command line,
and outputs the results into an XML file, which has been specified in the input
XML file. The output format is described
`here (PDF) <http://www.thegpm.org/docs/X_series_output_form.pdf>`_.

This module provides a minimalistic way to extract information from X!Tandem
output files. The main idea is the same as in :py:mod:`pyteomics.pepxml`:
the top-level function :py:func:`read` allows iterating over entries in
`<group>` elements, i.e. identifications for a certain spectrum.

Data access
-----------

  :py:func:`read` - iterate through peptide-spectrum matches in an X!Tandem
  output file. Data from a single PSM are converted to a human-readable dict.

  :py:func:`iterfind` - iterate over elements in an X!Tandem file.

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

import numpy as np
from . import auxiliary as aux

def _get_info_smart(source, element, **kw):
    return _get_info(source, element, **kw)

@aux._file_reader('rb')
def read(source):
    """Parse ``source`` and iterate through peptide-spectrum matches.

    Parameters
    ----------
    source : str or file
        A path to a target X!Tandem ooutput file or the file object itself.

    Returns
    -------
    out : iterator
       An iterator over dicts with PSM properties.
    """

    return iterfind(source, 'group[@type="model"]')

def _schema_info(source):
    """Stores defaults for X!Tandem output. Keys are: 'floats', 'ints',
    'bools', 'lists', 'intlists', 'floatlists', 'charlists'."""

    return {'ints': set(), 'floats': set(), 'bools': set(), 'lists': set(),
            'intlists': set(), 'floatlists': set(), 'charlists': set()}

_getinfo_env = {'keys': set(), 'schema_info': _schema_info,
    'get_info_smart': _get_info_smart}
_get_info = aux._make_get_info(_getinfo_env)

_iterfind_env = {'get_info_smart': _get_info_smart}
iterfind = aux._make_iterfind(_iterfind_env)
