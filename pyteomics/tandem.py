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

  :py:func:`chain` - read multiple files at once.

  :py:func:`chain.from_iterable` - read multiple files at once, using an
  iterable of files.

  :py:func:`filter` - iterate through peptide-spectrum matches in a chain of
  X!Tandem output files, yielding only top PSMs and keeping false discovery rate
  (FDR) at the desired level. The FDR is estimated using the target-decoy
  approach (TDA).

  :py:func:`filter.chain` - chain a series of filters applied independently to
  several files.

  :py:func:`filter.chain.from_iterable` - chain a series of filters applied
  independently to an iterable of files.

  :py:func:`iterfind` - iterate over elements in an X!Tandem file.

Auxiliary
---------

  :py:func:`is_decoy` - determine if a PSM is from the decoy database.

  :py:func:`fdr` - estimate the FDR in a data set using TDA.

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

import itertools as it
import operator
import numpy as np
from . import auxiliary as aux

def _get_info_smart(source, element, **kw):
    info = _get_info(source, element, **kw)
    # handy simplifications below
    if isinstance(info.get('note'), dict
            ) and set(info['note']) == {'label', 'note'}:
        info['note'] = info['note']['note']
    if 'protein' in info and 'label' in info:
        del info['label']
    if 'group' in info:
        for g in info['group']:
            label = g.pop('label')
            type_ = g.pop('type')
            info.setdefault(type_, {})[label] = g
        del info['group']
    if 'trace' in info:
        for t in info['trace']:
            info[t.pop('type')] = t
        del info['trace']
    if isinstance(info.get('values'), dict):
        info['values'] = info['values']['values']
    if isinstance(info.get('attribute'), list):
        for a in info.pop('attribute'):
            info[a['type']] = float(a['attribute'])
    if 'support' in info:
        for d in info['support']['supporting data'].values():
            for l in ['Xdata', 'Ydata']:
                d[l]['values'] = d[l]['values'].astype(int)
        fims = info['support']['fragment ion mass spectrum']
        fims.update(fims.pop('tandem mass spectrum'))
        for d in it.chain(info['support']['supporting data'].values(),
                (info['support']['fragment ion mass spectrum'],)):
            for l in ['Xdata', 'Ydata']:
                del d[l]['label']
    if 'charge' in info:
        info['charge'] = int(info['charge'])
    return info

@aux._file_reader('rb')
def read(source, read_schema=True):
    """Parse ``source`` and iterate through peptide-spectrum matches.

    Parameters
    ----------
    source : str or file
        A path to a target X!Tandem output file or the file object itself.

    read_schema : bool, optional
        If :py:const:`True`, attempt to extract information from the XML schema
        mentioned in the pepXML header (default). Otherwise, use default
        parameters. Disable this to avoid waiting on long network connections or
        if you don't like to get the related warnings.

    Returns
    -------
    out : iterator
       An iterator over dicts with PSM properties.
    """

    for g in iterfind(source, 'group[type=model]',
            recursive=True, read_schema=True):
        del g['type']
        yield g

def _schema_info(_, **kw):
    """Stores defaults for X!Tandem output. Keys are: 'floats', 'ints',
    'bools', 'lists', 'intlists', 'floatlists', 'charlists'."""

    return {'ints': {
        ('group', 'z'), ('aa', 'at')} | {('domain', k) for k in [
            'missed_cleavages', 'start', 'end', 'y_ions', 'b_ions',
            'a_ions', 'x_ions', 'c_ions', 'z_ions']},

            'floats': {('group', k) for k in [
                'fI', 'sumI', 'maxI', 'mh', 'expect', 'rt']} | {
                   ('domain', k) for k in [
                       'expect', 'hyperscore', 'b_score', 'y_score',
                       'a_score', 'x_score', 'c_score', 'z_score',
                       'nextscore', 'delta', 'mh']} | {
                   ('protein', 'expect'), ('protein', 'sumI'),
                   ('aa', 'modified')},

            'bools': set(),
            'lists': {'group', 'trace', 'attribute', 'protein', 'aa'},
            'floatlists': {('values', 'values')},
            'intlists': set(), 'charlists': set()}

_getinfo_env = {'keys': {'domain'}, 'schema_info': _schema_info,
    'get_info_smart': _get_info_smart}
_get_info = aux._make_get_info(_getinfo_env)

_iterfind_env = {'get_info_smart': _get_info_smart}
iterfind = aux._make_iterfind(_iterfind_env)

chain = aux._make_chain(read, 'read')

def is_decoy(psm, prefix='DECOY_'):
    """Given a PSM dict, return :py:const:`True` if all protein names for
    the PSM start with ``prefix``, and :py:const:`False` otherwise.

    Parameters
    ----------
    psm : dict
        A dict, as yielded by :py:func:`read`.
    prefix : str, optional
        A prefix used to mark decoy proteins. Default is "DECOY_".

    Returns
    -------
    out : bool
    """
    return all(prot['label'].startswith(prefix) for prot in psm['protein'])

filter = aux._make_filter(chain, is_decoy, operator.itemgetter('expect'))
fdr = aux._make_fdr(is_decoy)
filter.chain = aux._make_chain(filter, 'filter')
