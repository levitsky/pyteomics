"""
mzid - mzIdentML file reader
============================

Summary
-------

`mzIdentML <http://www.psidev.info/mzidentml>`_  is one of the standards
developed by the Proteomics Informatics working group of the HUPO Proteomics
Standard Initiative.

This module provides a minimalistic way to extract information from mzIdentML
files. The main idea is the same as in :py:mod:`pyteomics.pepxml`: the top-level
function :py:func:`read` allows iterating over entries in
`<SpectrumIdentificationResult>` elements, i.e. groups of identifications
for a certain spectrum. Note that each entry can contain more than one PSM
(peptide-spectrum match). They are accessible with "SpectrumIdentificationItem"
key.

Data access
-----------

  :py:func:`read` - iterate through peptide-spectrum matches in an mzIdentML
  file. Data from a single PSM group are converted to a human-readable dict.

  :py:func:`chain` - read multiple files at once.

  :py:func:`chain.from_iterable` - read multiple files at once, using an
  iterable of files.

  :py:func:`filter` - read a chain of mzIdentML files and filter to a certain
  FDR using TDA.

  :py:func:`filter.chain` - chain a series of filters applied independently to
  several files.

  :py:func:`filter.chain.from_iterable` - chain a series of filters applied
  independently to an iterable of files.

  :py:func:`get_by_id` - get an element by its ID and extract the data from it.

  :py:func:`iterfind` - iterate over elements in an mzIdentML file.

Miscellaneous
-------------

  :py:func:`version_info` - get information about mzIdentML version and schema.

  :py:func:`is_decoy` - determine if a "SpectrumIdentificationResult" should be
  consiudered decoy.

  :py:func:`fdr` - estimate the false discovery rate of a set of identifications
  using the target-decoy approach.

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

from lxml import etree
import numpy
from . import auxiliary as aux

def _get_info_smart(source, element, **kw):
    """Extract the info in a smart way depending on the element type"""
    name = aux._local_name(element)
    kwargs = dict(kw)
    rec = kwargs.pop('recursive', None)
    if name == 'MzIdentML':
        return _get_info(source, element, rec if rec is not None else False,
                **kwargs)
    else:
        return _get_info(source, element, rec if rec is not None else True,
                **kwargs)

@aux._keepstate
def get_by_id(source, elem_id):
    """Parse ``source`` and return the element with `id` attribute equal to
    ``elem_id``. Returns :py:const:`None` if no such element is found.

    Parameters
    ----------
    source : str or file
        A path to a target mzIdentML file of the file object itself.

    elem_id : str
        The value of the `id` attribute to match.

    Returns
    -------
    out : :py:class:`lxml.etree.Element` or :py:const:`None`
    """

    found = False
    for event, elem in etree.iterparse(source, events=('start', 'end'),
            remove_comments=True):
        if event == 'start':
            if elem.attrib.get('id') == elem_id:
                found = True
        else:
            if elem.attrib.get('id') == elem_id:
                return _get_info_smart(source, elem)
            if not found:
                elem.clear()
    return None

_version_info_env = {'format': 'mzIdentML', 'element': 'MzIdentML'}
version_info = aux._make_version_info(_version_info_env)

_schema_env = {'format': 'MzIdentML', 'version_info': version_info,
        'default_version': '1.1.0', 'defaults': aux._mzid_schema_defaults}
_schema_info = aux._make_schema_info(_schema_env)

# 'keys' should contain keys whose value is a dict
_get_info_env = {'keys':  {'Fragmentation',}, 'schema_info': _schema_info,
        'get_info_smart': _get_info_smart, 'get_by_id': get_by_id}
_get_info = aux._make_get_info(_get_info_env)

_iterfind_env = {'get_info_smart': _get_info_smart}
iterfind = aux._make_iterfind(_iterfind_env)

@aux._file_reader('rb')
def read(source, **kwargs):
    """Parse ``source`` and iterate through peptide-spectrum matches.

    Parameters
    ----------
    source : str or file
        A path to a target mzIdentML file or the file object itself.

    recursive : bool, optional
        If :py:const:`False`, subelements will not be processed when
        extracting info from elements. Default is :py:const:`True`.

    retrieve_refs : bool, optional
        If :py:const:`True`, additional information from references will be
        automatically added to the results. The file processing time will
        increase. Default is :py:const:`False`.

    read_schema : bool, optional
        If :py:const:`True`, attempt to extract information from the XML schema
        mentioned in the mzIdentML header (default). Otherwise, use default
        parameters. Disable this to avoid waiting on long network connections or
        if you don't like to get the related warnings.

    Returns
    -------
    out : iterator
       An iterator over the dicts with PSM properties.
    """

    return iterfind(source, 'SpectrumIdentificationResult', **kwargs)

chain = aux._make_chain(read, 'read')

def is_decoy(psm):
    """Given a PSM dict, return :py:const:`True` if all proteins in the dict
    are marked as decoy, and :py:const:`False` otherwise.

    Parameters
    ----------
    psm : dict
        A dict, as yielded by :py:func:`read`.

    Returns
    -------
    out : bool
    """
    return all(pe['isDecoy'] for sii in psm['SpectrumIdentificationItem']
            for pe in sii['PeptideEvidenceRef'])

fdr = aux._make_fdr(is_decoy)
filter = aux._make_filter(chain, is_decoy, lambda x: min(
    sii['mascot:expectation value'] for sii in x['SpectrumIdentificationItem']))
filter.chain = aux._make_chain(filter, 'filter')
