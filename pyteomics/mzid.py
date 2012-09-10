"""
mzidentml - reader for peptide identification data in mzIdentML format
======================================================================

Summary
-------

`mzIdentML <http://www.psidev.info/mzidentml>`_  is one of the standards
developed by the Proteomics Informatics working group of the HUPO Proteomics
Standard Initiative.

This module provides a minimalistic way to extract information from mzIdentML
files. The main idea is the same as in :py:mod:`pyteomics.pepxml`: the top-level
function :py:func:`read` allows iterating over entries in
`<SpectrumIdentificationResult>` tags, i.e. groups of identifications
for a certain spectrum. Note that each entry can contain more than one PSM
(peptide-spectrum match). They are accessible with `"search_hits"` key, just
like in :py:mod:`pyteomics.pepxml`.

Data access
-----------

  :py:func:`read` - iterate through peptide-spectrum matches in a pep.XML 
  file. Data from a single PSM group are converted to a human-readable dict. 

  :py:func:`get_by_id` - get an element by its ID and extract the data from it.

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

def _local_name(element):
    if element.tag.startswith('{'):
        return element.tag.rsplit('}', 1)[1]
    else:
        return element.tag

def _get_info(element, recursive=False):
    """Extract info from element's attributes, possibly recursive.
    <cvParam> elements are treated in a special way."""

    if _local_name(element) == 'cvParam':
        if 'name' in element.attrib:
            if 'value' in element.attrib:
                return {element.attrib['name']: element.attrib['value']}
            else:
                return element.attrib['name']
        else:
            return dict(element.attrib) 
    info = dict(element.attrib)
    if recursive:
        for child in element.iterchildren():
            info[_local_name(child)] = _get_info(child, True)
    return info

def _get_info_smart(element):
    """Extract the info in a smart way depending on the element type"""

    raise NotImplementedError

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
    for _, e in etree.iterparse(source):
        if e.attrib.get('id') == elem_id:
            return _get_info_smart(e)
    return None

def read(source):
    """Parse ``source`` and iterate through peptide-spectrum matches.

    Parameters
    ----------
    source : str or file
        A path to a target mzIdentML file or the file object itself.

    Returns
    -------
    out : iterator
       An iterator over the dicts with PSM properties.
    """

    for _, elem in etree.iterparse(source):
        if elem.tag == '{{{}}}SpectrumIdentificationResult'.format(
                elem.nsmap.get(None, '')):
            yield _get_info_smart(elem)
            elem.clear()
