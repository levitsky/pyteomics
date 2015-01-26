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
from . import auxiliary as aux
from . import xml

class MzIdentMLParser(xml.XMLParserBase):
    file_format = "mzIdentML"
    root_element = "MzIdentML"
    version_info_element = "MzIdentML"
    default_schema = xml._mzid_schema_defaults
    default_version = "1.1.0"
    default_iter_tag = "SpectrumIdentificationResult"
    structures_to_flatten = {'Fragmentation'}

    def __init__(self, source, **kwargs):
        super(MzIdentMLParser, self).__init__(source, **kwargs)
        if kwargs.pop('build_id_cache', False):
            self.build_id_cache()
        else:
            self.id_dict = {}

    def get_info_smart(self, element, **kwargs):
        """Extract the info in a smart way depending on the element type"""
        name = xml._local_name(element)
        kwargs = dict(kwargs)
        rec = kwargs.pop("recursive", None)

        # Try not to recursively unpack the root element
        # unless the user really wants to.
        if name == self.root_element:
            return self.get_info(element, recursive=(rec if rec is not None else False),
                                 **kwargs)
        else:
            return self.get_info(element, recursive=(rec if rec is not None else True),
                                 **kwargs)

    def _retrieve_refs(self, info, **kwargs):
        '''Retrieves and embeds the data for each attribute in `info` that
        ends in _ref. Removes the id attribute from `info`'''
        for k, v in dict(info).items():
            if k.endswith('_ref'):
                info.update(self.get_by_id(v, retrieve_refs=True))
                del info[k]
                info.pop('id', None)

    @xml._keepstate
    def build_id_cache(self):
        '''Constructs a cache for each element in the document, indexed by id
        attribute'''
        stack = 0
        id_dict = {}
        for event, elem in etree.iterparse(self.source, events=('start', 'end'),
                remove_comments=True):
            if event == 'start':
                if 'id' in elem.attrib:
                    stack += 1
            else:
                if 'id' in elem.attrib:
                    stack -= 1
                    id_dict[elem.attrib['id']] = elem
                elif stack == 0:
                    elem.clear()
        self.id_dict = id_dict

    def clear_id_cache(self):
        self.id_dict = {}

    @xml._keepstate
    def get_by_id(self, elem_id, **kwargs):
        """Parse `self.source` and return the element with `id` attribute equal
        to `elem_id`. Returns :py:const:`None` if no such element is found.

        Parameters
        ----------
        source : str or file
            A path to a target mzIdentML file of the file object itself.

        elem_id : str
            The value of the `id` attribute to match.

        Returns
        -------
        out : :py:class:`dict` or :py:const:`None`
        """
        if kwargs.get('iterative'):
            found = False
            for event, elem in etree.iterparse(self.source,
                    events=('start', 'end'), remove_comments=True):
                if event == 'start':
                    if elem.attrib.get('id') == elem_id:
                        found = True
                else:
                    if elem.attrib.get('id') == elem_id:
                        return self.get_info_smart(elem, retrieve_refs=True)
                    if not found:
                        elem.clear()
            return None
        # Otherwise do build and use the id_dict to cache elements
        else:
            self._build_id_cache()
            return self.get_info_smart(self.id_dict[elem_id], **kwargs)


def read(source, **kwargs):
    """Parse `source` and iterate through peptide-spectrum matches.

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

    iterative : bool, optional
        Specifies whether iterative XML parsing should be used. Iterative
        parsing significantly reduces memory usage and may be just a little
        slower. When `retrieve_refs` is :py:const:`True`, however, it is
        highly recommended to disable iterative parsing if possible.
        Default value is the opposite of `retrieve_refs`.

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
    return MzIdentMLParser(source, **kwargs)

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
_key = lambda x: min(
    sii['mascot:expectation value'] for sii in x['SpectrumIdentificationItem'])
local_fdr = aux._make_local_fdr(chain, is_decoy, _key)
filter = aux._make_filter(chain, is_decoy, _key, local_fdr)
filter.chain = aux._make_chain(filter, 'filter')
