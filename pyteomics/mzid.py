"""
mzid - mzIdentML file reader
============================

Summary
-------

`mzIdentML <http://www.psidev.info/mzidentml>`_  is one of the standards
developed by the Proteomics Informatics working group of the HUPO Proteomics
Standard Initiative.

This module provides a minimalistic way to extract information from mzIdentML
files. You can use the old functional interface (:py:func:`read`) or the new
object-oriented interface (:py:class:`MzIdentML`) to iterate over entries in
``<SpectrumIdentificationResult>`` elements, i.e. groups of identifications
for a certain spectrum. Note that each entry can contain more than one PSM
(peptide-spectrum match). They are accessible with "SpectrumIdentificationItem"
key.
:py:class:`MzIdentML` objects also support direct indexing by element ID.

Data access
-----------

  :py:class:`MzIdentML` - a class representing a single MzIdentML file.
  Other data access functions use this class internally.

  :py:func:`read` - iterate through peptide-spectrum matches in an mzIdentML
  file. Data from a single PSM group are converted to a human-readable dict.
  Basically creates an :py:class:`MzIdentML` object and reads it.

  :py:func:`chain` - read multiple files at once.

  :py:func:`chain.from_iterable` - read multiple files at once, using an
  iterable of files.

  :py:func:`DataFrame` - read MzIdentML files into a :py:class:`pandas.DataFrame`.

Target-decoy approach
---------------------

  :py:func:`filter` - read a chain of mzIdentML files and filter to a certain
  FDR using TDA.

  :py:func:`filter.chain` - chain a series of filters applied independently to
  several files.

  :py:func:`filter.chain.from_iterable` - chain a series of filters applied
  independently to an iterable of files.

  :py:func:`filter_df` - filter MzIdentML files and return a :py:class:`pandas.DataFrame`.

  :py:func:`is_decoy` - determine if a "SpectrumIdentificationResult" should be
  consiudered decoy.

  :py:func:`fdr` - estimate the false discovery rate of a set of identifications
  using the target-decoy approach.

  :py:func:`qvalues` - get an array of scores and local FDR values for a PSM
  set using the target-decoy approach.

Deprecated functions
--------------------

  :py:func:`version_info` - get information about mzIdentML version and schema.
  You can just read the corresponding attribute of the :py:class:`MzIdentML`
  object.

  :py:func:`get_by_id` - get an element by its ID and extract the data from it.
  You can just call the corresponding method of the :py:class:`MzIdentML`
  object.

  :py:func:`iterfind` - iterate over elements in an mzIdentML file.
  You can just call the corresponding method of the :py:class:`MzIdentML`
  object.

Dependencies
------------

This module requires :py:mod:`lxml`.

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

import warnings
warnings.formatwarning = lambda msg, *args: str(msg) + '\n'
from . import auxiliary as aux
from . import xml, _schema_defaults

class MzIdentML(xml.IndexSavingXML):
    """Parser class for MzIdentML files."""
    file_format = 'mzIdentML'
    _root_element = 'MzIdentML'
    _default_schema = _schema_defaults._mzid_schema_defaults
    _default_version = '1.1.0'
    _default_iter_tag = 'SpectrumIdentificationResult'
    _structures_to_flatten = {'Fragmentation'}
    _indexed_tags = {
        'PeptideEvidence', 'SpectrumIdentificationItem', 'SearchDatabase',
        'DBSequence', 'SpectraData', 'Peptide'}

    def __init__(self, *args, **kwargs):
        kwargs.setdefault('retrieve_refs', True)
        super(MzIdentML, self).__init__(*args, **kwargs)

    def _get_info_smart(self, element, **kwargs):
        """Extract the info in a smart way depending on the element type"""
        name = xml._local_name(element)
        kwargs = dict(kwargs)
        rec = kwargs.pop("recursive", None)

        # Try not to recursively unpack the root element
        # unless the user really wants to.
        if name == self._root_element:
            return self._get_info(element,
                    recursive=(rec if rec is not None else False),
                    **kwargs)
        else:
            return self._get_info(element,
                    recursive=(rec if rec is not None else True),
                    **kwargs)

    def _retrieve_refs(self, info, **kwargs):
        """Retrieves and embeds the data for each attribute in `info` that
        ends in _ref. Removes the id attribute from `info`"""
        for k, v in dict(info).items():
            if k.endswith('_ref'):
                by_id = self.get_by_id(v, retrieve_refs=True)
                if by_id is None:
                    warnings.warn('Ignoring unresolved reference: ' + v)
                else:
                    info.update(by_id)
                    del info[k]
                    info.pop('id', None)

def read(source, **kwargs):
    """Parse `source` and iterate through peptide-spectrum matches.

    .. note:: This function is provided for backward compatibility only.
        It simply creates an :py:class:`MzIdentML` instance using
        provided arguments and returns it.

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
        increase. Default is :py:const:`True`.

    iterative : bool, optional
        Specifies whether iterative XML parsing should be used. Iterative
        parsing significantly reduces memory usage and may be just a little
        slower. When `retrieve_refs` is :py:const:`True`, however, it is
        highly recommended to disable iterative parsing if possible.
        Default value is :py:const:`True`.

    read_schema : bool, optional
        If :py:const:`True`, attempt to extract information from the XML schema
        mentioned in the mzIdentML header (default). Otherwise, use default
        parameters. Disable this to avoid waiting on slow network connections or
        if you don't like to get the related warnings.

    build_id_cache : bool, optional
        Defines whether a cache of element IDs should be built and stored on the
        created :py:class:`MzIdentML` instance. Default value is the value of
        `retrieve_refs`.

        .. note:: This parameter is ignored when ``use_index`` is ``True`` (default).

    use_index : bool, optional
        Defines whether an index of byte offsets needs to be created for
        the indexed elements. If :py:const:`True` (default), `build_id_cache` is ignored.

    indexed_tags : container of bytes, optional
        Defines which elements need to be indexed. Empty set by default.

    Returns
    -------
    out : MzIdentML
       An iterator over the dicts with PSM properties.
    """
    kwargs = kwargs.copy()
    kwargs.setdefault('retrieve_refs', True)
    kwargs['build_id_cache'] = kwargs.get('build_id_cache', kwargs.get('retrieve_refs'))
    return MzIdentML(source, **kwargs)

def iterfind(source, path, **kwargs):
    """Parse `source` and yield info on elements with specified local
    name or by specified "XPath".

    .. note:: This function is provided for backward compatibility only.
        If you do multiple :py:func:`iterfind` calls on one file, you should
        create an :py:class:`MzIdentML` object and use its
        :py:meth:`!iterfind` method.

    Parameters
    ----------
    source : str or file
        File name or file-like object.

    path : str
        Element name or XPath-like expression. Only local names separated
        with slashes are accepted. An asterisk (`*`) means any element.
        You can specify a single condition in the end, such as:
        ``"/path/to/element[some_value>1.5]"``
        Note: you can do much more powerful filtering using plain Python.
        The path can be absolute or "free". Please don't specify
        namespaces.

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
        Default value is :py:const:`True`.

    read_schema : bool, optional
        If :py:const:`True`, attempt to extract information from the XML schema
        mentioned in the mzIdentML header (default). Otherwise, use default
        parameters. Disable this to avoid waiting on slow network connections or
        if you don't like to get the related warnings.

    build_id_cache : bool, optional
        Defines whether a cache of element IDs should be built and stored on the
        created :py:class:`MzIdentML` instance. Default value is the value of
        `retrieve_refs`.

    Returns
    -------
    out : iterator
    """
    kwargs = kwargs.copy()
    kwargs['build_id_cache'] = kwargs.get('build_id_cache',
            kwargs.get('retrieve_refs'))
    return MzIdentML(source, **kwargs).iterfind(path, **kwargs)

version_info = xml._make_version_info(MzIdentML)

def get_by_id(source, elem_id, **kwargs):
    """Parse `source` and return the element with `id` attribute equal
    to `elem_id`. Returns :py:const:`None` if no such element is found.

    .. note:: This function is provided for backward compatibility only.
        If you do multiple :py:func:`get_by_id` calls on one file, you should
        create an :py:class:`MzIdentML` object and use its
        :py:meth:`!get_by_id` method.

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
    return MzIdentML(source, **kwargs).get_by_id(elem_id, **kwargs)

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

def DataFrame(*args, **kwargs):
    """Read MzIdentML files into a :py:class:`pandas.DataFrame`.

    Requires :py:mod:`pandas`.

    .. warning :: Only the first 'SpectrumIdentificationItem' element is considered in every
                  'SpectrumIdentificationResult'.

    Parameters
    ----------
    *args, **kwargs : passed to :py:func:`chain`

    sep : str or None, optional
        Some values related to PSMs (such as protein information) are variable-length
        lists. If `sep` is a :py:class:`str`, they will be packed into single string using
        this delimiter. If `sep` is :py:const:`None`, they are kept as lists. Default is
        :py:const:`None`.

    Returns
    -------
    out : pandas.DataFrame
    """
    import pandas as pd
    data = []

    sep = kwargs.pop('sep', None)
    with chain(*args, **kwargs) as f:
        for item in f:
            info = {}
            for k, v in item.items():
                if isinstance(v, (str, int, float)):
                    info[k] = v
            sii = item.get('SpectrumIdentificationItem', [None])[0]
            if sii is not None:
                info.update((k, v) for k, v in sii.items() if isinstance(v, (str, int, float)))
                evref = sii.get('PeptideEvidenceRef')
                if evref:
                    prot_descr, accessions, isd, starts, ends, lengths = [], [], [], [], [], []
                    for d in evref:
                        prot_descr.append(d.get('protein description'))
                        accessions.append(d.get('accession'))
                        isd.append(d.get('isDecoy'))
                        starts.append(d.get('start'))
                        ends.append(d.get('end'))
                        lengths.append(d.get('length'))
                    isd = all(isd)
                    if sep is not None:
                        if all(isinstance(prd, str) for prd in prot_descr):
                            prot_descr = sep.join(prot_descr)

                        if all(isinstance(acc, str) for acc in accessions):
                            accessions = sep.join(accessions)

                    if all(prd is None for prd in prot_descr):
                        prot_descr = None
                    if all(acc is None for acc in accessions):
                        accessions = None

                    info.update((k, v) for k, v in evref[0].items() if isinstance(v, (str, int, float, list)))
                    info['protein description'] = prot_descr
                    info['accession'] = accessions
                    info['isDecoy'] = isd
                    info['start'] = starts
                    info['end'] = ends
                    info['length'] = lengths
            data.append(info)
    df = pd.DataFrame(data)
    return df

def filter_df(*args, **kwargs):
    """Read MzIdentML files or DataFrames and return a :py:class:`DataFrame` with filtered PSMs.
    Positional arguments can be MzIdentML files or DataFrames.

    Requires :py:mod:`pandas`.

    .. warning :: Only the first 'SpectrumIdentificationItem' element is considered in every
                  'SpectrumIdentificationResult'.

    Parameters
    ----------
    key : str / iterable / callable, optional
        Default is 'mascot:expectation value'.
    is_decoy : str / iterable / callable, optional
        Default is 'isDecoy'.
    *args, **kwargs : passed to :py:func:`auxiliary.filter` and/or :py:func:`DataFrame`.

    Returns
    -------
    out : pandas.DataFrame
    """
    import pandas as pd
    kwargs.setdefault('key', 'mascot:expectation value')
    kwargs.setdefault('is_decoy', 'isDecoy')
    if all(isinstance(arg, pd.DataFrame) for arg in args):
        df = pd.concat(args)
    else:
        df = DataFrame(*args, **kwargs)
    return aux.filter(df, **kwargs)

fdr = aux._make_fdr(is_decoy, None)
_key = lambda x: min(
    sii['mascot:expectation value'] for sii in x['SpectrumIdentificationItem'])
qvalues = aux._make_qvalues(chain, is_decoy, None, _key)
filter = aux._make_filter(chain, is_decoy, None, _key, qvalues)
filter.chain = aux._make_chain(filter, 'filter', True)