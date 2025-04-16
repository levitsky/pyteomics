"""
idxml - idXML file reader
=========================

Summary
-------

**idXML** is a format specified in the
`OpenMS <http://open-ms.sourceforge.net/about/>`_ project.
It defines a list of peptide identifications.

This module provides a minimalistic way to extract information from idXML
files. You can use the old functional interface (:py:func:`read`) or the new
object-oriented interface (:py:class:`IDXML`) to iterate over entries in
``<PeptideIdentification>`` elements. Note that each entry can contain more than one PSM
(peptide-spectrum match). They are accessible with ``'PeptideHit'`` key.
:py:class:`IDXML` objects also support direct indexing by element ID.

Data access
-----------

  :py:class:`IDXML` - a class representing a single idXML file.
  Other data access functions use this class internally.

  :py:func:`read` - iterate through peptide-spectrum matches in an idXML
  file. Data from a single PSM group are converted to a human-readable dict.
  Basically creates an :py:class:`IDXML` object and reads it.

  :py:func:`chain` - read multiple files at once.

  :py:func:`chain.from_iterable` - read multiple files at once, using an
  iterable of files.

  :py:func:`DataFrame` - read idXML files into a :py:class:`pandas.DataFrame`.

Target-decoy approach
---------------------

  :py:func:`filter` - read a chain of idXML files and filter to a certain
  FDR using TDA.

  :py:func:`filter.chain` - chain a series of filters applied independently to
  several files.

  :py:func:`filter.chain.from_iterable` - chain a series of filters applied
  independently to an iterable of files.

  :py:func:`filter_df` - filter idXML files and return a :py:class:`pandas.DataFrame`.

  :py:func:`is_decoy` - determine if a "SpectrumIdentificationResult" should be
  consiudered decoy.

  :py:func:`fdr` - estimate the false discovery rate of a set of identifications
  using the target-decoy approach.

  :py:func:`qvalues` - get an array of scores and local FDR values for a PSM
  set using the target-decoy approach.

Deprecated functions
--------------------

  :py:func:`version_info` - get information about idXML version and schema.
  You can just read the corresponding attribute of the :py:class:`IDXML`
  object.

  :py:func:`get_by_id` - get an element by its ID and extract the data from it.
  You can just call the corresponding method of the :py:class:`IDXML`
  object.

  :py:func:`iterfind` - iterate over elements in an idXML file.
  You can just call the corresponding method of the :py:class:`IDXML`
  object.

Dependencies
------------

This module requires :py:mod:`lxml`.

-------------------------------------------------------------------------------
"""

#   Copyright 2020 Lev Levitsky
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
from .. import auxiliary as aux
from .. import xml, _schema_defaults


class IDXML(xml.IndexedXML):
    """Parser class for idXML files."""
    file_format = 'idXML'
    _root_element = 'IdXML'
    _default_schema = _schema_defaults._idxml_schema_defaults
    _default_version = '1.5'
    _default_iter_tag = 'PeptideIdentification'
    _structures_to_flatten = {}
    _indexed_tags = {'ProteinHit'}
    _schema_location_param = 'noNamespaceSchemaLocation'

    def __init__(self, *args, **kwargs):
        kwargs.setdefault('retrieve_refs', True)
        super(IDXML, self).__init__(*args, **kwargs)

    def _get_info_smart(self, element, **kwargs):
        """Extract the info in a smart way depending on the element type"""
        name = xml._local_name(element)
        kwargs = dict(kwargs)
        rec = kwargs.pop("recursive", None)

        # Try not to recursively unpack the root element
        # unless the user really wants to.
        if name == self._root_element:
            info = self._get_info(element, recursive=(rec if rec is not None else False), **kwargs)
        else:
            info = self._get_info(element, recursive=(rec if rec is not None else True), **kwargs)
        for k in ['start', 'end']:
            v = info.get(k)
            if isinstance(v, list) and len(v) == 2:
                info[k] = [int(x) for x in v[0].split()]
        for k in ['aa_before', 'aa_after']:
            if k in info:
                info[k] = info[k].split()
        return info

    def _retrieve_refs(self, info, **kwargs):
        """Retrieves and embeds the data for each attribute in `info` that
        ends in _ref. Removes the id attribute from `info`"""
        for k, v in dict(info).items():
            if k[-5:] == '_refs':
                try:
                    by_id = [self.get_by_id(x, retrieve_refs=True) for x in v.split()]
                except KeyError:
                    warnings.warn('Ignoring unresolved reference: ' + v)
                else:
                    for x in by_id:
                        x.pop('id', None)
                    info[k[:-5]] = by_id
                    del info[k]


def read(source, **kwargs):
    """Parse `source` and iterate through peptide-spectrum matches.

    .. note:: This function is provided for backward compatibility only.
        It simply creates an :py:class:`IDXML` instance using
        provided arguments and returns it.

    Parameters
    ----------
    source : str or file
        A path to a target IDXML file or the file object itself.

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
        mentioned in the IDXML header (default). Otherwise, use default
        parameters. Disable this to avoid waiting on slow network connections or
        if you don't like to get the related warnings.

    build_id_cache : bool, optional
        Defines whether a cache of element IDs should be built and stored on the
        created :py:class:`IDXML` instance. Default value is the value of
        `retrieve_refs`.

        .. note:: This parameter is ignored when ``use_index`` is ``True`` (default).

    use_index : bool, optional
        Defines whether an index of byte offsets needs to be created for
        the indexed elements. If :py:const:`True` (default), `build_id_cache` is ignored.

    indexed_tags : container of bytes, optional
        Defines which elements need to be indexed. Empty set by default.

    Returns
    -------
    out : IDXML
       An iterator over the dicts with PSM properties.
    """
    kwargs = kwargs.copy()
    kwargs.setdefault('retrieve_refs', True)
    kwargs['build_id_cache'] = kwargs.get('build_id_cache', kwargs.get('retrieve_refs'))
    return IDXML(source, **kwargs)


def iterfind(source, path, **kwargs):
    """Parse `source` and yield info on elements with specified local
    name or by specified "XPath".

    .. note:: This function is provided for backward compatibility only.
        If you do multiple :py:func:`iterfind` calls on one file, you should
        create an :py:class:`IDXML` object and use its
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
        mentioned in the IDXML header (default). Otherwise, use default
        parameters. Disable this to avoid waiting on slow network connections or
        if you don't like to get the related warnings.

    build_id_cache : bool, optional
        Defines whether a cache of element IDs should be built and stored on the
        created :py:class:`IDXML` instance. Default value is the value of
        `retrieve_refs`.

    Returns
    -------
    out : iterator
    """
    kwargs = kwargs.copy()
    kwargs['build_id_cache'] = kwargs.get('build_id_cache', kwargs.get('retrieve_refs'))
    return IDXML(source, **kwargs).iterfind(path, **kwargs)


version_info = xml._make_version_info(IDXML)


def get_by_id(source, elem_id, **kwargs):
    """Parse `source` and return the element with `id` attribute equal
    to `elem_id`. Returns :py:const:`None` if no such element is found.

    .. note:: This function is provided for backward compatibility only.
        If you do multiple :py:func:`get_by_id` calls on one file, you should
        create an :py:class:`IDXML` object and use its
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
    return IDXML(source, **kwargs).get_by_id(elem_id, **kwargs)


chain = aux.ChainBase._make_chain(IDXML)


def is_decoy(psm, prefix=None):
    """Given a PSM dict, return :py:const:`True` if it is marked as decoy,
    and :py:const:`False` otherwise.

    Parameters
    ----------
    psm : dict
        A dict, as yielded by :py:func:`read`.
    prefix : ignored

    Returns
    -------
    out : bool
    """
    return psm['PeptideHit'][0]['target_decoy'] == 'decoy'


def DataFrame(*args, **kwargs):
    """Read idXML files into a :py:class:`pandas.DataFrame`.

    Requires :py:mod:`pandas`.

    .. warning :: Only the first 'PeptideHit' element is considered in every 'PeptideIdentification'.

    Parameters
    ----------
    *args
        Passed to :py:func:`chain`

    **kwargs
        Passed to :py:func:`chain`

    sep : str or None, keyword only, optional
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
            peptide_hit = item.get('PeptideHit', [None])[0]
            if peptide_hit is not None:
                info.update((k, v) for k, v in peptide_hit.items() if isinstance(v, (str, int, float)))
                protein = peptide_hit.get('protein')
                if protein:
                    accessions, isd, starts, ends, scores, aa_bs, aa_as = [], [], [], [], [], [], []
                    for d, start, end, aab, aaa in zip(protein, peptide_hit['start'], peptide_hit['end'], peptide_hit['aa_before'], peptide_hit['aa_after']):
                        accessions.append(d.get('accession'))
                        isd.append(d.get('target_decoy'))
                        scores.append(d.get('score'))
                        starts.append(start)
                        ends.append(end)
                        aa_bs.append(aab)
                        aa_as.append(aaa)

                    isd = all(x == 'decoy' for x in isd)
                    if sep is not None:
                        if all(isinstance(acc, str) for acc in accessions):
                            accessions = sep.join(accessions)
                        if all(isinstance(aaa, str) for aaa in aa_as):
                            aa_as = sep.join(aa_as)
                        if all(isinstance(aab, str) for aab in aa_bs):
                            aa_bs = sep.join(aa_bs)
                    if all(acc is None for acc in accessions):
                        accessions = None

                    info.update((k, v) for k, v in protein[0].items() if isinstance(v, (str, int, float, list)))
                    info['accession'] = accessions
                    info['is decoy'] = isd
                    info['start'] = starts
                    info['end'] = ends
                    info['aa_before'] = aa_bs
                    info['aa_after'] = aa_as
            data.append(info)
    df = pd.DataFrame(data)
    return df


def filter_df(*args, **kwargs):
    """Read idXML files or DataFrames and return a :py:class:`DataFrame` with filtered PSMs.
    Positional arguments can be idXML files or DataFrames.

    Requires :py:mod:`pandas`.

    .. warning :: Only the first 'PeptideHit' element is considered in every 'PeptideIdentification'.

    Parameters
    ----------
    key : str / iterable / callable, keyword only, optional
        Peptide identification score. Default is 'score'. You will probably need to change it.
    is_decoy : str / iterable / callable, keyword only, optional
        Default is 'is decoy'.
    *args
        Passed to :py:func:`auxiliary.filter` and/or :py:func:`DataFrame`.
    **kwargs
        Passed to :py:func:`auxiliary.filter` and/or :py:func:`DataFrame`.

    Returns
    -------
    out : pandas.DataFrame
    """
    import pandas as pd
    kwargs.setdefault('key', 'score')
    if all(isinstance(arg, pd.DataFrame) for arg in args):
        df = pd.concat(args)
    else:
        df = DataFrame(*args, **kwargs)
    if 'is_decoy' not in kwargs:
        kwargs['is_decoy'] = 'is decoy'
    return aux.filter(df, **kwargs)


fdr = aux._make_fdr(is_decoy, None)
_key = lambda x: x['PeptideHit'][0]['score']
qvalues = aux._make_qvalues(chain, is_decoy, None, _key)
filter = aux._make_filter(chain, is_decoy, None, _key, qvalues)
filter.chain = aux._make_chain(filter, 'filter', True)
