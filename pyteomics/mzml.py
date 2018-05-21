"""
mzml - reader for mass spectrometry data in mzML format
=======================================================

Summary
-------

mzML is a standard rich XML-format for raw mass spectrometry data storage.
Please refer to `psidev.info <http://www.psidev.info/index.php?q=node/257>`_
for the detailed specification of the format and structure of mzML files.

This module provides a minimalistic way to extract information from mzML
files. You can use the old functional interface (:py:func:`read`) or the new
object-oriented interface (:py:class:`MzML` or :py:class:`PreIndexedMzML`)
to iterate over entries in ``<spectrum>`` elements.
:py:class:`MzML` and :py:class:`PreIndexedMzML` also support direct indexing
with spectrum IDs.

Data access
-----------

  :py:class:`MzML` - a class representing a single mzML file.
  Other data access functions use this class internally.

  :py:class:`PreIndexedMzML` - a class representing a single mzML file.
  Uses byte offsets listed at the end of the file for quick access to spectrum elements.

  :py:func:`read` - iterate through spectra in mzML file. Data from a
  single spectrum are converted to a human-readable dict. Spectra themselves are
  stored under 'm/z array' and 'intensity array' keys.

  :py:func:`chain` - read multiple mzML files at once.

  :py:func:`chain.from_iterable` - read multiple files at once, using an
  iterable of files.

Deprecated functions
--------------------

  :py:func:`version_info` - get version information about the mzML file.
  You can just read the corresponding attribute of the :py:class:`MzML` object.

  :py:func:`iterfind` - iterate over elements in an mzML file.
  You can just call the corresponding method of the :py:class:`MzML` object.

Dependencies
------------

This module requires :py:mod:`lxml` and :py:mod:`numpy`.

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
import re
import warnings
from . import xml, auxiliary as aux, _schema_defaults
from .xml import etree

NON_STANDARD_DATA_ARRAY = 'non-standard data array'

STANDARD_ARRAYS = set([
    'm/z array',
    'intensity array',
    'charge array',
    'signal to noise array',
    'time array',
    'wavelength array',
    'flow rate array',
    'pressure array',
    'temperature array',
    'mean drift time array',
    'mean charge array',
    'resolution array',
    'baseline array'
])


class MzML(xml.ArrayConversionMixin, xml.IndexSavingXML):
    """Parser class for mzML files."""
    file_format = 'mzML'
    _root_element = 'mzML'
    _default_schema = _schema_defaults._mzml_schema_defaults
    _default_version = '1.1.0'
    _default_iter_tag = 'spectrum'
    _structures_to_flatten = {'binaryDataArrayList', "referenceableParamGroupRef"}
    _indexed_tags = {'spectrum', 'chromatogram'}

    def __init__(self, *args, **kwargs):
        self.decode_binary = kwargs.pop('decode_binary', True)
        xml.IndexSavingXML.__init__(self, *args, **kwargs)
        xml.ArrayConversionMixin.__init__(self, *args, **kwargs)

    def _detect_array_name(self, info):
        """Determine what the appropriate name for this
        array is by inspecting the available param-based
        keys.

        Parameters
        ----------
        info : dict
            The collapsed binary tag plus
            associated *Param  data

        Returns
        -------
        out : str
            The name for this array entry
        """
        # If this is a non-standard array, we hope the userParams
        # will conform to the same array suffix pattern.
        is_non_standard = False

        # Accumulate possible name candidates
        candidates = []
        for k in info:
            if k.endswith(' array') and not info[k]:
                if NON_STANDARD_DATA_ARRAY == k:
                    is_non_standard = True
                else:
                    candidates.append(k)
        if isinstance(info.get('name'), list):
            for val in info['name']:
                if val.endswith(' array'):
                    if NON_STANDARD_DATA_ARRAY == val:
                        is_non_standard = True
                    else:
                        candidates.append(val)

        # Name candidate resolution
        n_candidates = len(candidates)
        # Easy case, exactly one name given
        if n_candidates == 1:
            return candidates[0]
        # We are missing information, but at least
        # if we know the array is non-standard we
        # can report it as such. Otherwise fall back
        # to "binary". This fallback signals special
        # behavior elsewhere.
        elif n_candidates == 0:
            if is_non_standard:
                return NON_STANDARD_DATA_ARRAY
            else:
                return "binary"
        # Multiple choices means we need to make a decision which could
        # mask data from the user. This should never happen but stay safe.
        # There are multiple options to choose from. There is no way to
        # make a good choice here. We first prefer the standardized
        # arrays before falling back to just guessing.
        else:
            warnings.warn("Multiple options for naming binary array: %r" % candidates)
            standard_options = set(candidates) & STANDARD_ARRAYS
            if standard_options:
                return max(standard_options, key=len)
            else:
                return max(candidates, key=len)

    def _determine_array_dtype(self, info):
        dtype = None
        types = {'32-bit float': 'f', '64-bit float': 'd'}
        for t, code in types.items():
            if t in info:
                dtype = code
                del info[t]
                break
        # sometimes it's under 'name'
        else:
            if 'name' in info:
                for t, code in types.items():
                    if t in info['name']:
                        dtype = code
                        info['name'].remove(t)
                        break
        return dtype

    def _determine_compression(self, info):
        known_compression_types = set(self.compression_type_map)
        found_compression_types = known_compression_types & set(info)
        if found_compression_types:
            found_compression_types = tuple(found_compression_types)
            if len(found_compression_types) == 1:
                del info[found_compression_types[0]]
                return found_compression_types[0]
            else:
                warnings.warn("Multiple options for binary array compression: %r" % (
                    found_compression_types,))
                return found_compression_types[0]
        elif "name" in info:
            found_compression_types = known_compression_types & set(info['name'])
            if found_compression_types:
                found_compression_types = tuple(found_compression_types)
                if len(found_compression_types) == 1:
                    if not self._skip_empty_cvparam_values:
                        del info['name'][found_compression_types[0]]
                    return found_compression_types[0]
                else:
                    warnings.warn("Multiple options for binary array compression: %r" % (
                        found_compression_types,))
                    return found_compression_types[0]
        else:
            return 'no compression'

    def _handle_binary(self, info, **kwargs):
        """Special handling when processing and flattening
        a <binary> tag and its sibling *Param tags.

        Parameters
        ----------
        info : dict
            Unprocessed binary array data and metadata

        Returns
        -------
        out : dict
            The processed and flattened data array and metadata
        """
        dtype = self._determine_array_dtype(info)
        compressed = self._determine_compression(info)
        name = self._detect_array_name(info)
        binary = info.pop('binary')
        if not self.decode_binary:
            info[name] = self._make_record(binary, compressed, dtype, name)
            return info

        if binary:
            array = self.decode_data_array(binary, compressed, dtype)
        else:
            array = np.array([], dtype=dtype)

        if name == 'binary':
            info[name] = self._convert_array(None, array)
        else:
            info = {name: self._convert_array(name, array)}
        return info

    def _get_info_smart(self, element, **kw):
        name = xml._local_name(element)
        kwargs = dict(kw)
        rec = kwargs.pop('recursive', None)
        if name in {'indexedmzML', 'mzML'}:
            info =  self._get_info(element,
                    recursive=(rec if rec is not None else False),
                    **kwargs)
        else:
            info = self._get_info(element,
                    recursive=(rec if rec is not None else True),
                    **kwargs)
        if 'binary' in info and isinstance(info, dict):
            info = self._handle_binary(info, **kwargs)

        if 'binaryDataArray' in info and isinstance(info, dict):
            for array in info.pop('binaryDataArray'):
                info.update(array)
        intkeys = {'ms level'}
        for k in intkeys:
            if k in info:
                info[k] = int(info[k])
        return info

    def _retrieve_refs(self, info, **kwargs):
        """Retrieves and embeds the data for each attribute in `info` that
        ends in _ref. Removes the id attribute from `info`"""
        for k, v in dict(info).items():
            if k == 'ref':
                by_id = self.get_by_id(v, retrieve_refs=True)
                if by_id is None:
                    warnings.warn('Ignoring unresolved reference: ' + v)
                else:
                    info.update(by_id)
                    del info[k]
                    info.pop('id', None)


def read(source, read_schema=False, iterative=True, use_index=False, dtype=None, huge_tree=False):
    """Parse `source` and iterate through spectra.

    Parameters
    ----------
    source : str or file
        A path to a target mzML file or the file object itself.

    read_schema : bool, optional
        If :py:const:`True`, attempt to extract information from the XML schema
        mentioned in the mzML header. Otherwise, use default parameters.
        Not recommended without Internet connection or
        if you don't like to get the related warnings.

    iterative : bool, optional
        Defines whether iterative parsing should be used. It helps reduce
        memory usage at almost the same parsing speed. Default is
        :py:const:`True`.

    use_index : bool, optional
        Defines whether an index of byte offsets needs to be created for
        spectrum elements. Default is :py:const:`False`.

    dtype : type or dict, optional
        dtype to convert arrays to, one for both m/z and intensity arrays or one for each key.
        If :py:class:`dict`, keys should be 'm/z array' and 'intensity array'.

    decode_binary : bool, optional
        Defines whether binary data should be decoded and included in the output
        (under "m/z array", "intensity array", etc.).
        Default is :py:const:`True`.

    huge_tree : bool, optional
        This option is passed to the `lxml` parser and defines whether
        security checks for XML tree depth and node size should be disabled.
        Default is :py:const:`False`.
        Enable this option for trusted files to avoid XMLSyntaxError exceptions
        (e.g. `XMLSyntaxError: xmlSAX2Characters: huge text node`).

    Returns
    -------
    out : iterator
       An iterator over the dicts with spectrum properties.
    """

    return MzML(source, read_schema=read_schema, iterative=iterative,
        use_index=use_index, dtype=dtype, huge_tree=huge_tree)

def iterfind(source, path, **kwargs):
    """Parse `source` and yield info on elements with specified local
    name or by specified "XPath".

    .. note:: This function is provided for backward compatibility only.
        If you do multiple :py:func:`iterfind` calls on one file, you should
        create an :py:class:`MzML` object and use its
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

    iterative : bool, optional
        Specifies whether iterative XML parsing should be used. Iterative
        parsing significantly reduces memory usage and may be just a little
        slower. When `retrieve_refs` is :py:const:`True`, however, it is
        highly recommended to disable iterative parsing if possible.
        Default value is :py:const:`True`.

    read_schema : bool, optional
        If :py:const:`True`, attempt to extract information from the XML schema
        mentioned in the mzIdentML header. Otherwise, use default
        parameters. Not recommended without Internet connection or
        if you don't like to get the related warnings.

    decode_binary : bool, optional
        Defines whether binary data should be decoded and included in the output
        (under "m/z array", "intensity array", etc.).
        Default is :py:const:`True`.

    Returns
    -------
    out : iterator
    """
    return MzML(source, **kwargs).iterfind(path, **kwargs)

version_info = xml._make_version_info(MzML)

chain = aux._make_chain(read, 'read')


class PreIndexedMzML(MzML):
    """Parser class for mzML files, subclass of :py:class:`MzML`.
    Uses byte offsets listed at the end of the file for quick access to spectrum elements.
    """
    def _build_index(self):
        """
        Build up a `dict` of `dict` of offsets for elements. Calls :meth:`_find_index_list`
        and assigns the return value to :attr:`_offset_index`
        """
        self._offset_index = xml._flatten_map(self._find_index_list())

    @xml._keepstate
    def _iterparse_index_list(self, offset):
        index_map = {}
        index = {}
        self._source.seek(offset)
        try:
            for event, elem in etree.iterparse(self._source, events=('start', 'end'), remove_comments=True):
                if event == 'start':
                    if elem.tag == 'index':
                        index = {}
                        index_map[elem.attrib['name']] = index
                else:
                    if elem.tag == 'offset':
                        index[elem.attrib['idRef']] = int(elem.text)
                    elem.clear()
        except etree.XMLSyntaxError:
            # The iteration has reached the end of the indexList tag and the parser
            # encounters the later elements in the document.
            pass
        return index_map

    @xml._keepstate
    def _find_index_list_offset(self):
        """
        Search relative to the bottom of the file upwards to find the offsets
        of the index lists.

        Returns
        -------
        list of int
            A list of byte offsets for `<indexList>` elements
        """
        self._source.seek(-1024, 2)
        text = self._source.read(1024)
        index_offsets = list(map(int, re.findall(br'<indexListOffset>(\d+)</indexListOffset>', text)))
        return index_offsets

    @xml._keepstate
    def _find_index_list(self):
        """
        Extract lists of index offsets from the end of the file.

        Returns
        -------
        dict of str -> dict of str -> int
        """
        offsets = self._find_index_list_offset()
        index_list = {}
        for offset in offsets:
            # Sometimes the offset is at the very beginning of the file,
            # due to a bug in an older version of ProteoWizard. If this crude
            # check fails, don't bother searching the entire file, and fall back
            # on the base class's mechanisms.
            #
            # Alternative behavior here would be to start searching for the start
            # of the index from the bottom of the file, but this version of Proteowizard
            # also emits invalid offsets which do not improve retrieval time.
            if offset < 1024:
                continue
            index_list.update(self._iterparse_index_list(offset))
        return index_list
