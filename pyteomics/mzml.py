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
import zlib
import base64
import re
from . import xml, auxiliary as aux
from .xml import etree

def _decode_base64_data_array(source, dtype, is_compressed):
    """Read a base64-encoded binary array.

    Parameters
    ----------
    source : str
        A binary array encoded with base64.
    dtype : str
        The type of the array in numpy dtype notation.
    is_compressed : bool
        If True then the array will be decompressed with zlib.

    Returns
    -------
    out : numpy.array
    """

    decoded_source = base64.b64decode(source.encode('ascii'))
    if is_compressed:
        decoded_source = zlib.decompress(decoded_source)
    output = np.frombuffer(decoded_source, dtype=dtype)
    return output

class MzML(xml.IndexedXML):
    """Parser class for mzML files."""
    file_format = 'mzML'
    _root_element = 'mzML'
    _default_schema = xml._mzml_schema_defaults
    _default_version = '1.1.0'
    _default_iter_tag = 'spectrum'
    _structures_to_flatten = {'binaryDataArrayList'}
    _indexed_tags = {'spectrum', 'chromatogram'}

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
        if 'binary' in info:
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
            compressed = True
            if 'zlib compression' in info:
                del info['zlib compression']
            elif 'name' in info and 'zlib compression' in info['name']:
                info['name'].remove('zlib compression')
            else:
                compressed = False
                info.pop('no compression', None)
                try:
                    info['name'].remove('no compression')
                    if not info['name']: del info['name']
                except (KeyError, TypeError):
                    pass
            b = info.pop('binary')
            if b:
                array = _decode_base64_data_array(
                                b, dtype, compressed)
            else:
                array = np.array([], dtype=dtype)
            for k in info:
                if k.endswith(' array') and not info[k]:
                    info = {k: array}
                    break
            else:
                found = False
                # workaround for https://bitbucket.org/levitsky/pyteomics/issues/11
                if isinstance(info.get('name'), list):
                    knames = info['name']
                    for val in knames:
                        if val.endswith(' array'):
                            info = {val: array}
                            found = True
                            break
                # last fallback
                if not found:
                    info['binary'] = array
        if 'binaryDataArray' in info:
            for array in info.pop('binaryDataArray'):
                info.update(array)
        intkeys = {'ms level'}
        for k in intkeys:
            if k in info:
                info[k] = int(info[k])
        return info

def read(source, read_schema=True, iterative=True, use_index=False):
    """Parse `source` and iterate through spectra.

    Parameters
    ----------
    source : str or file
        A path to a target mzML file or the file object itself.

    read_schema : bool, optional
        If :py:const:`True`, attempt to extract information from the XML schema
        mentioned in the mzML header (default). Otherwise, use default
        parameters. Disable this to avoid waiting on long network connections or
        if you don't like to get the related warnings.

    iterative : bool, optional
        Defines whether iterative parsing should be used. It helps reduce
        memory usage at almost the same parsing speed. Default is
        :py:const:`True`.

    use_index : bool, optional
        Defines whether an index of byte offsets needs to be created for
        spectrum elements. Default is :py:const:`False`.

    Returns
    -------
    out : iterator
       An iterator over the dicts with spectra properties.
    """

    return MzML(source, read_schema=read_schema, iterative=iterative, use_index=use_index)

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
        mentioned in the mzIdentML header (default). Otherwise, use default
        parameters. Disable this to avoid waiting on long network connections or
        if you don't like to get the related warnings.

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