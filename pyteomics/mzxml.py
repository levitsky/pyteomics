"""
mzxml - reader for mass spectrometry data in mzXML format
=========================================================

Summary
-------

**mzXML** is a (formerly) standard XML-format for raw mass spectrometry data storage,
intended to be replaced with **mzML**.

This module provides a minimalistic way to extract information from mzXML
files. You can use the old functional interface (:py:func:`read`) or the new
object-oriented interface (:py:class:`MzXML`)
to iterate over entries in ``<scan>`` elements.
:py:class:`MzXML` also supports direct indexing with scan IDs.

Data access
-----------

  :py:class:`MzXML` - a class representing a single mzXML file.
  Other data access functions use this class internally.

  :py:func:`read` - iterate through spectra in mzXML file. Data from a
  single scan are converted to a human-readable dict. Spectra themselves are
  stored under 'm/z array' and 'intensity array' keys.

  :py:func:`chain` - read multiple mzXML files at once.

  :py:func:`chain.from_iterable` - read multiple files at once, using an
  iterable of files.

Deprecated functions
--------------------

  :py:func:`version_info` - get version information about the mzXML file.
  You can just read the corresponding attribute of the :py:class:`MzXML` object.

  :py:func:`iterfind` - iterate over elements in an mzXML file.
  You can just call the corresponding method of the :py:class:`MzXML` object.

Dependencies
------------

This module requires :py:mod:`lxml` and :py:mod:`numpy`.

-------------------------------------------------------------------------------
"""

#   Copyright 2016 Joshua Klein, Lev Levitsky
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

import heapq

from . import xml, auxiliary as aux, _schema_defaults
import numpy as np


def _decode_peaks(info, peaks_data):
    """Decode the interleaved base 64 encoded, potentially
    compressed, raw data points.

    Parameters
    ----------
    info : dict
        The current context
    peaks_data : str
        The textually encoded peak data

    Returns
    -------
    tuple of np.array
        A pair of NumPy arrays containing
        m/z and intensity values.
    """
    compressed = (info.get('compressionType') == 'zlib')
    dt = np.float32 if info['precision'] == '32' else np.float64
    dtype = np.dtype([('m/z array', dt), ('intensity array', dt)]).newbyteorder('>')
    data = aux._decode_base64_data_array(peaks_data, dtype, compressed)
    return data


class IteratorQueue(object):
    def __init__(self, iterator):
        q = list()
        heapq.heapify(q)
        self.queue = q
        self.iterator = iterator
        self.last_index = -1
        self.producer = self.consume(iterator)

    def insert_item(self, scan):
        heapq.heappush(self.queue, (int(scan['num']), scan))

    def __iter__(self):
        return self.producer

    def consume(self, iterator):
        for scan in iterator:
            scan.pop("scan", None)
            if scan['msLevel'] != 1:
                self.insert_item(scan)
            else:
                self.insert_item(scan)
                barrier = int(scan['num'])
                while True:
                    idx, item = heapq.heappop(self.queue)
                    if idx >= barrier:
                        self.insert_item(item)
                        break
                    yield item
        while self.queue:
            idx, item = heapq.heappop(self.queue)
            yield item


class MzXML(xml.ArrayConversionMixin, xml.IndexSavingXML):
    """Parser class for mzXML files."""
    _root_element = 'mzXML'
    _default_iter_tag = 'scan'
    _indexed_tags = {'scan'}
    _indexed_tag_keys = {'scan': 'num'}
    _default_version = None
    _default_schema = _schema_defaults._mzxml_schema_defaults
    _default_id_attr = 'num'

    def __init__(self, *args, **kwargs):
        self.decode_binary = kwargs.pop('decode_binary', True)
        xml.IndexSavingXML.__init__(self, *args, **kwargs)
        xml.ArrayConversionMixin.__init__(self, *args, **kwargs)

    def _get_info_smart(self, element, **kw):
        name = xml._local_name(element)

        kwargs = dict(kw)
        rec = kwargs.pop('recursive', None)
        if name in {'mzXML'}:
            info = self._get_info(element,
                                  recursive=(
                                      rec if rec is not None else False),
                                  **kwargs)
        else:
            info = self._get_info(element,
                                  recursive=(rec if rec is not None else True),
                                  **kwargs)
        if 'num' in info and isinstance(info, dict):
            info['id'] = info['num']
        if 'peaks' in info and isinstance(info, dict):
            self._decode_peaks(info)
        return info

    def _determine_compression(self, info):
        if info.get('compressionType') == 'zlib':
            return 'zlib compression'
        return "no compression"

    def _determine_dtype(self, info):
        dt = np.float32 if info['precision'] == '32' else np.float64
        endianess = ">" if info['byteOrder'] in ('network', "big") else "<"
        dtype = np.dtype(
            [('m/z array', dt), ('intensity array', dt)]).newbyteorder(endianess)
        return dtype

    def _finalize_record_conversion(self, array, record):
        key = record.key
        return self._convert_array(key, array[key])

    def _decode_peaks(self, info):
        # handle cases where peaks is the encoded binary data which must be
        # unpacked
        if not isinstance(info['peaks'], (dict, list)):
            compression_type = self._determine_compression(info)
            dtype = self._determine_dtype(info)
            binary = info.pop('peaks')
            if not self.decode_binary:
                for k in self._array_keys:
                    record = self._make_record(binary, compression_type, dtype, k)
                    info[k] = record
            else:
                peak_data = self.decode_data_array(binary, compression_type, dtype)
                for k in self._array_keys:
                    info[k] = self._convert_array(k, peak_data[k])
        # otherwise we've already decoded the arrays and we're just passing
        # them up the hierarchy
        else:
            if not self.decode_binary:
                arrays = info.pop('peaks')[0]
                for k in self._array_keys:
                    info[k] = arrays[k]
            else:
                peak_data = info.pop('peaks')[0]
                for k in self._array_keys:
                    info[k] = self._convert_array(k, peak_data.get(k, np.array([])))

    def iterfind(self, path, **kwargs):
        if path == 'scan':
            generator = super(MzXML, self).iterfind(path, **kwargs)
            for item in IteratorQueue(generator):
                yield item
        else:
            for item in super(MzXML, self).iterfind(path, **kwargs):
                yield item

def read(source, read_schema=False, iterative=True, use_index=False, dtype=None, huge_tree=False):
    """Parse `source` and iterate through spectra.

    Parameters
    ----------
    source : str or file
        A path to a target mzML file or the file object itself.

    read_schema : bool, optional
        If :py:const:`True`, attempt to extract information from the XML schema
        mentioned in the mzML header. Otherwise, use default
        parameters. Not recommended without Internet connection or
        if you don't like to get the related warnings.

    iterative : bool, optional
        Defines whether iterative parsing should be used. It helps reduce
        memory usage at almost the same parsing speed. Default is
        :py:const:`True`.

    use_index : bool, optional
        Defines whether an index of byte offsets needs to be created for
        spectrum elements. Default is :py:const:`False`.

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

    return MzXML(source, read_schema=read_schema, iterative=iterative,
        use_index=use_index, dtype=dtype, huge_tree=huge_tree)


def iterfind(source, path, **kwargs):
    """Parse `source` and yield info on elements with specified local
    name or by specified XPath.

    .. note:: This function is provided for backward compatibility only.
        If you do multiple :py:func:`iterfind` calls on one file, you should
        create an :py:class:`MzXML` object and use its
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
        parameters. Disable this to avoid waiting on slow network connections or
        if you don't like to get the related warnings.

    decode_binary : bool, optional
        Defines whether binary data should be decoded and included in the output
        (under "m/z array", "intensity array", etc.).
        Default is :py:const:`True`.

    Returns
    -------
    out : iterator
    """
    return MzXML(source, **kwargs).iterfind(path, **kwargs)

version_info = xml._make_version_info(MzXML)

chain = aux._make_chain(read, 'read')
