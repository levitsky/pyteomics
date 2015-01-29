"""
mzml - reader for mass spectrometry data in mzML format
=======================================================

Summary
-------

mzML is a standard rich XML-format for raw mass spectrometry data storage.
Please refer to http://www.psidev.info/index.php?q=node/257 for the detailed
specification of the format and the structure of mzML files.

This module provides a minimalistic way to extract information from mzIdentML
files. You can use the old functional interface (:py:func:`read`) or the new
object-oriented interface (:py:class:`MzML`) to iterate over entries in
``<spectrum>`` elements.

Data access
-----------

  :py:class:`MzML` - a class representing a single mzML file.
  Other data access functions use this class internally.

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
from . import xml, auxiliary as aux

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

class MzML(xml.XML):
    """Parser class for mzML files."""
    file_format = 'mzML'
    _root_element = 'mzML'
    _default_schema = xml._mzml_schema_defaults
    _default_version = '1.1.0'
    _default_iter_tag = 'spectrum'
    _structures_to_flatten = {'binaryDataArrayList'}

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
                info['binary'] == array
        if 'binaryDataArray' in info:
            for array in info.pop('binaryDataArray'):
                info.update(array)
        intkeys = {'ms level'}
        for k in intkeys:
            if k in info:
                info[k] = int(info[k])
        return info

def read(source, read_schema=True, iterative=True):
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

    Returns
    -------
    out : iterator
       An iterator over the dicts with spectra properties.
    """

    return MzML(source, read_schema=read_schema, iterative=iterative)

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
