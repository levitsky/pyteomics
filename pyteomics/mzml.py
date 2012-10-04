"""
mzml - reader for mass spectrometry data in mzML format
=======================================================

Summary
-------

mzML is a standard rich XML-format for raw mass spectrometry data storage.
Please refer to http://www.psidev.info/index.php?q=node/257 for the detailed 
specification of the format and the structure of mzML files.

This module provides minimalistic infrastructure for access to data stored in
mzML files. The most important function is :py:func:`read`, which 
reads spectra and related information as saves them into human-readable dicts.
The rest of data can be obtained via a combination of :py:func:`get_node` and
:py:func:`read_params` functions. These functions rely on the terminology of 
the underlying `lxml library <http://lxml.de/>`_. 

Data access
-----------

  :py:func:`read` - iterate through spectra in mzML file. Data from a
  single spectrum are converted to a human-readable dict. Spectra themselves are 
  stored under 'm/z array' and 'intensity array' keys.

  :py:func:`version_info` - get version information about the mzML file

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

import numpy
import zlib
import base64
from lxml import etree
from . import auxiliary as aux

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
    output = numpy.frombuffer(decoded_source, dtype=dtype)
    return output

def read(source):
    """Parse ``source`` and iterate through spectra.

    Parameters
    ----------
    source : str or file
        A path to a target mzML file or the file object itself.

    Returns
    -------
    out : iterator
       An iterator over the dicts with spectra properties.
    """
    
    return _itertag(source, 'spectrum')

def _get_info_smart(source, element, **kw):
    name = aux._local_name(element)
    kwargs = dict(kw)
    rec = kwargs.pop('recursive', None)
    if name in {'indexedmzML', 'mzML'}:
        info =  _get_info(source, element, rec if rec is not None else False,
                **kwargs)
    else:
        info = _get_info(source, element, rec if rec is not None else True,
                **kwargs)
    if 'binary' in info:
        types = {'32-bit float': 'f', '64-bit float': 'd'}
        for t, code in types.items():
            if t in info:
                dtype = code
                del info[t]
        if 'zlib compression' in info:
            compressed = True
            del info['zlib compression']
        else:
            compressed = False
            info.pop('no compression', None)
        b = info.pop('binary')
        if b:
            array = _decode_base64_data_array(
                            b, dtype, compressed)
        else:
            array = numpy.array([], dtype=dtype)
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

_version_info_env = {'format': 'mzML', 'element': 'mzML'}
version_info = aux._make_version_info(_version_info_env)

_schema_defaults = {'ints': {
    ('spectrum', 'index'),
     ('instrumentConfigurationList', 'count'),
     ('binaryDataArray', 'encodedLength'),
     ('cvList', 'count'),
     ('binaryDataArray', 'arrayLength'),
     ('scanWindowList', 'count'),
     ('componentList', 'count'),
     ('sourceFileList', 'count'),
     ('productList', 'count'),
     ('referenceableParamGroupList', 'count'),
     ('scanList', 'count'),
     ('spectrum', 'defaultArrayLength'),
     ('dataProcessingList', 'count'),
     ('sourceFileRefList', 'count'),
     ('scanSettingsList', 'count'),
     ('selectedIonList', 'count'),
     ('chromatogram', 'defaultArrayLength'),
     ('precursorList', 'count'),
     ('chromatogram', 'index'),
     ('processingMethod', 'order'),
     ('targetList', 'count'),
     ('sampleList', 'count'),
     ('softwareList', 'count'),
     ('binaryDataArrayList', 'count'),
     ('spectrumList', 'count'),
     ('chromatogramList', 'count')},
        'floats': {},
        'bools': {},
        'lists': {'scan', 'spectrum', 'sample', 'cv', 'dataProcessing',
            'cvParam', 'source', 'userParam', 'detector', 'product',
            'referenceableParamGroupRef', 'selectedIon', 'sourceFileRef',
            'binaryDataArray', 'analyzer', 'scanSettings',
            'instrumentConfiguration', 'chromatogram', 'target',
            'processingMethod', 'precursor', 'sourceFile',
            'referenceableParamGroup', 'contact', 'scanWindow', 'software'},
        'intlists': {},
        'floatlists': {},
        'charlists': {}}
_schema_env = {'format': 'mzML', 'version_info': version_info,
        'default_version': '1.1.0', 'defaults': _schema_defaults}
_schema_info = aux._make_schema_info(_schema_env)

_getinfo_env = {'keys': {'binaryDataArrayList'}, 'schema_info': _schema_info,
        'get_info_smart': _get_info_smart}
_get_info = aux._make_get_info(_getinfo_env)

_itertag_env = {'get_info_smart': _get_info_smart}
_itertag = aux._make_itertag(_itertag_env)
