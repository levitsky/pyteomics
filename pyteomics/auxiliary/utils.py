from __future__ import print_function

import base64
import zlib

from functools import wraps
from collections import namedtuple


try:
    basestring
except NameError:
    basestring = (str, bytes)

try:
    import numpy as np
except ImportError:
    np = None


def print_tree(d, indent_str=' -> ', indent_count=1):
    """Read a nested dict (with strings as keys) and print its structure.
    """
    def structure(d):
        out = {}
        for k, v in d.items():
            if isinstance(v, dict):
                out[k] = structure(v)
            elif isinstance(v, list) and v and isinstance(v[0], dict):
                out['{} [list]'.format(k)] = structure(v[0])
            else:
                out[k] = None
        return out

    def _print(d, level=0):
        for k, v in d.items():
            print('{}{}'.format(indent_str * indent_count * level, k))
            if v is not None:
                _print(v, level + 1)
    _print(structure(d))


def memoize(maxsize=1000):
    """Make a memoization decorator. A negative value of `maxsize` means
    no size limit."""
    def deco(f):
        """Memoization decorator. Items of `kwargs` must be hashable."""
        memo = {}

        @wraps(f)
        def func(*args, **kwargs):
            key = (args, frozenset(kwargs.items()))
            if key not in memo:
                if len(memo) == maxsize:
                    memo.popitem()
                memo[key] = f(*args, **kwargs)
            return memo[key]
        return func
    return deco


def _decode_base64_data_array(source, dtype, is_compressed):
    """Read a base64-encoded binary array.

    Parameters
    ----------
    source : str
        A binary array encoded with base64.
    dtype : dtype
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


class BinaryDataArrayTransformer(object):
    """A base class that provides methods for reading
    base64-encoded binary arrays.

    Attributes
    ----------
    compression_type_map : dict
        Maps compressor type name to decompression function
    """

    compression_type_map = {
        'no compression': lambda x: x,
        'zlib compression': zlib.decompress,
    }

    class binary_array_record(namedtuple(
            "binary_array_record", ("data", "compression", "dtype", "source", "key"))):
        """Hold all of the information about a base64 encoded array needed to
        decode the array.
        """

        def decode(self):
            """Decode :attr:`data` into a numerical array

            Returns
            -------
            np.ndarray
            """
            return self.source._decode_record(self)

    def _make_record(self, data, compression, dtype, key=None):
        return self.binary_array_record(data, compression, dtype, self, key)

    def _decode_record(self, record):
        array = self.decode_data_array(
            record.data, record.compression, record.dtype)
        return self._finalize_record_conversion(array, record)

    def _finalize_record_conversion(self, array, record):
        return array

    def _base64_decode(self, source):
        decoded_source = base64.b64decode(source.encode('ascii'))
        return decoded_source

    def _decompress(self, source, compression_type=None):
        if compression_type is None:
            return source
        decompressor = self.compression_type_map.get(compression_type)
        decompressed_source = decompressor(source)
        return decompressed_source

    def _transform_buffer(self, binary, dtype):
        output = np.frombuffer(binary, dtype=dtype)
        return output

    def decode_data_array(self, source, compression_type=None, dtype=np.float64):
        """Decode a base64-encoded, compressed bytestring into a numerical
        array.

        Parameters
        ----------
        source : bytes
            A base64 string encoding a potentially compressed numerical
            array.
        compression_type : str, optional
            The name of the compression method used before encoding the
            array into base64.
        dtype : type, optional
            The data type to use to decode the binary array from the
            decompressed bytes.

        Returns
        -------
        np.ndarray
        """
        binary = self._base64_decode(source)
        binary = self._decompress(binary, compression_type)
        array = self._transform_buffer(binary, dtype)
        return array
