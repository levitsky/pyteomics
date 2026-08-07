import base64
import zlib
import struct
from functools import wraps
from collections import namedtuple
from urllib.parse import urlparse
import pathlib

try:
    import numpy as np
except ImportError:
    np = None

try:
    import pynumpress
except ImportError:
    pynumpress = None


zstd = None
try:
    import pyzstd as zstd
except ImportError:
    try:
        from compression import zstd
    except ImportError:
        pass


from .structures import PyteomicsError


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
    output = np.frombuffer(bytearray(decoded_source), dtype=dtype)
    return output


_default_compression_map = {
        'no compression': lambda x: x,
        'zlib compression': zlib.decompress,
}


def _pynumpressDecompress(decoder):
    def decode(data):
        return decoder(np.frombuffer(data, dtype=np.uint8))
    return decode


def _zlibNumpress(decoder):
    def decode(data):
        return decoder(np.frombuffer(zlib.decompress(data), dtype=np.uint8))
    return decode


def _zstdNumpress(decoder):
    def decode(data):
        return decoder(np.frombuffer(zstd.decompress(data), dtype=np.uint8))

    return decode




def byte_shuffle_decode(byte_buffer, dtype):
    """
    Decode a byte-shuffled array into a specified data type.

    A byte-shuffled array of a multi-byte type stores the 1st byte of each value, then
    the second byte of each value, and so on in a sequence until all bytes are stored.
    The decoding process involves reconstructing the original value by putting the bytes
    from the ``i``th value back together from their disparate locations in the array.

    For sequential values, byteshuffling can make the values easier to compress because
    for small deltas most of the bytes in a value do not change.

    Parameters
    ----------
    byte_buffer : np.ndarray
        A numpy array of bytes whose size is evenly divisible by ``dtype.itemsize``
    dtype : np.dtype
        The real physical data type encoded in the array

    Returns
    -------
    np.ndarray :
        The array after shuffling
    """
    return byte_buffer.reshape((dtype().itemsize, -1)).T.flatten().view(dtype)


def dict_decode_with_shuffle(buffer, dtype) -> bytes:
    """
    Decode a dictionary encoded array with byte shuffling into a specified data type.

    A dictionary-encoded array efficiently stored repetitive values by storing unique values
    once in a dictionary, and then storing the positions of their repetitions as unique indices
    into that dictionary. Depending upon the number of unique values, these indices can be much
    smaller than the physical values stored in the dictionary. A byte-shuffling transform is applied
    to both the dictionary values and the indices to make them easier to compress as well.

    The dictionary is prefixed with an uint64 offset from the start of the byte buffer to the dictionary
    keys and a uint64 count of the number of values in the original array, so the overall size of the encoded
    array is ``8 + 8 + <size of unique values in bytes> + <size of keys in bytes>``.

    Parameters
    ----------
    byte_buffer : np.ndarray
        A numpy array of bytes which encode the dictionary and key array
    dtype : np.dtype
        The real physical data type encoded in the array

    Returns
    -------
    np.ndarray :
        The array after decoding
    """
    bytes_views = {
        8: np.uint64,
        4: np.uint32,
        2: np.uint16,
        1: np.uint8,
    }


    index_widths = [
        (2**8, np.uint8),
        (2**16, np.uint16),
        (2**32, np.uint32),
        (2**64, np.uint64),
    ]

    index_tp = None
    byte_view_tp = bytes_views[dtype.itemsize]

    (offset_to_data,) = struct.unpack("<Q", buffer[:8])
    (n_values,) = struct.unpack("<Q", buffer[8:16])
    for size, idx_type in index_widths:
        if size > n_values:
            index_tp = idx_type
            break
    else:
        raise ValueError(f"Cannot find index dtype for {n_values} indices")
    values = buffer[16:offset_to_data].view(byte_view_tp)
    indices = buffer[offset_to_data:].view(index_tp)

    indices = byte_shuffle_decode(indices.view(np.uint8), index_tp)
    values = byte_shuffle_decode(values.view(np.uint8), byte_view_tp)
    return values[indices].view(dtype)


DICTIONARY_ENCODED_ZSTD = "dictionary-encoded zstd compression"
BYTE_SHUFFLE_ZSTD = "byte-shuffled zstd compression"


if pynumpress:
    _default_compression_map.update(
        {
            'MS-Numpress short logged float compression': _pynumpressDecompress(pynumpress.decode_slof),
            'MS-Numpress positive integer compression':   _pynumpressDecompress(pynumpress.decode_pic),
            'MS-Numpress linear prediction compression':  _pynumpressDecompress(pynumpress.decode_linear),
            'MS-Numpress short logged float compression followed by zlib compression': _zlibNumpress(pynumpress.decode_slof),
            'MS-Numpress positive integer compression followed by zlib compression':   _zlibNumpress(pynumpress.decode_pic),
            'MS-Numpress linear prediction compression followed by zlib compression':  _zlibNumpress(pynumpress.decode_linear),
        })

if zstd:
    _default_compression_map.update(
        {
            "zstd compression": zstd.decompress,
            DICTIONARY_ENCODED_ZSTD: zstd.decompress,
            BYTE_SHUFFLE_ZSTD: zstd.decompress,
        }
    )

    if pynumpress:
        _default_compression_map.update(
            {
                "MS-Numpress short logged float compression followed by zstd compression": _zstdNumpress(
                    pynumpress.decode_slof
                ),
                "MS-Numpress positive integer compression followed by zstd compression": _zstdNumpress(
                    pynumpress.decode_pic
                ),
                "MS-Numpress linear prediction compression followed by zstd compression": _zlibNumpress(
                    pynumpress.decode_linear
                ),
            }
        )
else:
    def _zstd_missing(name: str):
        def _zstd_codec(*args, **kwargs):
            raise ModuleNotFoundError("pyzstd or the compression.zstd standard library module were absent, cannot decode {!r}".format(name))
        return _zstd_codec

    for c in [
        "zstd compression",
        DICTIONARY_ENCODED_ZSTD,
        BYTE_SHUFFLE_ZSTD,
        "MS-Numpress short logged float compression followed by zstd compression",
        "MS-Numpress positive integer compression followed by zstd compression",
        "MS-Numpress linear prediction compression followed by zstd compression",
    ]:

        _default_compression_map[c] = _zstd_missing(c)


class ArrayConversionMixin(object):
    _dtype_dict = {}
    _array_keys = ['m/z array', 'intensity array']

    def __init__(self, *args, **kwargs):
        self._dtype_dict = {None: None}
        dtype = kwargs.pop('dtype', None)
        if isinstance(dtype, dict):
            self._dtype_dict.update(dtype)
        elif dtype:
            self._dtype_dict = {k: dtype for k in self._array_keys}
            self._dtype_dict[None] = dtype
        self._convert_arrays = kwargs.pop('convert_arrays', 1)
        if self._convert_arrays and np is None:
            raise PyteomicsError('numpy is required for array conversion')
        super(ArrayConversionMixin, self).__init__(*args, **kwargs)

    def __getstate__(self):
        state = super(ArrayConversionMixin, self).__getstate__()
        state['_dtype_dict'] = self._dtype_dict
        state['_convert_arrays'] = self._convert_arrays
        state['_array_keys'] = self._array_keys
        return state

    def __setstate__(self, state):
        super(ArrayConversionMixin, self).__setstate__(state)
        self._dtype_dict = state['_dtype_dict']
        self._convert_arrays = state['_convert_arrays']
        self._array_keys = state['_array_keys']

    def _build_array(self, k, data):
        dtype = self._dtype_dict.get(k)
        return np.array(data, dtype=dtype)

    def _convert_array(self, k, array):
        dtype = self._dtype_dict.get(k)
        if dtype is not None:
            return array.astype(dtype)
        return array

    def _build_all_arrays(self, info):
        if self._convert_arrays:
            for k in self._array_keys:
                if k in info:
                    info[k] = self._build_array(k, info[k])


class MaskedArrayConversionMixin(ArrayConversionMixin):
    _masked_array_keys = ['charge array']
    _mask_value = 0

    def __init__(self, *args, **kwargs):
        self._convert_arrays = kwargs.pop('convert_arrays', 2)
        kwargs['convert_arrays'] = self._convert_arrays
        super(MaskedArrayConversionMixin, self).__init__(*args, **kwargs)

    def __getstate__(self):
        state = super(MaskedArrayConversionMixin, self).__getstate__()
        state['_masked_array_keys'] = self._masked_array_keys
        state['_mask_value'] = self._mask_value
        return state

    def __setstate__(self, state):
        super(MaskedArrayConversionMixin, self).__setstate__(state)
        self._masked_array_keys = state['_masked_array_keys']
        self._mask_value = state['_mask_value']

    def _build_masked_array(self, k, data):
        array = self._build_array(k, data)
        return self._convert_masked_array(k, array)

    def _convert_masked_array(self, k, array):
        return np.ma.masked_equal(array, self._mask_value)

    def _ensure_masked_array(self, k, data):
        if isinstance(data, np.ndarray):
            return self._convert_masked_array(k, data)
        return self._build_masked_array(self, k, data)

    def _build_all_arrays(self, info):
        super(MaskedArrayConversionMixin, self)._build_all_arrays(info)
        if self._convert_arrays == 2:
            for k in self._masked_array_keys:
                if k in info:
                    info[k] = self._ensure_masked_array(k, info[k])


if np is not None:
    class BinaryDataArrayTransformer(object):
        """A base class that provides methods for reading
        base64-encoded binary arrays.

        Attributes
        ----------
        compression_type_map : dict
            Maps compressor type name to decompression function
        """

        compression_type_map = _default_compression_map

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

        def _transform_buffer(self, binary, dtype, compression_type=None):
            if compression_type == BYTE_SHUFFLE_ZSTD:
                if not isinstance(binary, np.ndarray):
                    binary = np.frombuffer(binary, dtype=np.uint8)
                return byte_shuffle_decode(binary, dtype=dtype)
            elif compression_type == DICTIONARY_ENCODED_ZSTD:
                if not isinstance(binary, np.ndarray):
                    binary = np.frombuffer(binary, dtype=np.uint8)
                return dict_decode_with_shuffle(binary, dtype=dtype)
            if isinstance(binary, np.ndarray):
                return binary.astype(dtype, copy=False)
            return np.frombuffer(binary, dtype=dtype)

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
            if isinstance(binary, bytes):
                binary = bytearray(binary)
            array = self._transform_buffer(binary, dtype, compression_type)
            return array

    class BinaryArrayConversionMixin(ArrayConversionMixin, BinaryDataArrayTransformer):
        def _finalize_record_conversion(self, array, record):
            key = record.key
            return self._convert_array(key, array)


else:
    BinaryDataArrayTransformer = None
    BinaryArrayConversionMixin = None


def add_metaclass(metaclass):
    """Class decorator for creating a class with a metaclass."""
    def wrapper(cls):
        orig_vars = cls.__dict__.copy()
        slots = orig_vars.get('__slots__')
        if slots is not None:
            if isinstance(slots, str):
                slots = [slots]
            for slots_var in slots:
                orig_vars.pop(slots_var)
        orig_vars.pop('__dict__', None)
        orig_vars.pop('__weakref__', None)
        if hasattr(cls, '__qualname__'):
            orig_vars['__qualname__'] = cls.__qualname__
        return metaclass(cls.__name__, cls.__bases__, orig_vars)
    return wrapper


_supported_url_schemes = {'file', 'ftp', 'http', 'https'}


def ensure_url_prefix(path_or_url):
    """
    Ensure that a string path has a URL schema prefix, add "file://" if it is missing.
    Only expects a subset of URL schemes that can be urlopen-ed for reading.

    Parameters
    ----------
    path_or_url : str or Path

    Returns
    -------
    out : str
    """
    if isinstance(path_or_url, pathlib.Path):
        path = path_or_url
    else:
        parsed = urlparse(path_or_url)
        if parsed.scheme in _supported_url_schemes:
            return path_or_url
        if parsed.params or parsed.query or parsed.fragment:
            return path_or_url
        # for all other cases assume it's a local path
        path = pathlib.Path(path_or_url)
    return path.absolute().as_uri()
