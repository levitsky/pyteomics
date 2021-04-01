# -*- coding: utf8 -*-
"""
mzmlb - reader for mass spectrometry data in mzMLb format
=========================================================

.. warning::
    This is a **Provisional Implementation**. The mzMLb format has been published
    but is not yet broadly available.

Summary
-------
mzMLb is an HDF5 container format wrapping around the standard rich XML-format
for raw mass spectrometry data storage. Please refer to [1]_ for more information
about mzMLb and its features. Please refer to
`psidev.info <http://www.psidev.info/index.php?q=node/257>`_ for the detailed
specification of the format and structure of mzML files.

This module provides a minimalistic way to extract information from mzML
files. You can use the old functional interface (:py:func:`read`) or the new
object-oriented interface (:py:class:`MzMLb` to iterate over entries in ``<spectrum>`` elements.
:py:class:`MzMLb` also support direct indexing with spectrum IDs or indices.

Data access
-----------

  :py:class:`MzMLb` - a class representing a single mzMLb file.
  Other data access functions use this class internally.

  :py:func:`read` - iterate through spectra in mzMLb file. Data from a
  single spectrum are converted to a human-readable dict. Spectra themselves are
  stored under 'm/z array' and 'intensity array' keys.


References
----------
.. [1] Bhamber, R. S., Jankevics, A., Deutsch, E. W., Jones, A. R., & Dowsey, A. W. (2021).
    MzMLb: A Future-Proof Raw Mass Spectrometry Data Format Based on Standards-Compliant
    mzML and Optimized for Speed and Storage Requirements. Journal of Proteome Research,
    20(1), 172–183. https://doi.org/10.1021/acs.jproteome.0c00192
"""

import io
import warnings
import logging
from collections import namedtuple

import h5py
try:
    logging.getLogger("hdf5plugin").addHandler(logging.NullHandler())
    import hdf5plugin
except ImportError:
    hdf5plugin = None

import numpy as np

from pyteomics.mzml import MzML as _MzML
from pyteomics.auxiliary.file_helpers import HierarchicalOffsetIndex, TaskMappingMixin, TimeOrderedIndexedReaderMixin, FileReader


def delta_predict(data, copy=True):
    '''Reverse the lossy transformation of the delta compression
    helper.

    Parameters
    ----------
    data : :class:`numpy.ndarray`
        The data to transform
    copy : bool
        Whether to make a copy of the data array or transform it in-place.

    Returns
    -------
    :class:`numpy.ndarray`
        The transformed data array
    '''
    if copy:
        out = data.copy()
    else:
        out = data
    for i in range(2, len(data)):
        out[i] = out[i] + out[i - 1] - out[0]
    return out


def linear_predict(data, copy=True):
    '''Reverse the lossy transformation of the linear interpolation compression
    helper.

    Parameters
    ----------
    data : :class:`numpy.ndarray`
        The data to transform
    copy : bool
        Whether to make a copy of the data array or transform it in-place.

    Returns
    -------
    :class:`numpy.ndarray`
        The transformed data array
    '''
    if copy:
        out = data.copy()
    else:
        out = data
    for i in range(2, len(data)):
        out[i] = out[i] + 2 * out[i - 1] - out[i - 2] - out[1]
    return out


class HDF5ByteBuffer(io.RawIOBase):
    '''Helper class that looks file-like so that we can pass a HDF5 byte dataset to
    an arbitrary XML parser.

    Implements :class:`~io.RawIOBase` for reading.
    '''
    def __init__(self, buffer, offset=None):
        if offset is None:
            offset = 0
        self.buffer = buffer
        self.offset = offset
        self.size = self.buffer.size

    def seek(self, offset, whence=0):
        if whence == io.SEEK_SET:
            self.offset = offset
        elif whence == io.SEEK_CUR:
            self.offset += offset
        elif whence == io.SEEK_END:
            self.offset = self.size - offset
        else:
            raise ValueError("Bad whence %r" % whence)

    def tell(self):
        return self.offset

    def close(self):
        return

    def closed(self):
        return False

    def readinto(self, b):
        n = len(b)
        temp = self._read(n)
        m = len(temp)
        b[:m] = temp[:]
        return m

    def readall(self):
        return bytes(self._read(-1))

    def _read(self, n=-1):
        if n == -1:
            n = self.size + 1
        dat = bytearray(self.buffer[self.offset:self.offset + n])
        self.offset += n
        return dat


class external_array_slice(namedtuple('external_array_slice',
                           ['array_name', 'offset', 'length', 'source', 'transform', 'key', 'dtype'])):
    def decode(self):
        """Decode :attr:`data` into a numerical array

        Returns
        -------
        np.ndarray
        """
        return self.source._decode_record(self)


class ExternalDataMzML(_MzML):
    '''An MzML parser that reads data arrays from an external provider.

    This is an implementation detail of :class:`MzMLb`.
    '''
    def __init__(self, *args, **kwargs):
        self._external_data_registry = kwargs.pop("external_data_registry", None)
        super(ExternalDataMzML, self).__init__(*args, **kwargs)

    def _make_record(self, array_name, offset, length, transform, name, dtype):
        return external_array_slice(array_name, offset, length, self, transform, name, dtype)

    def _transform_array(self, array, transform):
        if transform is None:
            return array
        elif "linear prediction" == transform:
            return linear_predict(array, copy=False)
        elif "delta prediction" == transform:
            return delta_predict(array, copy=False)
        else:
            raise ValueError("Transformation not recognized")

    def _retrieve_external_array(self, array_name, length, offset):
        array = self._external_data_registry.get(array_name, length, offset)
        return array

    def decode_data_array(self, array_name, offset, length, transform=None, dtype=np.float64):
        array = self._retrieve_external_array(array_name, length, offset)
        array = self._transform_array(array, transform)
        return array

    def _decode_record(self, record):
        array = self.decode_data_array(
            record.array_name, record.offset, record.length, record.transform, record.dtype)
        return self._finalize_record_conversion(array, record)

    def _handle_binary(self, info, **kwargs):
        if not self.decode_binary:
            self.decode_binary = True
            # Binary decoding works totally differently here, not supporting the previous signatures
            # that the parent method will use. Pretend we are decoding because it is a no-op in the
            # parent method.
            result = super(ExternalDataMzML, self)._handle_binary(info, **kwargs)
            self.decode_binary = False
        else:
            result = super(ExternalDataMzML, self)._handle_binary(info, **kwargs)
        try:
            array_name = info['external HDF5 dataset']
        except KeyError:
            array_name = info['external dataset']
        offset = int(info['external offset'])
        length = int(info['external array length'])

        transform = None
        # The zlib compression in these two terms happens automatically during HDF5 encoding and
        # the reader needn't even know about it. Need an example of how Numpress will be signaled.
        if "linear prediction" in info or "truncation, linear prediction and zlib compression" in info:
            transform = 'linear prediction'
        elif "delta prediction" in info or "truncation, delta prediction and zlib compression" in info:
            transform = 'delta prediction'

        if not self.decode_binary:
            name = self._detect_array_name(info)
            result[name] = self._make_record(
                array_name, offset, length, transform, name,
                self._external_data_registry.dtype_of(array_name))
            return result

        array = self._retrieve_external_array(array_name, length, offset)

        if len(result) == 1:
            name = next(iter(result))
        else:
            name = self._detect_array_name(info)
        result[name] = self._convert_array(name, array)
        return result


class ExternalArrayRegistry(object):
    '''Read chunks out of a single long array

    This is an implementation detail of :class:`MzMLb`
    '''
    def __init__(self, registry):
        self.registry = registry

    def get(self, array_name, length, offset=0):
        return self.registry[array_name][offset:offset + length]

    def dtype_of(self, array_name):
        return self.registry[array_name].dtype

    def __call__(self, array_name, length, offset=0):
        return self.get(array_name, length, offset)


class MzMLb(TimeOrderedIndexedReaderMixin, TaskMappingMixin):
    '''A parser for mzMLb [1]_

    Provides an identical interface to :class:`~pyteomics.mzml.MzML`

    Attributes
    ----------
    path : str, Path-like, or file-like object
        The mzMLb file path or a file-like object providing it.
    handle : :class:`h5py.File`
        The raw HDF5 file container.
    mzml_parser : :class:`~.ExternalDataMzML`
        The mzML parser for the XML stream inside the HDF5 file with
        special behavior for retrieving the out-of-band data arrays
        from their respective storage locations.
    schema_version : str
        The mzMLb HDF5 schema version, distinct from the mzML schema inside it


    References
    ----------
    [1] Bhamber, R. S., Jankevics, A., Deutsch, E. W., Jones, A. R., & Dowsey, A. W. (2021).
        MzMLb: A Future-Proof Raw Mass Spectrometry Data Format Based on Standards-Compliant
        mzML and Optimized for Speed and Storage Requirements. Journal of Proteome Research,
        20(1), 172–183. https://doi.org/10.1021/acs.jproteome.0c00192
    '''
    _default_iter_tag = ExternalDataMzML._default_iter_tag

    file_format = "mzMLb"

    def __init__(self, path, hdfargs=None, mzmlargs=None, allow_updates=False,
                 use_index=True, **kwargs):
        if hdfargs is None:
            hdfargs = {}
        if mzmlargs is None:
            mzmlargs = {}
        mzmlargs.update(kwargs)

        self.path = path
        self._hdfargs = hdfargs
        self._mzmlargs = mzmlargs
        self._allow_updates = allow_updates
        self.handle = h5py.File(self.path, 'r+' if self._allow_updates else 'r', **hdfargs)
        self.schema_version = self.handle['mzML'].attrs.get('version')
        self._check_compressor()

        self._xml_buffer = HDF5ByteBuffer(self.handle['mzML'])
        self._array_registry = ExternalArrayRegistry(self.handle)
        self._make_mzml_parser(mzmlargs)

        super(MzMLb, self).__init__(**kwargs)

    def _check_compressor(self):
        for key in self.handle.keys():
            if "spectrum_MS_" in key or "chromatogram_MS_":
                data = self.handle[key]
                try:
                    filts = data._filters
                except AttributeError:
                    continue
                if '32001' in filts:
                    if hdf5plugin is None:
                        warnings.warn(
                            ("Blosc meta-compressor detected, but hdf5plugin is "
                             "not installed, may not be able to access %r") % (key))

    def _make_mzml_parser(self, kwargs):
        self._mzml_parser = ExternalDataMzML(
            self._xml_buffer, external_data_registry=self._array_registry,
            use_index=False, **kwargs)
        self._mzml_parser._offset_index = self._build_index()
        self._mzml_parser._use_index = True

    @property
    def name(self):
        if hasattr(self.path, 'name'):
            return self.path.name
        return self.path

    def _build_index(self):
        index = HierarchicalOffsetIndex()
        for label in [u'spectrum', u'chromatogram']:
            sub = index[label]
            ids = bytearray(self.handle['mzML_{}Index_idRef'.format(label)]).split(b"\x00")
            offsets = self.handle["mzML_{}Index".format(label)][:-1]
            for i, o in enumerate(offsets):
                sub[ids[i].decode('utf8')] = o
        return index

    def get_by_id(self, id):
        """Parse the file and return the element with `id` attribute equal
        to `elem_id`. Returns :py:const:`None` if no such element is found.

        Parameters
        ----------
        elem_id : str
            The value of the `id` attribute to match.

        Returns
        -------
        out : :py:class:`dict` or :py:const:`None`
        """
        return self._mzml_parser.get_by_id(id)

    def get_by_ids(self, ids):
        return self._mzml_parser.get_by_ids(ids)

    def get_by_index(self, i):
        return self._mzml_parser.get_by_index(i)

    def get_by_indexes(self, indexes):
        return self._mzml_parser.get_by_indexes(indexes)

    def get_by_index_slice(self, s):
        return self._mzml_parser.get_by_index_slice(s)

    def get_by_key_slice(self, s):
        return self._mzml_parser.get_by_key_slice(s)

    def __contains__(self, key):
        return key in self.index

    def __getitem__(self, i):
        return self._mzml_parser[i]

    def __len__(self):
        return len(self._mzml_parser)

    def __iter__(self):
        return iter(self._mzml_parser)

    def __next__(self):
        return next(self._mzml_parser)

    def next(self):
        return self.__next__()

    def __reduce__(self):
        return self.__class__, (self.path, self._hdfargs, self._mzmlargs, self._allow_updates)

    def close(self):
        self.handle.close()

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def iterfind(self, *args, **kwargs):
        iterf = self._mzml_parser.iterfind(*args, **kwargs)
        iterf.parser = self
        return iterf

    def _iterfind_impl(self, path, *args, **kwargs):
        return self._mzml_parser._iterfind_impl(path, *args, **kwargs)

    @property
    def index(self):
        return self._mzml_parser.index

    @property
    def _offset_index(self):
        return self._mzml_parser._offset_index

    @property
    def default_index(self):
        return self._mzml_parser.default_index

    def _get_time(self, scan):
        return self._mzml_parser._get_time(scan)

    @property
    def mzml_parser(self):
        return self._mzml_parser

    def _task_map_iterator(self):
        """Returns the :class:`Iteratable` to use when dealing work items onto the input IPC
        queue used by :meth:`map`

        Returns
        -------
        :class:`Iteratable`
        """
        return iter(self.index[self._default_iter_tag])

    def read(self, n=-1):
        return self._mzml_parser.read(n)

    def reset(self):
        self._mzml_parser.reset()

    def seek(self, offset, whence=0):
        self._mzml_parser.seek(offset, whence)

    def tell(self):
        return self._mzml_parser.tell()

    def get_dataset(self, name):
        '''Get an HDF5 dataset by its name or path relative to
        the root node.

        .. warning::
            Because this accesses HDF5 data directly, it may be possible to mutate
            the underlying file if :attr:`allow_updates` is :const:`True`.

        Parameters
        ----------
        name : :class:`str`
            The dataset name or path.

        Returns
        -------
        :class:`h5py.Dataset` or :class:`h5py.Group`

        Raises
        ------
        KeyError :
            The name is not found.
        '''
        return self.handle[name]


def read(source, dtype=None):
    """Parse `source` and iterate through spectra.

    Parameters
    ----------
    source : str or file
        A path to a target mzML file or the file object itself.
    dtype : type or dict, optional
        dtype to convert arrays to, one for both m/z and intensity arrays or one for each key.
        If :py:class:`dict`, keys should be 'm/z array' and 'intensity array'.

    Returns
    -------
    out : iterator
       An iterator over the dicts with spectrum properties.
    """
    reader = MzMLb(source, dtype=dtype)
    return reader


# The MzMLb class is detatched from the normal :class:`FileReader`-based inheritance tree,
# this grafts it back on for :func:`isinstance` and :func:`issubclass` tests at least.
FileReader.register(MzMLb)