"""
ms2 - read and write MS/MS data in MS2 format
=============================================

Summary
-------

`MS2 <http://dx.doi.org/10.1002/rcm.1603>`_ is a simple
human-readable format for MS2 data. It allows storing MS2 peak lists and
exprimental parameters.

This module provides minimalistic infrastructure for access to data stored in
MS2 files.
Two main classes are :py:class:`MS2`, which provides an iterative, text-mode parser,
and :py:class:`IndexedMS2`, which is a binary-mode parser that supports random access using scan IDs
and retention times.
The function :py:func:`read` helps dispatch between the two classes.
Also, common parameters can be read from MS2 file header with
:py:func:`read_header` function.

Classes
-------

  :py:class:`MS2` - a text-mode MS2 parser. Suitable to read spectra from a file consecutively.
  Needs a file opened in text mode (or will open it if given a file name).

  :py:class:`IndexedMS2` - a binary-mode MS2 parser. When created, builds a byte offset index
  for fast random access by spectrum ID. Sequential iteration is also supported.
  Needs a seekable file opened in binary mode (if created from existing file object).

  :py:class:`MS2Base` - abstract class, the common ancestor of the two classes above.
  Can be used for type checking.

Functions
---------

  :py:func:`read` - an alias for :py:class:`MS2` or :py:class:`IndexedMS1`.

  :py:func:`chain` - read multiple files at once.

  :py:func:`chain.from_iterable` - read multiple files at once, using an
  iterable of files.

  :py:func:`read_header` - get a dict with common parameters for all spectra
  from the beginning of MS2 file.

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

from pyteomics import auxiliary as aux
from pyteomics.ms1 import MS1, IndexedMS1, MS1Base


class MS2Base(aux.MaskedArrayConversionMixin, MS1Base):
    """Abstract class representing an MS2 file. Subclasses implement different approaches to parsing."""
    _array_keys = ['m/z array', 'intensity array', 'charge array', 'resolution array']
    _float_keys = ['RTime', 'RetTime', 'IonInjectionTime', 'PrecursorInt']

    def __init__(self, source=None, use_header=False, convert_arrays=2, dtype=None, read_charges=True, read_resolutions=True, encoding=None, **kwargs):
        """
        Create an instance of a :py:class:`MS2Base` parser.

        Parameters
        ----------

        source : str or file or None, optional
            A file object (or file name) with data in MS1 format. Default is
            :py:const:`None`, which means read standard input.

        use_header : bool, optional
            Add the info from file header to each dict. Spectrum-specific parameters
            override those from the header in case of conflict.
            Default is :py:const:`False`.

        convert_arrays : one of {0, 1, 2}, optional
            If `0`, m/z, intensities and (possibly) charges will be returned as regular lists.
            If `1`, they will be converted to regular :py:class:`numpy.ndarray`'s.
            If `2`, charges will be reported as a masked array (default).
            The default option is the slowest. `1` and `2` require :py:mod:`numpy`.

        read_charges : bool, optional
            If `True` (default), fragment charges are reported. Disabling it improves performance.
            Charge is expected to be the **third** number on the line, after peak *m/z* and intensity.

        read_resolutions : bool, optional
            If `True` (default), fragment peak resolutions are reported. Disabling it improves performance.
            Resolution is expected to be the **fourth** number on the line, after peak *m/z*, intensity, and charge.

        dtype : type or str or dict, optional
            dtype argument to :py:mod:`numpy` array constructor, one for all arrays or one for each key.
            Keys should be 'm/z array', 'intensity array', 'charge array', 'resolution array'.

        encoding : str, optional
            File encoding.
        """
        super(MS2Base, self).__init__(source=source, use_header=use_header, convert_arrays=convert_arrays, dtype=dtype,
            encoding=encoding, **kwargs)
        self._read_charges = read_charges
        self._read_resolutions = read_resolutions

    def _handle_peak(self, line, sline, info):
        super(MS2Base, self)._handle_peak(line, sline, info)
        if self._read_charges:
            if len(sline) > 2:
                sline = line.strip().split()
                try:
                    info['charge array'].append(int(sline[2]))
                except ValueError:
                    raise aux.PyteomicsError("Error parsing fragment charge on line: " + line)
            else:
                info['charge array'].append(0)
        if self._read_resolutions:
            if len(sline) > 2:
                sline = line.strip().split()
                try:
                    info['resolution array'].append(int(sline[3]))
                except ValueError:
                    raise aux.PyteomicsError("Error parsing fragment peak resolution on line: " + line)
            else:
                info['resolution array'].append(0)

    def _make_scan(self, info):
        if not self._read_charges:
            del info['charge array']
        if not self._read_resolutions:
            del info['resolution array']
        return super(MS2Base, self)._make_scan(info)

    def __reduce_ex__(self, protocol):
        return (self.__class__,
            (self._source_init, False, self._convert_arrays, None, self._read_charges, self._read_resolutions, self.encoding),
            self.__getstate__())


class MS2(MS2Base, MS1):
    """
    A class representing an MS2 file. Supports the `with` syntax and direct iteration for sequential
    parsing.

    :py:class:`MS2` object behaves as an iterator, **yielding** spectra one by one.
    Each 'spectrum' is a :py:class:`dict` with three keys: 'm/z array',
    'intensity array', and 'params'. 'm/z array' and
    'intensity array' store :py:class:`numpy.ndarray`'s of floats,
    and 'params' stores a :py:class:`dict` of parameters.

    Attributes
    ----------

    header : dict
        The file header.

    """
    def __init__(self, *args, **kwargs):
        """
        Create an :py:class:`MS2` (text-mode) reader for a given MS2 file.

        Parameters
        ----------

        source : str or file or None, optional
            A file object (or file name) with data in MS2 format. Default is
            :py:const:`None`, which means read standard input.

            .. note :: If a file object is given, it must be opened in text mode.

        use_header : bool, optional
            Add the info from file header to each dict. Spectrum-specific parameters
            override those from the header in case of conflict.
            Default is :py:const:`False`.

        convert_arrays : one of {0, 1, 2}, optional
            If `0`, m/z, intensities and (possibly) charges will be returned as regular lists.
            If `1`, they will be converted to regular :py:class:`numpy.ndarray`'s.
            If `2`, charges will be reported as a masked array (default).
            The default option is the slowest. `1` and `2` require :py:mod:`numpy`.

        read_charges : bool, optional
            If `True` (default), fragment charges are reported. Disabling it improves performance.
            Charge is expected to be the **third** number on the line, after peak *m/z* and intensity.

        read_resolutions : bool, optional
            If `True` (default), fragment peak resolutions are reported. Disabling it improves performance.
            Resolution is expected to be the **fourth** number on the line, after peak *m/z*, intensity, and charge.

        dtype : type or str or dict, optional
            dtype argument to :py:mod:`numpy` array constructor, one for all arrays or one for each key.
            Keys should be 'm/z array', 'intensity array', 'charge array'.

        encoding : str, optional
            File encoding.

        Returns
        -------

        out : MS2
            The reader object.
        """
        super(MS2, self).__init__(*args, **kwargs)


class IndexedMS2(IndexedMS1, MS2Base):
    """
    A class representing an MS2 file. Supports the `with` syntax and direct iteration for sequential
    parsing. Specific spectra can be accessed by title using the indexing syntax in constant time.
    If created using a file object, it needs to be opened in binary mode.

    When iterated, :py:class:`IndexedMS2` object yields spectra one by one.
    Each 'spectrum' is a :py:class:`dict` with four keys: 'm/z array',
    'intensity array', 'charge array' and 'params'. 'm/z array' and
    'intensity array' store :py:class:`numpy.ndarray`'s of floats,
    'charge array' is a masked array (:py:class:`numpy.ma.MaskedArray`) of ints,
    and 'params' stores a :py:class:`dict` of parameters (keys and values are
    :py:class:`str`, keys corresponding to MS2).

    .. warning ::
        Labels for scan objects are constructed as the first number in the S line, as follows:
        for a line ``S  0   1   123.4`` the label is `'0'`. If these labels are not unique
        for the scans in the file, the indexed parser will not work correctly. Consider using
        :py:class:`MS2` instead.

    Attributes
    ----------

    header : dict
        The file header.
    time : RTLocator
        A property used for accessing spectra by retention time.
    """
    def __init__(self, source=None, use_header=False, convert_arrays=2, dtype=None, read_charges=True, read_resolutions=True,
                 encoding='utf-8', _skip_index=False, **kwargs):
        """
        Create an :py:class:`IndexedMS2` (binary-mode) reader for a given MS2 file.

        Parameters
        ----------

        source : str or file or None, optional
            A file object (or file name) with data in MS2 format. Default is
            :py:const:`None`, which means read standard input.

            .. note :: If a file object is given, it must be opened in binary mode.

        use_header : bool, optional
            Add the info from file header to each dict. Spectrum-specific parameters
            override those from the header in case of conflict.
            Default is :py:const:`True`.

        convert_arrays : one of {0, 1, 2}, optional
            If `0`, m/z, intensities and (possibly) charges will be returned as regular lists.
            If `1`, they will be converted to regular :py:class:`numpy.ndarray`'s.
            If `2`, charges will be reported as a masked array (default).
            The default option is the slowest. `1` and `2` require :py:mod:`numpy`.

        read_charges : bool, optional
            If `True` (default), fragment charges are reported. Disabling it improves performance.
            Charge is expected to be the **third** number on the line, after peak *m/z* and intensity.

        read_resolutions : bool, optional
            If `True` (default), fragment peak resolutions are reported. Disabling it improves performance.
            Resolution is expected to be the **fourth** number on the line, after peak *m/z*, intensity, and charge.

        dtype : type or str or dict, optional
            dtype argument to :py:mod:`numpy` array constructor, one for all arrays or one for each key.
            Keys should be 'm/z array', 'intensity array', 'charge array'.

        encoding : str, optional
            File encoding.

        block_size : int, optinal
            Size of the chunk (in bytes) used to parse the file when creating the byte offset index.

        Returns
        -------

        out : IndexedMS2
            The reader object.
        """
        super(IndexedMS2, self).__init__(source, use_header=use_header, convert_arrays=convert_arrays, dtype=dtype,
            read_charges=read_charges, read_resolutions=read_resolutions, encoding=encoding, _skip_index=_skip_index, **kwargs)

    def __reduce_ex__(self, protocol):
        return (self.__class__,
            (self._source_init, False, self._convert_arrays, None, self._read_charges, self._read_resolutions, self.encoding, True),
            self.__getstate__())


def read_header(source, *args, **kwargs):
    """
    Read the specified MS2 file, get the parameters specified in the header
    as a :py:class:`dict`.

    Parameters
    ----------

    source : str or file
        File name or file object representing an file in MS2 format.

    Returns
    -------

    header : dict
    """
    kwargs['use_header'] = True
    return read(source, *args, **kwargs).header


def read(*args, **kwargs):
    """Read an MS2 file and return entries iteratively.

    Read the specified MS2 file, **yield** spectra one by one.
    Each 'spectrum' is a :py:class:`dict` with three keys: 'm/z array',
    'intensity array', and 'params'. 'm/z array' and
    'intensity array' store :py:class:`numpy.ndarray`'s of floats,
    and 'params' stores a :py:class:`dict` of parameters.

    Parameters
    ----------

    source : str or file or None, optional
        A file object (or file name) with data in MS2 format. Default is
        :py:const:`None`, which means read standard input.

    use_header : bool, optional
        Add the info from file header to each dict. Spectrum-specific parameters
        override those from the header in case of conflict.
        Default is :py:const:`False`.

    convert_arrays : bool, optional
        If :py:const:`False`, m/z and intensities will be returned as regular lists.
        If :py:const:`True` (default), they will be converted to regular :py:class:`numpy.ndarray`'s.
        Conversion requires :py:mod:`numpy`.

    read_charges : bool, optional
        If `True` (default), fragment charges are reported. Disabling it improves performance.
        Charge is expected to be the **third** number on the line, after peak *m/z* and intensity.

    read_resolutions : bool, optional
        If `True` (default), fragment peak resolutions are reported. Disabling it improves performance.
        Resolution is expected to be the **fourth** number on the line, after peak *m/z*, intensity, and charge.

    dtype : type or str or dict, optional
        dtype argument to :py:mod:`numpy` array constructor, one for all arrays or one for each key.
        Keys should be 'm/z array' and/or 'intensity array'.

    encoding : str, optional
        File encoding.

    use_index : bool, optional
        Determines which parsing method to use. If :py:const:`True`, an instance of
        :py:class:`IndexedMS2` is created. This facilitates random access by scan titles.
        If an open file is passed as `source`, it needs to be open in binary mode.

        .. warning ::
            Labels for scan objects are constructed as the first number in the S line, as follows:
            for a line ``S  0   1   123.4`` the label is `'0'`. If these labels are not unique
            for the scans in the file, the indexed parser will not work correctly.

        If :py:const:`False` (default), an instance of :py:class:`MS2` is created. It reads
        `source` in text mode and is suitable for iterative parsing.

    block_size : int, optinal
        Size of the chunk (in bytes) used to parse the file when creating the byte offset index.
        (Accepted only for :py:class:`IndexedMS2`.)

    Returns
    -------

    out :
        An instance of :py:class:`MS2` or :py:class:`IndexedMS2`, depending on `use_index` and `source`.
    """
    if args:
        source = args[0]
    else:
        source = kwargs.get('source')
    use_index = kwargs.pop('use_index', None)
    use_index = aux._check_use_index(source, use_index, False)
    tp = IndexedMS2 if use_index else MS2

    return tp(*args, **kwargs)


chain = aux._make_chain(read, 'read')
