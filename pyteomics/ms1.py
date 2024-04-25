"""
ms1 - read and write MS/MS data in MS1 format
=============================================

Summary
-------

`MS1 <http://dx.doi.org/10.1002/rcm.1603>`_ is a simple
human-readable format for MS1 data. It allows storing MS1 peak lists and exprimental parameters.

This module provides minimalistic infrastructure for access to data stored in MS1 files.
Two main classes are :py:class:`MS1`, which provides an iterative, text-mode parser,
and :py:class:`IndexedMS1`, which is a binary-mode parser that supports random access using scan IDs
and retention times.
The function :py:func:`read` helps dispatch between the two classes.
Also, common parameters can be read from MS1 file header with :py:func:`read_header` function.

Classes
-------

  :py:class:`MS1` - a text-mode MS1 parser. Suitable to read spectra from a file consecutively.
  Needs a file opened in text mode (or will open it if given a file name).

  :py:class:`IndexedMS1` - a binary-mode MS1 parser. When created, builds a byte offset index
  for fast random access by spectrum ID. Sequential iteration is also supported.
  Needs a seekable file opened in binary mode (if created from existing file object).

  :py:class:`MS1Base` - abstract class, the common ancestor of the two classes above.
  Can be used for type checking.

Functions
---------

  :py:func:`read` - an alias for :py:class:`MS1` or :py:class:`IndexedMS1`.

  :py:func:`chain` - read multiple files at once.

  :py:func:`chain.from_iterable` - read multiple files at once, using an
  iterable of files.

  :py:func:`read_header` - get a dict with common parameters for all spectra
  from the beginning of MS1 file.

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

from . import auxiliary as aux
try:
    import numpy as np
except ImportError:
    np = None


class MS1Base(aux.ArrayConversionMixin):
    """Abstract class representing an MS1 file. Subclasses implement different approaches to parsing."""
    _array_keys = ['m/z array', 'intensity array']
    _float_keys = ['RTime', 'RetTime']

    def __init__(self, source=None, use_header=False, convert_arrays=True, dtype=None, encoding=None, **kwargs):
        """
        Create an instance of a :py:class:`MS1Base` parser.

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

        dtype : type or str or dict, optional
            dtype argument to :py:mod:`numpy` array constructor, one for all arrays or one for each key.
            Keys should be 'm/z array', 'intensity array', 'charge array'.

        encoding : str, optional
            File encoding.
        """
        super(MS1Base, self).__init__(source, use_header=use_header, convert_arrays=convert_arrays, dtype=dtype, encoding=encoding, **kwargs)
        if convert_arrays and np is None:
            raise aux.PyteomicsError('numpy is required for array conversion')
        self._use_header = use_header
        if use_header:
            self._header = self._read_header()
        else:
            self._header = None
        self._source_name = getattr(source, 'name', str(source))

    def reset(self):
        super(MS1Base, self).reset()
        self._pending_line = None

    @property
    def header(self):
        return self._header

    def _read_header_lines(self, lines):
        header = {}
        for line in lines:
            if line[0] != 'H':
                break
            tokens = line.split('\t', 2)
            if len(tokens) < 3:
                tokens = line.split(None, 2)
            key = tokens[1]
            val = tokens[2].strip()
            header[key] = val
        return header

    def _make_scan(self, info):
        for key in self._float_keys:
            if key in info['params']:
                info['params'][key] = float(info['params'][key])
        self._build_all_arrays(info)
        return info

    def _handle_S(self, line, sline, params):
        sline = line.strip().split(None, 3)
        params['scan'] = tuple(sline[1:3])
        if len(sline) == 4:  # in MS2 the S line contains the precursor m/z as a 4th column
            params['precursor m/z'] = float(sline[3])

    def _handle_I(self, line, sline, params):
        params[sline[1]] = sline[2] if len(sline) > 2 else ''

    def _handle_Z(self, line, sline, params):
        params.setdefault('charge', []).append(float(sline[1]))
        params.setdefault('neutral mass', []).append(float(sline[2]))

    def _handle_D(self, line, sline, params):
        params.setdefault('analyzer', []).append(sline[1:])

    def _handle_peak(self, line, sline, info):
        try:
            info['m/z array'].append(float(sline[0]))            # this may cause
            info['intensity array'].append(float(sline[1]))      # exceptions...
        except ValueError:
            raise aux.PyteomicsError(
                'Error when parsing %s. Line: %s' % (self._source_name, line))
        except IndexError:
            pass

    def _read_spectrum_lines(self, lines):
        params = {}
        info = {'params': params}
        for k in self._array_keys:
            info[k] = []
        if self._use_header:
            params.update(self.header)
        if self._pending_line:
            reading_spectrum = True
            self._handle_S(self._pending_line, None, params)
        else:
            reading_spectrum = False
        line_count = 0
        for i, line in enumerate(lines):
            line_count = i
            sline = line.strip().split(None, 2)
            if not sline:
                continue
            if not reading_spectrum:
                if sline[0] == 'S':
                    reading_spectrum = True
                    self._handle_S(line, sline, params)
                # otherwise we are not interested; do nothing, just move along
            else:
                if not sline:
                    pass
                elif sline[0] == 'S':
                    self._pending_line = line
                    return self._make_scan(info)

                else:
                    if sline[0] == 'I':  # spectrum-specific parameters!
                        self._handle_I(line, sline, params)
                    elif sline[0] == 'Z':  # MS2-specific charge state guess
                        self._handle_Z(line, sline, params)
                    elif sline[0] == 'D':  # MS2-specific analyzer annotation
                        self._handle_D(line, sline, params)
                    else:  # this must be a peak list
                        self._handle_peak(line, sline, info)
        self._pending_line = None
        if line_count == 0:
            return
        return self._make_scan(info)

    def __getstate__(self):
        state = super(MS1Base, self).__getstate__()
        state['use_header'] = self._use_header
        state['header'] = self._header
        return state

    def __setstate__(self, state):
        super(MS1Base, self).__setstate__(state)
        self._use_header = state['use_header']
        self._header = state['header']

    def __reduce_ex__(self, protocol):
        return (self.__class__,
            (self._source_init, False, self._convert_arrays, None, self.encoding),
            self.__getstate__())


class MS1(MS1Base, aux.FileReader):
    """
    A class representing an MS1 file. Supports the `with` syntax and direct iteration for sequential
    parsing.

    :py:class:`MS1` object behaves as an iterator, **yielding** spectra one by one.
    Each 'spectrum' is a :py:class:`dict` with three keys: 'm/z array',
    'intensity array', and 'params'. 'm/z array' and
    'intensity array' store :py:class:`numpy.ndarray`'s of floats,
    and 'params' stores a :py:class:`dict` of parameters.

    Attributes
    ----------

    header : dict
        The file header.

    """
    def __init__(self, source=None, use_header=False, convert_arrays=True, dtype=None, encoding=None, **kwargs):
        """
        Create an :py:class:`MS1` (text-mode) reader for a given MS1 file.

        Parameters
        ----------

        source : str or file or None, optional
            A file object (or file name) with data in MS1 format. Default is
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

        dtype : type or str or dict, optional
            dtype argument to :py:mod:`numpy` array constructor, one for all arrays or one for each key.
            Keys should be 'm/z array', 'intensity array', 'charge array'.

        encoding : str, optional
            File encoding.

        Returns
        -------

        out : MS1
            The reader object.
        """
        super(MS1, self).__init__(source, use_header=use_header, convert_arrays=convert_arrays, dtype=dtype, encoding=encoding,
            mode='r', parser_func=self._read, pass_file=False, args=(), kwargs={})

    @aux._keepstate_method
    def _read_header(self):
        return self._read_header_lines(self._source)

    def _read(self):
        def get_next_spectrum():
            return self._read_spectrum_lines(self._source)

        for spectrum in iter(get_next_spectrum, None):
            yield spectrum


class IndexedMS1(MS1Base, aux.TaskMappingMixin, aux.TimeOrderedIndexedReaderMixin, aux.IndexedTextReader):
    """
    A class representing an MS1 file. Supports the `with` syntax and direct iteration for sequential
    parsing. Specific spectra can be accessed by title using the indexing syntax in constant time.
    If created using a file object, it needs to be opened in binary mode.

    When iterated, :py:class:`IndexedMS1` object yields spectra one by one.
    Each 'spectrum' is a :py:class:`dict` with three keys: 'm/z array', 'intensity array' and 'params'.
    'm/z array' and 'intensity array' store :py:class:`numpy.ndarray`'s of floats,
    and 'params' stores a :py:class:`dict` of parameters (keys and values are
    :py:class:`str`, keys corresponding to MS1).

    .. warning ::
        Labels for scan objects are constructed as the first number in the S line, as follows:
        for a line ``S  0   1`` the label is `'0'`. If these labels are not unique
        for the scans in the file, the indexed parser will not work correctly. Consider using
        :py:class:`MS1` instead.

    Attributes
    ----------

    header : dict
        The file header.
    time : RTLocator
        A property used for accessing spectra by retention time.
    """

    delimiter = '\nS'
    label = r'^[\n]?S\s+(\S+)'

    def __init__(self, source=None, use_header=False, convert_arrays=True, dtype=None, encoding='utf-8', _skip_index=False, **kwargs):
        """
        Create an :py:class:`IndexedMS1` (binary-mode) reader for a given MS1 file.

        Parameters
        ----------

        source : str or file or None, optional
            A file object (or file name) with data in MS1 format. Default is
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

        dtype : type or str or dict, optional
            dtype argument to :py:mod:`numpy` array constructor, one for all arrays or one for each key.
            Keys should be 'm/z array', 'intensity array', 'charge array'.

        encoding : str, optional
            File encoding.

        block_size : int, optinal
            Size of the chunk (in bytes) used to parse the file when creating the byte offset index.

        Returns
        -------

        out : IndexedMS1
            The reader object.
        """
        super(IndexedMS1, self).__init__(source, use_header=use_header, convert_arrays=convert_arrays, dtype=dtype, encoding=encoding,
            parser_func=self._read, pass_file=False, args=(), kwargs={}, _skip_index=_skip_index, **kwargs)

    def __reduce_ex__(self, protocol):
        return (self.__class__,
            (self._source_init, False, self._convert_arrays, None, self.encoding, True),
            self.__getstate__())

    @aux._keepstate_method
    def _read_header(self):
        try:
            first = next(v for v in self._offset_index.values())[0]
        except StopIteration: # the index is empty, no spectra in file
            first = -1
        header_lines = self.read(first).decode(self.encoding).split('\n')
        return self._read_header_lines(header_lines)

    def _item_from_offsets(self, offsets):
        start, end = offsets
        lines = self._read_lines_from_offsets(start, end)
        return self._read_spectrum_lines(lines)

    def _read(self, **kwargs):
        for _, offsets in self._offset_index.items():
            spectrum = self._item_from_offsets(offsets)
            yield spectrum

    def get_spectrum(self, key):
        return self.get_by_id(key)

    def _get_time(self, spectrum):
        try:
            return spectrum['params']['RTime']
        except KeyError:
            raise aux.PyteomicsError('RT information not found.')


def read_header(source, *args, **kwargs):
    """
    Read the specified MS1 file, get the parameters specified in the header
    as a :py:class:`dict`.

    Parameters
    ----------

    source : str or file
        File name or file object representing an file in MS1 format.

    Returns
    -------

    header : dict
    """
    kwargs['use_header'] = True
    return read(source, *args, **kwargs).header


def read(*args, **kwargs):
    """Read an MS1 file and return entries iteratively.

    Read the specified MS1 file, **yield** spectra one by one.
    Each 'spectrum' is a :py:class:`dict` with three keys: 'm/z array',
    'intensity array', and 'params'. 'm/z array' and
    'intensity array' store :py:class:`numpy.ndarray`'s of floats,
    and 'params' stores a :py:class:`dict` of parameters.

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

    dtype : type or str or dict, optional
        dtype argument to :py:mod:`numpy` array constructor, one for all arrays or one for each key.
        Keys should be 'm/z array' and/or 'intensity array'.

    encoding : str, optional
        File encoding.

    use_index : bool, optional
        Determines which parsing method to use. If :py:const:`True`, an instance of
        :py:class:`IndexedMS1` is created. This facilitates random access by scan titles.
        If an open file is passed as `source`, it needs to be open in binary mode.

        If :py:const:`False` (default), an instance of :py:class:`MS1` is created. It reads
        `source` in text mode and is suitable for iterative parsing.

        .. warning ::
            Labels for scan objects are constructed as the first number in the S line, as follows:
            for a line ``S  0   1`` the label is `'0'`. If these labels are not unique
            for the scans in the file, the indexed parser will not work correctly.

    block_size : int, optinal
        Size of the chunk (in bytes) used to parse the file when creating the byte offset index.
        (Accepted only for :py:class:`IndexedMS1`.)

    Returns
    -------

    out : :py:class:`MS1Base`
        An instance of :py:class:`MS1` or :py:class:`IndexedMS1`, depending on `use_index` and `source`.
    """
    if args:
        source = args[0]
    else:
        source = kwargs.get('source')
    use_index = kwargs.pop('use_index', None)
    use_index = aux._check_use_index(source, use_index, False)
    tp = IndexedMS1 if use_index else MS1

    return tp(*args, **kwargs)


chain = aux._make_chain(read, 'read')
