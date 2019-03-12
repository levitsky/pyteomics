"""
mgf - read and write MS/MS data in Mascot Generic Format
========================================================

Summary
-------

`MGF <http://www.matrixscience.com/help/data_file_help.html>`_ is a simple
human-readable format for MS/MS data. It allows storing MS/MS peak lists and
exprimental parameters.

This module provides classes and functions for access to data stored in
MGF files.
Parsing is done using :py:class:`MGF` and :py:class:`IndexedMGF` classes.
The :py:func:`read` function can be used as an entry point.
MGF spectra are converted to dictionaries. MS/MS data points are
(optionally) represented as :py:mod:`numpy` arrays.
Also, common parameters can be read from MGF file header with
:py:func:`read_header` function.
:py:func:`write` allows creation of MGF files.

Classes
-------

  :py:class:`MGF` - a text-mode MGF parser. Suitable to read spectra from a file consecutively.
  Needs a file opened in text mode (or will open it if given a file name).

  :py:class:`IndexedMGF` - a binary-mode MGF parser. When created, builds a byte offset index
  for fast random access by spectrum titles. Sequential iteration is also supported.
  Needs a seekable file opened in binary mode (if created from existing file object).

  :py:class:`MGFBase` - abstract class, the common ancestor of the two classes above.
  Can be used for type checking.

Functions
---------

  :py:func:`read` - iterate through spectra in MGF file. Data from a
  single spectrum are converted to a human-readable dict.

  :py:func:`get_spectrum` - read a single spectrum with given title from a file.

  :py:func:`chain` - read multiple files at once.

  :py:func:`chain.from_iterable` - read multiple files at once, using an
  iterable of files.

  :py:func:`read_header` - get a dict with common parameters for all spectra
  from the beginning of MGF file.

  :py:func:`write` - write an MGF file.

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

try:
    import numpy as np
except ImportError:
    np = None
import itertools as it
import sys
from . import auxiliary as aux

class MGFBase():
    """Abstract class representing an MGF file. Subclasses implement different approaches to parsing."""
    _comments = set('#;!/')
    _array = (lambda x, dtype: np.array(x, dtype=dtype)) if np is not None else None
    _ma = (lambda x, dtype: np.ma.masked_equal(np.array(x, dtype=dtype), 0)) if np is not None else None
    _identity = lambda x, **kw: x
    _array_converters = {
        'm/z array': [_identity, _array, _array],
        'intensity array': [_identity, _array, _array],
        'charge array': [_identity, _array, _ma]
    }
    _array_keys = ['m/z array', 'intensity array', 'charge array']
    _array_keys_unicode = [u'm/z array', u'intensity array', u'charge array']
    encoding = None


    def __init__(self, source=None, use_header=True, convert_arrays=2, read_charges=True, dtype=None):
        """Create an MGF file object, set MGF-specific parameters.

        Parameters
        ----------

        source : str or file or None, optional
            A file object (or file name) with data in MGF format. Default is
            :py:const:`None`, which means read standard input.

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

        dtype : type or str or dict, optional
            dtype argument to :py:mod:`numpy` array constructor, one for all arrays or one for each key.
            Keys should be 'm/z array', 'intensity array' and/or 'charge array'.

        encoding : str, optional
            File encoding.
        """

        self._use_header = use_header
        self._convert_arrays = convert_arrays
        if self._convert_arrays and np is None:
            raise aux.PyteomicsError('numpy is required for array conversion')
        self._read_charges = read_charges
        self._dtype_dict = dtype if isinstance(dtype, dict) else {k: dtype for k in self._array_keys}
        if self._use_header:
            self._read_header()
        else:
            self._header = None

    def parse_charge(self, charge_text, list_only=False):
        return aux._parse_charge(charge_text, list_only=list_only)

    @property
    def header(self):
        if self._header is None:
            self._read_header()
        return self._header

    def _read_header_lines(self, header_lines):
        header = {}
        for line in header_lines:
            if line.strip() == 'BEGIN IONS':
                break
            l = line.split('=')
            if len(l) == 2:
                key = l[0].lower()
                val = l[1].strip()
                header[key] = val
        if 'charge' in header:
            header['charge'] = self.parse_charge(header['charge'], True)
        self._header = header

    def _read_spectrum_lines(self, lines):
        """Read a single spectrum from ``self._source``.

        Returns
        -------
        out : dict
        """

        masses = []
        intensities = []
        charges = []

        params = self.header.copy() if self._use_header else {}

        for i, line in enumerate(lines):
            sline = line.strip()
            if sline == 'BEGIN IONS':
                if i == 0:
                    continue
                else:
                    raise aux.PyteomicsError('Error when parsing MGF: unexpected start of spectrum.')
            if not sline or sline[0] in self._comments:
                pass
            elif sline == 'END IONS':
                if 'pepmass' in params:
                    try:
                        pepmass = tuple(map(float, params['pepmass'].split()))
                    except ValueError:
                        raise aux.PyteomicsError('MGF format error: cannot parse '
                                'PEPMASS = {}'.format(params['pepmass']))
                    else:
                        params['pepmass'] = pepmass + (None,) * (2-len(pepmass))
                if isinstance(params.get('charge'), aux.basestring):
                    params['charge'] = self.parse_charge(params['charge'], True)
                if 'rtinseconds' in params:
                    params['rtinseconds'] = aux.unitfloat(params['rtinseconds'], 'second')
                out = {'params': params}
                data = {'m/z array': masses, 'intensity array': intensities}
                if self._read_charges:
                    data['charge array'] = charges
                for key, values in data.items():
                    out[key] = self._array_converters[key][self._convert_arrays](values, dtype=self._dtype_dict.get(key))
                if self.encoding and sys.version_info.major == 2:
                    for key, ukey in zip(self._array_keys + ['params'], self._array_keys_unicode + [u'params']):
                        if key in out:
                            out[ukey] = out.pop(key)
                return out

            else:
                if '=' in sline: # spectrum-specific parameters!
                    l = sline.split('=', 1)
                    params[l[0].lower()] = l[1].strip()
                else: # this must be a peak list
                    l = sline.split()
                    try:
                        masses.append(float(l[0]))
                        intensities.append(float(l[1]))
                        if self._read_charges:
                            charges.append(self.parse_charge(l[2]) if len(l) > 2 else 0)
                    except ValueError:
                        raise aux.PyteomicsError(
                             'Error when parsing %s. Line:\n%s' % (getattr(self._source, 'name', 'MGF file'), line))
                    except IndexError:
                        pass

    def get_spectrum(self, title):
        raise NotImplementedError

    def __getitem__(self, key):
        return self.get_spectrum(key)


class IndexedMGF(aux.TaskMappingMixin, aux.TimeOrderedIndexedReaderMixin, aux.IndexSavingTextReader, MGFBase):
    """
    A class representing an MGF file. Supports the `with` syntax and direct iteration for sequential
    parsing. Specific spectra can be accessed by title using the indexing syntax in constant time.
    If created using a file object, it needs to be opened in binary mode.

    When iterated, :py:class:`IndexedMGF` object yields spectra one by one.
    Each 'spectrum' is a :py:class:`dict` with four keys: 'm/z array',
    'intensity array', 'charge array' and 'params'. 'm/z array' and
    'intensity array' store :py:class:`numpy.ndarray`'s of floats,
    'charge array' is a masked array (:py:class:`numpy.ma.MaskedArray`) of ints,
    and 'params' stores a :py:class:`dict` of parameters (keys and values are
    :py:class:`str`, keys corresponding to MGF, lowercased).


    Attributes
    ----------

    header : dict
        The file header.
    time : RTLocator
        A property used for accessing spectra by retention time.
    """

    delimiter = 'BEGIN IONS'
    label = r'TITLE=([^\n]*\S)\s*'

    def __init__(self, source=None, use_header=True, convert_arrays=2, read_charges=True,
        dtype=None, encoding='utf-8', block_size=1000000, _skip_index=False):
        aux.TimeOrderedIndexedReaderMixin.__init__(self, source, self._read, False, (), {}, encoding,
            block_size, _skip_index=_skip_index)
        MGFBase.__init__(self, source, use_header, convert_arrays, read_charges, dtype)

    def __reduce_ex__(self, protocol):
        return (self.__class__,
            (self._source_init, False, self._convert_arrays,
                self._read_charges, self._dtype_dict, self.encoding, self.block_size, True),
            self.__getstate__())

    def __getstate__(self):
        state = super(IndexedMGF, self).__getstate__()
        state['use_header'] = self._use_header
        state['header'] = self._header
        return state

    def __setstate__(self, state):
        super(IndexedMGF, self).__setstate__(state)
        self._use_header = state['use_header']
        self._header = state['header']

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
            return spectrum['params']['rtinseconds']
        except KeyError:
            raise aux.PyteomicsError('RT information not found.')


class MGF(aux.FileReader, MGFBase):
    """
    A class representing an MGF file. Supports the `with` syntax and direct iteration for sequential
    parsing. Specific spectra can be accessed by title using the indexing syntax (if the file is seekable),
    but it takes linear time to search through the file. Consider using :py:class:`IndexedMGF` for
    constant-time access to spectra.

    :py:class:`MGF` object behaves as an iterator, **yielding** spectra one by one.
    Each 'spectrum' is a :py:class:`dict` with four keys: 'm/z array',
    'intensity array', 'charge array' and 'params'. 'm/z array' and
    'intensity array' store :py:class:`numpy.ndarray`'s of floats,
    'charge array' is a masked array (:py:class:`numpy.ma.MaskedArray`) of ints,
    and 'params' stores a :py:class:`dict` of parameters (keys and values are
    :py:class:`str`, keys corresponding to MGF, lowercased).

    Attributes
    ----------

    header : dict
        The file header.

    """

    def __init__(self, source=None, use_header=True, convert_arrays=2, read_charges=True,
        dtype=None, encoding=None):
        aux.FileReader.__init__(self, source, 'r', self._read, False, (), {}, encoding)
        MGFBase.__init__(self, source, use_header, convert_arrays, read_charges, dtype)
        self.encoding = encoding

    @aux._keepstate_method
    def _read_header(self):
        return self._read_header_lines(self._source)

    def _read_spectrum(self):
        return self._read_spectrum_lines(self._source)

    def _read(self):
        for line in self._source:
            if line.strip() == 'BEGIN IONS':
                yield self._read_spectrum()

    @aux._keepstate_method
    def get_spectrum(self, title):
        for line in self._source:
            sline = line.strip()
            if sline[:5] == 'TITLE' and sline.split('=', 1)[1].strip() == title:
                spectrum = self._read_spectrum()
                spectrum['params']['title'] = title
                return spectrum


def read(*args, **kwargs):
    """Returns a reader for a given MGF file. Most of the parameters repeat the
    instantiation signature of :py:class:`MGF` and :py:class:`IndexedMGF`.
    Additional parameter `use_index` helps decide which class to instantiate
    for given `source`.

    Parameters
    ----------

    source : str or file or None, optional
        A file object (or file name) with data in MGF format. Default is
        :py:const:`None`, which means read standard input.

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

    dtype : type or str or dict, optional
        dtype argument to :py:mod:`numpy` array constructor, one for all arrays or one for each key.
        Keys should be 'm/z array', 'intensity array' and/or 'charge array'.

    encoding : str, optional
        File encoding.

    use_index : bool, optional
        Determines which parsing method to use. If :py:const:`True` (default), an instance of
        :py:class:`IndexedMGF` is created. This facilitates random access by spectrum titles.
        If an open file is passed as `source`, it needs to be open in binary mode.

        If :py:const:`False`, an instance of :py:class:`MGF` is created. It reads
        `source` in text mode and is suitable for iterative parsing. Access by spectrum title
        requires linear search and thus takes linear time.

    block_size : int, optinal
        Size of the chunk (in bytes) used to parse the file when creating the byte offset index.
        (Accepted only for :py:class:`IndexedMGF`.)

    Returns
    -------

    out : MGFBase
        Instance of :py:class:`MGF` or :py:class:`IndexedMGF`.
    """
    if args:
        source = args[0]
    else:
        source = kwargs.get('source')
    use_index = kwargs.pop('use_index', None)
    use_index = aux._check_use_index(source, use_index, True)
    tp = IndexedMGF if use_index else MGF

    return tp(*args, **kwargs)


def get_spectrum(source, title, *args, **kwargs):
    """Read one spectrum (with given `title`) from `source`.

    See :py:func:`read` for explanation of parameters affecting the output.

    .. note :: Only the key-value pairs after the "TITLE =" line will be included in the output.

    Parameters
    ----------

    source : str or file or None
        File to read from.
    title : str
        Spectrum title.

    The rest of the arguments are the same as for :py:func:`read`.

    Returns
    -------
    out : dict or None
        A dict with the spectrum, if it is found, and None otherwise.

    """
    with read(source, *args, **kwargs) as f:
        return f[title]


@aux._keepstate
def read_header(source):
    """
    Read the specified MGF file, get search parameters specified in the header
    as a :py:class:`dict`, the keys corresponding to MGF format (lowercased).

    Parameters
    ----------

    source : str or file
        File name or file object representing an file in MGF format.

    Returns
    -------

    header : dict
    """
    with aux._file_obj(source, 'r') as source:
        header = {}
        for line in source:
            if line.strip() == 'BEGIN IONS':
                break
            l = line.split('=')
            if len(l) == 2:
                key = l[0].lower()
                val = l[1].strip()
                header[key] = val
        if 'charge' in header:
            header['charge'] = aux._parse_charge(header['charge'], True)
        return header


_default_key_order = ['title', 'pepmass', 'rtinseconds', 'charge']

def _pepmass_repr(k, pepmass):
    outstr = k.upper() + '='
    if not isinstance(pepmass, (str, int, float)): # assume iterable
        try:
            outstr += ' '.join(str(x) for x in pepmass if x is not None)
        except TypeError:
            raise aux.PyteomicsError(
                    'Cannot handle parameter: PEPMASS = {}'.format(pepmass))
    else:
        outstr += str(pepmass)
    return outstr

def _charge_repr(k, charge):
    return '{}={}'.format(k.upper(), aux._parse_charge(str(charge)))

def _default_repr(key, val):
    return '{}={}'.format(key.upper(), val)

_default_value_formatters = {'pepmass': _pepmass_repr, 'charge': _charge_repr}

@aux._file_writer()
def write(spectra, output=None, header='', key_order=_default_key_order,
    fragment_format=None, write_charges=True, use_numpy=None,
    param_formatters=_default_value_formatters):
    """
    Create a file in MGF format.

    Parameters
    ----------

    spectra : iterable
        A sequence of dictionaries with keys 'm/z array', 'intensity array',
        and 'params'. 'm/z array' and 'intensity array' should be sequences of
        :py:class:`int`, :py:class:`float`, or :py:class:`str`. Strings will
        be written 'as is'. The sequences should be of equal length, otherwise
        excessive values will be ignored.

        'params' should be a :py:class:`dict` with keys corresponding to MGF
        format. Keys must be strings, they will be uppercased and used as is,
        without any format consistency tests. Values can be of any type allowing
        string representation.

        'charge array' can also be specified.

    output : str or file or None, optional
        Path or a file-like object open for writing. If an existing file is
        specified by file name, it will be opened for appending. In this case
        writing with a header can result in violation of format conventions.
        Default value is :py:const:`None`, which means using standard output.

    header : dict or (multiline) str or list of str, optional
        In case of a single string or a list of strings, the header will be
        written 'as is'. In case of dict, the keys (must be strings) will be
        uppercased.

    write_charges : bool, optional
        If :py:const:`False`, fragment charges from 'charge array' will not be written.
        Default is :py:const:`True`.

    fragment_format : str, optional
        Format string for m/z, intensity and charge of a fragment. Useful to set
        the number of decimal places, e.g.:
        ``fragment_format='%.4f %.0f'``. Default is ``'{} {} {}'``.

        .. note::
            The supported format syntax differs depending on other parameters.
            If `use_numpy` is :py:const:`True` and :py:mod:`numpy` is available,
            fragment peaks will be written using :py:func:`numpy.savetxt`. Then,
            `fragment_format` must be recognized by that function.

            Otherwise, plain Python string formatting is done.
            See `the docs
            <https://docs.python.org/library/string.html#format-specification-mini-language>`_
            for details on writing the format string.
            If some or all charges are missing, an empty string is substituted
            instead, so formatting as :py:class:`!float` or :py:class:`!int` will raise an exception.
            Hence it is safer to just use ``{}`` for charges.

    key_order : list, optional
        A list of strings specifying the order in which params will be written in
        the spectrum header. Unlisted keys will be in arbitrary order.
        Default is :py:data:`_default_key_order`.

        .. note:: This does not affect the order of lines in the global header.

    param_formatters : dict, optional
        A dict mapping parameter names to functions. Each function must accept
        two arguments (key and value) and return a string.
        Default is :py:data:`_default_value_formatters`.

    use_numpy : bool, optional
        Controls whether fragment peak arrays are written using :py:func:`numpy.savetxt`.
        Using :py:func:`numpy.savetxt` is faster, but cannot handle sparse arrays of fragment charges.
        You may want to disable this if you need to save spectra with 'charge arrays' with missing values.

        If not specified, will be set to the opposite of `write_chrages`.
        If :py:mod:`numpy` is not available, this parameter has no effect.

    file_mode : str, keyword only, optional
        If `output` is a file name, defines the mode the file will be opened in.
        Otherwise will be ignored. Default is 'a'.

    encoding : str, keyword only, optional
        Output file encoding (if `output` is specified by name).

    Returns
    -------

    output : file
    """
    def key_value_line(key, val):
        return param_formatters.get(key, _default_repr)(key, val) + '\n'

    nones = (None, np.nan, np.ma.masked) if np is not None else (None,)

    if fragment_format is None:
        fragment_format = '{} {} {}'
        np_format_2 = '%.5f %.1f'
        np_format_3 = '%.5f %.1f %d'
    else:
        np_format_2 = np_format_3 = fragment_format
    format_str = fragment_format + '\n'

    if use_numpy is None:
        use_numpy = not write_charges

    if isinstance(header, dict):
        head_dict = header.copy()
        head_lines = [key_value_line(k, v) for k, v in header.items()]
        head_str = '\n'.join(head_lines)
    else:
        if isinstance(header, str):
            head_str = header
            head_lines = header.split('\n')
        else:
            head_lines = list(header)
            head_str = '\n'.join(header)
        head_dict = {}
        for line in head_lines:
            if not line.strip() or any(
                line.startswith(c) for c in MGF._comments):
               continue
            l = line.split('=')
            if len(l) == 2:
                head_dict[l[0].lower()] = l[1].strip()
    if head_str:
        output.write(head_str + '\n\n')

    for spectrum in spectra:
        output.write('BEGIN IONS\n')
        found = set()
        for key in it.chain(key_order, spectrum['params']):
            if key not in found and key in spectrum['params']:
                found.add(key)
                val = spectrum['params'][key]
                if val != head_dict.get(key):
                    output.write(key_value_line(key, val))

        try:
            success = True
            if np is not None and use_numpy:
                if not write_charges or 'charge array' not in spectrum:
                    X = np.empty((len(spectrum['m/z array']), 2))
                    X[:, 0] = spectrum['m/z array']
                    X[:, 1] = spectrum['intensity array']
                    np.savetxt(output, X, fmt=np_format_2)
                elif isinstance(spectrum.get('charge array'), np.ndarray):
                    X = np.empty((len(spectrum['m/z array']), 3))
                    X[:, 0] = spectrum['m/z array']
                    X[:, 1] = spectrum['intensity array']
                    X[:, 2] = spectrum['charge array']
                    np.savetxt(output, X, fmt=np_format_3)
                else:
                    success = False
            else:
                success = False

            if not success:
                for m, i, c in zip(spectrum['m/z array'],
                        spectrum['intensity array'],
                        spectrum.get('charge array', it.cycle((None,)))):
                    output.write(format_str.format(
                        m, i,
                        (c if c not in nones else '')))
        except KeyError:
            raise aux.PyteomicsError("'m/z array' and 'intensity array' must be present in all spectra.")
        output.write('END IONS\n\n')
    return output

chain = aux._make_chain(read, 'read')