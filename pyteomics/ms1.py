"""
ms1 - read and write MS/MS data in MS1 format
=============================================

Summary
-------

`MS1 <http://dx.doi.org/10.1002/rcm.1603>`_ is a simple
human-readable format for MS1 data. It allows storing MS1 peak lists and
exprimental parameters.

This module provides minimalistic infrastructure for access to data stored in
MS1 files. The most important function is :py:func:`read`, which
reads spectra and related information as saves them into human-readable
:py:class:`dicts`.
Also, common parameters can be read from MS1 file header with
:py:func:`read_header` function.

Functions
---------

  :py:func:`read` - iterate through spectra in MS1 file. Data from a
  single spectrum are converted to a human-readable dict.

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
_array_keys = ['m/z array', 'intensity array']

@aux._file_reader()
def read(source=None, use_header=False, convert_arrays=2, read_charges=True, dtype=None):
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

    convert_arrays : bool, optional
        If :py:const:`False`, m/z and intensities will be returned as regular lists.
        If :py:const:`True` (default), they will be converted to regular :py:class:`numpy.ndarray`'s.
        Conversion requires :py:mod:`numpy`.
    
    dtype : type or str or dict, optional
        dtype argument to :py:mod:`numpy` array constructor, one for all arrays or one for each key.
        Keys should be 'm/z array' and/or 'intensity array'.

    Returns
    -------

    out : FileReader
    """
    if convert_arrays and np is None:
        raise aux.PyteomicsError('numpy is required for array conversion')
    dtype_dict = dtype if isinstance(dtype, dict) else {k: dtype for k in _array_keys}
    header = read_header(source)
    reading_spectrum = False
    params = {}
    masses = []
    intensities = []
    if use_header: params.update(header)

    def make_out():
        out = {'params': params}
        if convert_arrays:
            data = {'m/z array': masses, 'intensity array': intensities}
            for key, values in data.items():
                out[key] = np.array(values, dtype=dtype_dict.get(key))
        else:
            out['m/z array'] = masses
            out['intensity array'] = intensities
        return out

    for line in source:
        sline = line.strip().split(None, 2)
        if not reading_spectrum:
            if sline[0] == 'S':
                reading_spectrum = True
                params['scan'] = tuple(sline[1:])
            # otherwise we are not interested; do nothing, just move along
        else:
            if not sline:
                pass
            elif sline[0] == 'S':
                yield make_out()
                params = dict(header) if use_header else {}
                masses = []
                intensities = []
            else:
                if sline[0] == 'I': # spectrum-specific parameters!
                    params[sline[1]] = sline[2]
                else: # this must be a peak list
                    try:
                        masses.append(float(sline[0]))            # this may cause
                        intensities.append(float(sline[1]))       # exceptions...\
                    except ValueError:
                        raise aux.PyteomicsError(
                             'Error when parsing %s. Line: %s' %
                             (source.name, line))
                    except IndexError:
                        pass

    yield make_out()

@aux._keepstate
def read_header(source):
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
    with aux._file_obj(source, 'r') as source:
        header = {}
        for line in source:
            if line[0] != 'H':
                break
            l = line.split('\t', 2)
            if len(l) < 3:
                l = line.split(None, 2)
            key = l[1]
            val = l[2].strip()
            header[key] = val
        return header

chain = aux._make_chain(read, 'read')