"""
mgf - read and write MS/MS data in Mascot Generic Format. 
=========================================================

Summary
-------

`MGF <http://www.matrixscience.com/help/data_file_help.html>`_ is a simple
human-readable format for MS/MS data. It allows storing MS/MS peak lists and
exprimental parameters.

This module provides minimalistic infrastructure for access to data stored in
MGF files. The most important function is :py:func:`iter_spectrum`, which 
reads spectra and related information as saves them into human-readable
:py:class:`dict`'s.
Also, common parameters can be read from MGF file header with
:py:func:`read_header` function. :py:func:`write_mgf` allows creation of MGF
files.

Functions
---------

  :py:func:`iter_spectrum` - iterate through spectra in MGF file. Data from a
  single spectrum are converted to a human-readable dict.

  :py:func:`read_header` - get a dict with common parameters for all spectra
  from the beginning of MGF file.

  :py:func:`write_mgf` - write an MGF file.

-------------------------------------------------------------------------------

"""
# Licensed under the MIT license:
# http://www.opensource.org/licenses/mit-license.php

import sys
import os
from auxiliary import PyteomicsError
import numpy


_comments = '#;!/'

def iter_spectrum(source, use_header=True):
    """Read an MGF file and return entries iteratively.
    
    Read the specified MGF file, **yield** spectra one by one.
    Each 'spectrum' is a :py:class:`dict` with three keys: 'masses', 'intensities'
    and 'params'. 'masses' and 'intensities' store :py:class:`numpy.array`'s of
    floats, and 'params' stores a :py:class:`dict` of parameters (keys and values are
    :py:class:`str`, keys corresponding to MGF, lowercased).

    Parameters
    ----------

    source : str or file
        A file object (or file name) with data in MGF format.
    
    use_header : bool, optional
        Add the info from file header to each dict. Spectrum-specific parameters
        override those from the header in case of conflict.
        Default is True.

    Returns
    -------
    dict : {'masses': mumpy.array, 'intensities': numpy.asrray, 'params': dict} 
    """
    if type(source) == file:
        pos = source.tell()
        header = read_header(source, False)
        source.seek(pos)
        MGF = source
    elif type(source) == str:
        header = read_header(source)
        MGF = open(source)
    else:
        raise PyteomicsError("Unsupported argument type in `pyteomics.mgf.iter_spectrum`."
                "'source' must be a file object or a path (string), %s given." % type(source))
    reading_spectrum = False
    params = {}
    masses = []
    intensities = []
    if use_header: params.update(header)
    for line in MGF:
        if not reading_spectrum:
            if line.strip() == 'BEGIN IONS':
                reading_spectrum = True
            # otherwise we are not interested; do nothing, just move along
        else:
            if line.strip() == '' or \
                any([line.startswith(c) for c in _comments]):
                    pass
            elif line.strip() == 'END IONS':
                reading_spectrum = False
                yield {'params': params, 'masses': numpy.array(masses),
                        'intensities': numpy.array(intensities)}
                params = {}
                if use_header: params.update(header)
                masses = []
                intensities = []
                continue
            else: 
                l = line.split('=')
                if len(l) == 2:
                    params[l[0].lower()] = l[1].strip()
                elif len(l) == 1: # this must be a peak list
                    l = line.split()
                    if len(l) == 2:
                        masses.append(float(l[0]))            # this may cause
                        intensities.append(float(l[1]))       # exceptions...
                        # ... if this is not a peak list

def read_header(source, close=True):
    """
    Read the specified MGF file, get search parameters specified in the header
    as a :py:class:`dict`, the keys corresponding to MGF format (lowercased).
    
    Parameters
    ----------

    source : str or file
        File name or file object representing an file in MGF format.

    close : bool, optional
        Defines whether the output file should be closed in the end. 

    Returns
    -------

    header : dict
    """
    if type(source) == file:
        MGF = source
    elif type(source) == str:
        MGF = open(source)
    else:
        raise PyteomicsError("Unsupported argument `mgf` in `pyteomics.mgf.read_header`."
                "Must be a file object or a path (str), %s given" % type(source))
    header = {}
    for line in MGF:
        if line.strip() == 'BEGIN IONS':
            break
        l = line.split('=')
        if len(l) == 2:
            header[l[0].lower()] = l[1].strip()
    if close: MGF.close()
    return header

def write_mgf(output=None, spectra=None, header='', close=True):
    """
    Create a file in MGF format.

    Parameters
    ----------

    output : str or file, optional
        Path or a file-like object open for writing. If an existing file is
        specified by file name, it will be opened for appending. In this case
        writing with a header can result in violation of format conventions.
        Default value is None, which means using standard output.

    spectra : a list of dicts, optional
        A list of dictionaries with keys 'masses', 'intensities', and 'params'.
        'masses' and 'intensities' should be sequences of :py:class:`int`,
        :py:class:`float`, or :py:class:`str`. Strings will be written 'as is'.
        The sequences should be of equal length, otherwise excessive values will
        be ignored.

        'params' should be a :py:class:`dict` with keys corresponding to MGF
        format. Keys must be strings, they will be uppercased and used as is,
        without any format consistency tests. Values can be of any type allowing
        string representation.

        Default value is :py:const:`None`, which means no spectra will be written, only a
        header.

    header : dict or (multiline) str or list of str, optional
        In case of a single string or a list of strings, the header will be
        written 'as is'. In case of dict, the keys (must be strings) will be
        uppercased, the values will be written as strings, also uppercased.

    close : bool, optional
        Defines whether the output file should be closed in the end. If `output`
        is None, the `close` option is ignored (standard output will not be
        closed).

    Returns
    -------

    output : file
    """
    if type(output) == file:
        MGF = output
    elif type(output) == str:
        if not os.path.isdir(os.path.split(mgf)[0]):
            os.makedirs(os.path.split(mgf)[0])
        MGF = open(output, 'a')
    elif output == None:
        MGF = sys.stdout
    else:
        raise PyteomicsError("Unsupported argument `output` in `pyteomics.mgf.write_mgf`."
                "Must be a file object or a path (string) or None, not %s" % type(output))

    if type(header) == dict:
        head_dict = header
        head_str = '\n'.join(
                ['%s=%s' % (x.upper(), str(header[x])) for x in header])
    else:
        if type(header) == str:
            head_str = header
            head_lines = header.split('\n')
        elif type(header) == list:
            head_lines = header
            head_str = '\n'.join(header)
        head_dict = {}
        for line in head_lines:
            if line.strip() == '' or \
                any([line.startswith(c) for c in _comments]):
                   continue 
        l = line.split('=')
        if len(l) == 2:
            head_dict[l[0].lower()] = l[1].strip()
    MGF.write(head_str)

    if spectra:
        for spectrum in spectra:
            MGF.write('\n\nBEGIN IONS\n')
            MGF.write('\n'.join(['%s=%s' % (x.upper(), str(spectrum['params'][x]))
                for x in spectrum['params'] if not
                (x in head_dict and spectrum['params'][x] == head_dict[x])]))
            if 'masses' in spectrum and 'intensities' in spectrum:
                for i in range(min(map(len,
                    (spectrum['masses'], spectrum['intensities'])))):
                        MGF.write('\n%s %s' % (str(spectrum['masses'][i]),
                            str(spectrum['intensities'][i])))

            MGF.write('\nEND IONS')
    MGF.write('\n')
    if close and output:
        MGF.close()
    return MGF
