"""
tandem - X!Tandem output file reader
====================================

Summary
-------

`X!Tandem <http://thegpm.org/tandem/>`_  is an open-source proteomic search
engine with a very simple, sophisticated application programming interface
(API): it simply takes an XML file of instructions on its command line,
and outputs the results into an XML file, which has been specified in the input
XML file. The output format is described
`here (PDF) <http://www.thegpm.org/docs/X_series_output_form.pdf>`_.

This module provides a minimalistic way to extract information from X!Tandem
output files. You can use the old functional interface (:py:func:`read`) or the
new object-oriented interface (:py:class:`TandemXML`) to iterate over entries in
`<group>` elements, i.e. identifications for a certain spectrum.

Data access
-----------

  :py:class:`TandemXML` - a class representing a single X!Tandem output file.
  Other data access functions use this class internally.

  :py:func:`read` - iterate through peptide-spectrum matches in an X!Tandem
  output file. Data from a single PSM are converted to a human-readable dict.

  :py:func:`chain` - read multiple files at once.

  :py:func:`chain.from_iterable` - read multiple files at once, using an
  iterable of files.

  :py:func:`DataFrame` - read X!Tandem output files into a :py:class:`pandas.DataFrame`.

  :py:func:`filter` - iterate through peptide-spectrum matches in a chain of
  X!Tandem output files, yielding only top PSMs and keeping false discovery rate
  (FDR) at the desired level. The FDR is estimated using the target-decoy
  approach (TDA).

  :py:func:`filter.chain` - chain a series of filters applied independently to
  several files.

  :py:func:`filter.chain.from_iterable` - chain a series of filters applied
  independently to an iterable of files.

  :py:func:`filter_df` - filter X!Tandem output files and return a :py:class:`pandas.DataFrame`.

Miscellaneous
-------------

  :py:func:`is_decoy` - determine if a PSM is from the decoy database.

  :py:func:`fdr` - estimate the FDR in a data set using TDA.

  :py:func:`qvalues` - get an array of scores and local FDR values for a PSM
  set using the target-decoy approach.

Deprecated functions
--------------------

  :py:func:`iterfind` - iterate over elements in an X!Tandem file.
  You can just call the corresponding method of the :py:class:`TandemXML`
  object.

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

import itertools as it
import operator
from . import xml, auxiliary as aux

class TandemXML(xml.XML):
    """Parser class for TandemXML files."""
    file_format = "TandemXML"
    _root_element = "bioml"
    _default_schema = xml._tandem_schema_defaults
    _default_iter_tag = 'group[type="model"]'
    _structures_to_flatten = {'domain'}

    def __init__(self, *args, **kwargs):
        if 'recursive' not in kwargs:
            super(TandemXML, self).__init__(*args, recursive=True, **kwargs)
        else:
            super(TandemXML, self).__init__(*args, **kwargs)

    __init__.__doc__ = xml.XML.__init__.__doc__

    def _get_info_smart(self, element, **kw):
        info = self._get_info(element, **kw)
        # handy simplifications below
        if isinstance(info.get('note'), list
                ) and len(info['note']) == 1 and set(
                        info['note'][0]) == {'label', 'note'}:
            info['note'] = info['note'][0]['note']
        if 'protein' in info and 'label' in info:
            del info['label']
        if 'group' in info:
            for g in info['group']:
                label = g.pop('label')
                type_ = g.pop('type')
                info.setdefault(type_, {})[label] = g
            del info['group']
        if 'trace' in info:
            for t in info['trace']:
                info[t.pop('type')] = t
            del info['trace']
        if isinstance(info.get('values'), dict):
            info['values'] = info['values']['values']
        if isinstance(info.get('attribute'), list):
            for a in info.pop('attribute'):
                info[a['type']] = float(a['attribute'])
        if 'support' in info:
            for d in info['support'].get('supporting data', {}).values():
                for l in ['Xdata', 'Ydata']:
                    d[l]['values'] = d[l]['values'].astype(int)
            fims = info['support']['fragment ion mass spectrum']
            fims.update(fims.pop('tandem mass spectrum'))
            for d in it.chain(
                    info['support'].get('supporting data', {}).values(),
                    (info['support']['fragment ion mass spectrum'],)):
                for l in ['Xdata', 'Ydata']:
                    del d[l]['label']
        if 'charge' in info:
            info['charge'] = int(info['charge'])
        return info

    def _get_schema_info(self, read_schema):
        return self._default_schema

    def __next__(self):
        n = super(TandemXML, self).__next__()
        del n['type']
        return n

    next = __next__

def read(source, iterative=True):
    """Parse `source` and iterate through peptide-spectrum matches.

    Parameters
    ----------
    source : str or file
        A path to a target X!Tandem output file or the file object itself.

    iterative : bool, optional
        Defines whether iterative parsing should be used. It helps reduce
        memory usage at almost the same parsing speed. Default is
        :py:const:`True`.

    Returns
    -------
    out : iterator
       An iterator over dicts with PSM properties.
    """
    return TandemXML(source, read_schema=False,
            recursive=True, iterative=iterative)

def iterfind(source, path, **kwargs):
    """Parse `source` and yield info on elements with specified local
    name or by specified "XPath".

    .. note:: This function is provided for backward compatibility only.
        If you do multiple :py:func:`iterfind` calls on one file, you should
        create a :py:class:`TandemXML` object and use its
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
    return TandemXML(source, **kwargs).iterfind(path, **kwargs)

chain = aux._make_chain(read, 'read')

def is_decoy(psm, prefix='DECOY_'):
    """Given a PSM dict, return :py:const:`True` if all protein names for
    the PSM start with ``prefix``, and :py:const:`False` otherwise.

    Parameters
    ----------
    psm : dict
        A dict, as yielded by :py:func:`read`.
    prefix : str, optional
        A prefix used to mark decoy proteins. Default is `'DECOY_'`.

    Returns
    -------
    out : bool
    """
    return all(prot['label'].startswith(prefix) for prot in psm['protein'])

qvalues = aux._make_qvalues(chain, is_decoy, operator.itemgetter('expect'))
filter = aux._make_filter(chain, is_decoy, operator.itemgetter('expect'),
        qvalues)
fdr = aux._make_fdr(is_decoy)
filter.chain = aux._make_chain(filter, 'filter', True)

def DataFrame(*args, **kwargs):
    """Read X!Tandem output files into a :py:class:`pandas.DataFrame`.

    Requires :py:mod:`pandas`.

    Parameters
    ----------
    *args, **kwargs : passed to :py:func:`chain`

    sep : str or None, optional
        Some values related to PSMs (such as protein information) are variable-length
        lists. If `sep` is a :py:class:`str`, they will be packed into single string using
        this delimiter. If `sep` is :py:const:`None`, they are kept as lists. Default is
        :py:const:`None`.

    Returns
    -------
    out : pandas.DataFrame
    """
    import pandas as pd
    data = []
    prot_keys = ['id', 'uid', 'label', 'expect']
    pep_keys = ['id', 'pre', 'post', 'start', 'end']
    sep = kwargs.pop('sep', None)
    with chain(*args, **kwargs) as f:
        for item in f:
            info = {}
            for k, v in item.items():
                if isinstance(v, (str, int, float)):
                    info[k] = v
            protein = item['protein'][0]
            
            for key in prot_keys:
                vals = [prot.get(key) for prot in item['protein']]
                if sep is not None:
                    vals = sep.join(str(val) if val is not None else '' for val in vals)
                info['protein_' + key] = vals
            for key in pep_keys:
                vals = [prot['peptide'].get(key) for prot in item['protein']]
                if sep is not None:
                    vals = sep.join(str(val) if val is not None else '' for val in vals)
                info['peptide_' + key] = vals
            aa = protein['peptide'].pop('aa', [])
            info['modifications'] = ','.join('{0[modified]:.3f}@{0[type]}'.format(x) for x in aa)
            for k in prot_keys:
                protein.pop(k, None)
            for k in pep_keys:
                protein['peptide'].pop(k, None)
            del protein['peptide']['peptide']
            info.update(protein['peptide'])
            info['scan'] = item['support']['fragment ion mass spectrum']['note']
            data.append(info)
    return pd.DataFrame(data)

def filter_df(*args, **kwargs):
    """Read X!Tandem output files or DataFrames and return a :py:class:`DataFrame` with filtered PSMs.
    Positional arguments can be X!Tandem output files or DataFrames.

    Requires :py:mod:`pandas`.

    Parameters
    ----------
    key : str / iterable / callable, optional
        Default is 'expect'.
    is_decoy : str / iterable / callable, optional
        Default is to check if all strings in the "protein" column start with `'DECOY_'`
    *args, **kwargs : passed to :py:func:`auxiliary.filter` and/or :py:func:`DataFrame`.

    Returns
    -------
    out : pandas.DataFrame
    """
    import pandas as pd
    sep = kwargs.get('sep')
    kwargs.setdefault('key', 'expect')
    if all(isinstance(arg, pd.DataFrame) for arg in args):
        df = pd.concat(args)
    else:
        read_kw = {k: kwargs.pop(k) for k in ['iterative', 'read_schema', 'sep'] if k in kwargs}
        df = DataFrame(*args, **read_kw)
    if sep is not None:
        kwargs.setdefault('is_decoy',
            df['protein_label'].str.split(sep).apply(lambda s: all(x.startswith('DECOY') for x in s)))
    else:
        kwargs.setdefault('is_decoy',
            df['protein_label'].apply(lambda s: all(x.startswith('DECOY') for x in s)))
    return aux.filter(df, **kwargs)