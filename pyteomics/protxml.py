"""
protxml - parsing of ProteinProphet output files
================================================

Summary
-------

**protXML** is the output format of the `ProteinProphet software <http://proteinprophet.sourceforge.net/>`_.
It contains information about identified proteins and their statistical significance.

This module provides minimalistic infrastructure for access to data stored in
protXML files. The central class is :py:class:`ProtXML`, which
reads protein entries and related information and saves them into
Python dicts.

Data access
-----------

  :py:class:`ProtXML` - a class representing a single protXML file.
  Other data access functions use this class internally.

  :py:func:`read` - iterate through peptide-spectrum matches in a protXML
  file. Calling the function is synonymous to instantiating the :py:class:`ProtXML` class.

  :py:func:`chain` - read multiple files at once.

  :py:func:`chain.from_iterable` - read multiple files at once, using an
  iterable of files.

  :py:func:`DataFrame` - read protXML files into a :py:class:`pandas.DataFrame`.

Target-decoy approach
---------------------

  :py:func:`filter` - filter protein groups from a chain of protXML files to a specific FDR
  using TDA.

  :py:func:`filter.chain` - chain a series of filters applied independently to
  several files.

  :py:func:`filter.chain.from_iterable` - chain a series of filters applied
  independently to an iterable of files.

  :py:func:`filter_df` - filter protXML files and return a :py:class:`pandas.DataFrame`.

  :py:func:`fdr` - estimate the false discovery rate of a set of protein groups using the
  target-decoy approach.

  :py:func:`qvalues` - get an array of scores and *q* values for protein groups using the target-decoy approach.

  :py:func:`is_decoy` - determine whether a protein group is decoy or not. This function may not suit your use case.

Dependencies
------------

This module requres :py:mod:`lxml`.

--------------------------------------------------------------------------------
"""

#   Copyright 2018 Lev Levitsky
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

from . import xml, auxiliary as aux, _schema_defaults
import operator as op

class ProtXML(xml.XML):
    """Parser class for protXML files."""
    file_format = 'protXML'
    _root_element = 'protein_summary'
    _default_schema = _schema_defaults._protxml_schema_defaults
    # _default_version = None
    _default_iter_tag = 'protein_group'
    _structures_to_flatten = {'annotation'}
    # attributes which contain unconverted values
    _convert_items = {'float':  {'pct_spectrum_ids'},
        'int': {'group_number', 'prot_length'},
        'bool': {'is_contributing_evidence', 'is_nondegenerate_evidence'}
        }.items()

    def _get_info_smart(self, element, **kwargs):
        """Extract the info in a smart way depending on the element type"""
        try:
            name = kwargs.pop('ename')
        except KeyError:
            name = xml._local_name(element)
        rec = kwargs.pop('recursive', None)
        if name == self._root_element:
            info = self._get_info(element, ename=name,
                    recursive=(rec if rec is not None else False),
                    **kwargs)
        else:
            info = self._get_info(element, ename=name,
                    recursive=(rec if rec is not None else True),
                    **kwargs)

        converters = {'float': float, 'int': int,
                'bool': lambda x: x.lower() in {'1', 'true', 'y'}}
        for k, v in dict(info).items():
            for t, s in self._convert_items:
                if k in s:
                    del info[k]
                    info[k] = converters[t](v)
        p = info.get('parameter')
        if isinstance(p, list) and len(p) == 1 and isinstance(p[0], dict):
            info.update(info.pop('parameter')[0])

        if 'modification_info' in info:
            # this is a list with one element
            info.update(info.pop('modification_info')[0])

        if 'unique_stripped_peptides' in info:
            info['unique_stripped_peptides'] = info['unique_stripped_peptides'].split('+')
        return info

def read(source, read_schema=False, iterative=True, **kwargs):
    """Parse `source` and iterate through protein groups.

    Parameters
    ----------
    source : str or file
        A path to a target protXML file or the file object itself.

    read_schema : bool, optional
        If :py:const:`True`, attempt to extract information from the XML schema
        mentioned in the protXML header. Otherwise, use default parameters.
        Not recommended without Internet connection or
        if you don't like to get the related warnings.

    iterative : bool, optional
        Defines whether iterative parsing should be used. It helps reduce
        memory usage at almost the same parsing speed. Default is
        :py:const:`True`.

    Returns
    -------
    out : ProtXML
       An iterator over dicts with protein group properties.
    """

    return ProtXML(source, read_schema=read_schema, iterative=iterative)

chain = aux._make_chain(read, 'read')

def _is_decoy_prefix(pg, prefix='DECOY_'):
    """Determine if a protein group should be considered decoy.

    This function checks that all protein names in a group start with `prefix`.
    You may need to provide your own function for correct filtering and FDR estimation.

    Parameters
    ----------

    pg : dict
        A protein group dict produced by the :py:class:`ProtXML` parser.
    prefix : str, optional
        A prefix used to mark decoy proteins. Default is `'DECOY_'`.

    Returns
    -------

    out : bool
    """
    return all(p['protein_name'].startswith(prefix) for p in pg['protein'])

def _is_decoy_suffix(pg, suffix='_DECOY'):
    """Determine if a protein group should be considered decoy.

    This function checks that all protein names in a group end with `suffix`.
    You may need to provide your own function for correct filtering and FDR estimation.

    Parameters
    ----------

    pg : dict
        A protein group dict produced by the :py:class:`ProtXML` parser.
    suffix : str, optional
        A suffix used to mark decoy proteins. Default is `'_DECOY'`.

    Returns
    -------

    out : bool
    """
    return all(p['protein_name'].endswith(suffix) for p in pg['protein'])

is_decoy = _is_decoy_prefix

fdr = aux._make_fdr(_is_decoy_prefix, _is_decoy_suffix)
_key = op.itemgetter('probability')
qvalues = aux._make_qvalues(chain, _is_decoy_prefix, _is_decoy_suffix, _key)
filter = aux._make_filter(chain, _is_decoy_prefix, _is_decoy_suffix, _key, qvalues)
filter.chain = aux._make_chain(filter, 'filter', True)

def DataFrame(*args, **kwargs):
    """Read protXML output files into a :py:class:`pandas.DataFrame`.

    .. note :: Rows in the DataFrame correspond to individual proteins, not protein groups.

    Requires :py:mod:`pandas`.

    Parameters
    ----------
    *args, **kwargs : passed to :py:func:`chain`

    sep : str or None, optional
        Some values related to protein groups are variable-length lists.
        If `sep` is a :py:class:`str`, they will be packed into single string using
        this delimiter. If `sep` is :py:const:`None`, they are kept as lists. Default is
        :py:const:`None`.

    pd_kwargs : dict, optional
        Keyword arguments passed to the :py:class:`pandas.DataFrame` constructor.

    Returns
    -------
    out : pandas.DataFrame
    """
    import pandas as pd
    kwargs = kwargs.copy()
    sep = kwargs.pop('sep', None)
    pd_kwargs = kwargs.pop('pd_kwargs', {})
    def gen_items():
        with chain(*args, **kwargs) as f:
            for item in f:
                info = {}
                for k, v in item.items():
                    if isinstance(v, (str, int, float)):
                        info[k] = v
                if 'protein' in item:
                    for prot in item['protein']:
                        out = dict(info)
                        out.update(prot)
                        if 'unique_stripped_peptides' in out:
                            if sep is not None:
                                out['unique_stripped_peptides'] = sep.join(out['unique_stripped_peptides'])
                        if 'indistinguishable_protein' in out:
                            if sep is None:
                                out['indistinguishable_protein'] = [p['protein_name'] for p in out['indistinguishable_protein']]
                            else:
                                out['indistinguishable_protein'] = sep.join(p['protein_name'] for p in out['indistinguishable_protein'])
                        yield out
    return pd.DataFrame(gen_items(), **pd_kwargs)

def filter_df(*args, **kwargs):
    """Read protXML files or DataFrames and return a :py:class:`DataFrame` with filtered PSMs.
    Positional arguments can be protXML files or DataFrames.

    .. note :: Rows in the DataFrame correspond to individual proteins, not protein groups.

    Requires :py:mod:`pandas`.

    Parameters
    ----------
    key : str / iterable / callable, optional
        Default is 'probability'.
    is_decoy : str / iterable / callable, optional
        Default is to check that "protein_name" starts with `'DECOY_'`.
    reverse : bool, optional
        Should be :py:const:`True` if higher score is better.
        Default is :py:const:`True` (because the default key is 'probability').

    *args, **kwargs : passed to :py:func:`auxiliary.filter` and/or :py:func:`DataFrame`.

    Returns
    -------
    out : pandas.DataFrame
    """
    import pandas as pd
    kwargs.setdefault('key', 'probability')
    kwargs.setdefault('reverse', True)
    if all(isinstance(arg, pd.DataFrame) for arg in args):
        if len(args) > 1:
            df = pd.concat(args)
        else:
            df = args[0]
    else:
        read_kw = {k: kwargs.pop(k) for k in ['iterative', 'read_schema', 'sep'] if k in kwargs}
        df = DataFrame(*args, **read_kw)
    if 'is_decoy' not in kwargs:
        if 'decoy_suffix' in kwargs:
            kwargs['is_decoy'] = df['protein_name'].str.endswith(kwargs['decoy_suffix'])
        else:
            kwargs['is_decoy'] = df['protein_name'].str.startswith(kwargs.get('decoy_prefix', 'DECOY_'))
    return aux.filter(df, **kwargs)