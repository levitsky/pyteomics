"""
pepxml - pepXML file reader
===========================

Summary
-------

`pepXML <http://tools.proteomecenter.org/wiki/index.php?title=Formats:pepXML>`_
was the first widely accepted format for proteomics search engines' output.
Even though it is to be replaced by a community standard
`mzIdentML <http://www.psidev.info/index.php?q=node/454>`_, it is still used
commonly.

This module provides minimalistic infrastructure for access to data stored in
pepXML files. The most important function is :py:func:`read`, which
reads peptide-spectum matches and related information and saves them into
human-readable dicts. This function relies on the terminology of the underlying
`lxml library <http://lxml.de/>`_.

Data access
-----------

  :py:class:`PepXML` - a class representing a single pepXML file.
  Other data access functions use this class internally.

  :py:func:`read` - iterate through peptide-spectrum matches in a pepXML
  file. Data for a single spectrum are converted to an easy-to-use dict.

  :py:func:`chain` - read multiple files at once.

  :py:func:`chain.from_iterable` - read multiple files at once, using an
  iterable of files.

  :py:func:`DataFrame` - read pepXML files into a :py:class:`pandas.DataFrame`.

  :py:func:`filter` - filter PSMs from a chain of pepXML files to a specific FDR
  using TDA.

  :py:func:`filter.chain` - chain a series of filters applied independently to
  several files.

  :py:func:`filter.chain.from_iterable` - chain a series of filters applied
  independently to an iterable of files.

  :py:func:`filter_df` - filter pepXML files and return a :py:class:`pandas.DataFrame`.

Miscellaneous
-------------

  :py:func:`fdr` - estimate the false discovery rate of a PSM set using the
  target-decoy approach.

  :py:func:`qvalues` - get an array of scores and local FDR values for a PSM
  set using the target-decoy approach.

  :py:func:`is_decoy` - determine whether a PSM is decoy or not.

  :py:func:`roc_curve` - get a receiver-operator curve (min peptideprophet
  probability is a sample vs. false discovery rate) of peptideprophet analysis.

Deprecated functions
--------------------

  :py:func:`iterfind` - iterate over elements in a pepXML file.
  You can just call the corresponding method of the :py:class:`PepXML`
  object.

  :py:func:`version_info` - get information about pepXML version and schema.
  You can just read the corresponding attribute of the :py:class:`PepXML`
  object.

Dependencies
------------

This module requires :py:mod:`lxml`.

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

from lxml import etree
from . import xml, auxiliary as aux

class PepXML(xml.XML):
    """Parser class for pepXML files."""
    file_format = 'pepXML'
    _root_element = 'msms_pipeline_analysis'
    _default_schema = xml._pepxml_schema_defaults
    _default_version = '1.15'
    _default_iter_tag = 'spectrum_query'
    _structures_to_flatten = {'search_score_summary', 'modification_info'}

    def _get_info_smart(self, element, **kw):
        """Extract the info in a smart way depending on the element type"""
        name = xml._local_name(element)
        kwargs = dict(kw)
        rec = kwargs.pop('recursive', None)
        if name == 'msms_pipeline_analysis':
            info = self._get_info(element,
                    recursive=(rec if rec is not None else False),
                    **kwargs)
        else:
            info = self._get_info(element,
                    recursive=(rec if rec is not None else True),
                    **kwargs)

        # attributes which contain unconverted values
        convert = {'float':  {'calc_neutral_pep_mass', 'massdiff'},
            'int': {'start_scan', 'end_scan', 'index'},
            'bool': {'is_rejected'},
            'floatarray': {'all_ntt_prob'}}
        def safe_float(s):
            try:
                return float(s)
            except ValueError:
                if s.startswith('+-0'): return 0
                return None

        converters = {'float': safe_float, 'int': int,
                'bool': lambda x: x.lower() in {'1', 'true'},
                'floatarray': lambda x: list(map(float, x[1:-1].split(',')))}
        for k, v in dict(info).items():
            for t, s in convert.items():
                if k in s:
                    del info[k]
                    info[k] = converters[t](v)
        for k in {'search_score', 'parameter'}:
            if k in info and isinstance(info[k], list) and all(
                    isinstance(x, dict) and len(x) == 1 for x in info[k]):
                scores = {}
                for score in info[k]:
                    name, value = score.popitem()
                    try:
                        scores[name] = float(value)
                    except ValueError:
                        scores[name] = value
                info[k] = scores
        if 'search_result' in info and len(info['search_result']) == 1:
            info.update(info['search_result'][0])
            del info['search_result']
        if 'protein' in info and 'peptide' in info:
            info['proteins'] = [{'protein': info.pop('protein'),
                'protein_descr': info.pop('protein_descr', None)}]
            for add_key in {'peptide_prev_aa', 'peptide_next_aa', 'protein_mw'}:
                if add_key in info:
                    info['proteins'][0][add_key] = info.pop(add_key)
            info['proteins'][0]['num_tol_term'] = info.pop('num_tol_term', 0)
            if 'alternative_protein' in info:
                info['proteins'].extend(info['alternative_protein'])
                del info['alternative_protein']
        if 'peptide' in info and not 'modified_peptide' in info:
            info['modified_peptide'] = info['peptide']
        if 'peptide' in info:
            info['modifications'] = info.pop('mod_aminoacid_mass', [])
            if 'mod_nterm_mass' in info:
                info['modifications'].insert(0, {'position': 0,
                    'mass': float(info.pop('mod_nterm_mass'))})
            if 'mod_cterm_mass' in info:
                info['modifications'].append({'position': 1 + len(info['peptide']),
                    'mass': float(info.pop('mod_cterm_mass'))})
        if 'modified_peptide' in info and info['modified_peptide'] == info.get(
                'peptide'):
            if not info.get('modifications'):
                info['modifications'] = []
            else:
                mp = info['modified_peptide']
                for mod in sorted(info['modifications'],
                        key=lambda m: m['position'],
                        reverse=True):
                    if mod['position'] not in {0, 1+len(info['peptide'])}:
                        p = mod['position']
                        mp = mp[:p] + '[{}]'.format(int(mod['mass'])) + mp[p:]
                info['modified_peptide'] = mp
        if 'search_hit' in info:
            info['search_hit'].sort(key=lambda x: x['hit_rank'])
        return info

def read(source, read_schema=True, iterative=True):
    """Parse `source` and iterate through peptide-spectrum matches.

    Parameters
    ----------
    source : str or file
        A path to a target pepXML file or the file object itself.

    read_schema : bool, optional
        If :py:const:`True`, attempt to extract information from the XML schema
        mentioned in the pepXML header (default). Otherwise, use default
        parameters. Disable this to avoid waiting on long network connections or
        if you don't like to get the related warnings.

    iterative : bool, optional
        Defines whether iterative parsing should be used. It helps reduce
        memory usage at almost the same parsing speed. Default is
        :py:const:`True`.

    Returns
    -------
    out : PepXML
       An iterator over dicts with PSM properties.
    """

    return PepXML(source, read_schema=read_schema, iterative=iterative)

def iterfind(source, path, **kwargs):
    """Parse `source` and yield info on elements with specified local
    name or by specified "XPath".

    .. note:: This function is provided for backward compatibility only.
        If you do multiple :py:func:`iterfind` calls on one file, you should
        create an :py:class:`PepXML` object and use its
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
    return PepXML(source, **kwargs).iterfind(path, **kwargs)

version_info = xml._make_version_info(PepXML)

def roc_curve(source):
    """Parse source and return a ROC curve for peptideprophet analysis.

    Parameters
    ----------
    source : str or file
        A path to a target pepXML file or the file object itself.

    Returns
    -------
    out : list
        A list of ROC points, sorted by ascending min prob.
    """

    parser = etree.XMLParser(remove_comments=True, ns_clean=True)
    tree = etree.parse(source, parser=parser)

    roc_curve = []
    for roc_element in tree.xpath(
        "/*[local-name()='msms_pipeline_analysis']"
        "/*[local-name()='analysis_summary and @analysis='peptideprophet']"
        "/*[local-name()='peptideprophet_summary']"
        "/*[local-name()='roc_data_point']"):

        roc_data_point = dict(roc_element.attrib)
        for key in roc_data_point:
            roc_data_point[key] = float(roc_data_point[key])
        roc_curve.append(roc_data_point)

    return sorted(roc_curve, key=lambda x: x['min_prob'])

chain = aux._make_chain(read, 'read')

def is_decoy(psm, prefix='DECOY_'):
    """Given a PSM dict, return :py:const:`True` if all protein names for
    the PSM start with ``prefix``, and :py:const:`False` otherwise. This
    function might not work for some pepXML flavours. Use the source to get the
    idea and suit it to your needs.

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
    return all(protein['protein'].startswith(prefix)
            for protein in  psm['search_hit'][0]['proteins'])

fdr = aux._make_fdr(is_decoy)
_key = lambda x: min(
    sh['search_score']['expect'] for sh in x['search_hit'])
qvalues = aux._make_qvalues(chain, is_decoy, _key)
filter = aux._make_filter(chain, is_decoy, _key, qvalues)
filter.chain = aux._make_chain(filter, 'filter', True)

def DataFrame(*args, **kwargs):
    """Read pepXML output files into a :py:class:`pandas.DataFrame`.

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
    sep = kwargs.pop('sep', None)
    with chain(*args, **kwargs) as f:
        for item in f:
            info = {}
            for k, v in item.items():
                if isinstance(v, (str, int, float)):
                    info[k] = v
            sh = item['search_hit'][0]
            proteins = sh.pop('proteins')
            prot_dict = {}
            for p in proteins:
                for k in p:
                    prot_dict[k] = []
            for p in proteins:
                for k, v in prot_dict.items():
                    v.append(p.get(k))
            if sep is None:
                info.update(prot_dict)
            else:
                for k, v in prot_dict.items():
                    info[k] = sep.join(str(val) if val is not None else '' for val in v)
            info.update(sh.pop('search_score'))
            mods = sh.pop('modifications', [])
            info['modifications'] = ','.join('{0[mass]:.3f}@{0[position]}'.format(x) for x in mods)
            for k, v in sh.items():
                if isinstance(v, (str, int, float)):
                    info[k] = v
            data.append(info)
    return pd.DataFrame(data)

def filter_df(*args, **kwargs):
    """Read pepXML files or DataFrames and return a :py:class:`DataFrame` with filtered PSMs.
    Positional arguments can be pepXML files or DataFrames.

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
            df['protein'].str.split(';').apply(lambda s: all(x.startswith('DECOY') for x in s)))
    else:
        kwargs.setdefault('is_decoy',
            df['protein'].apply(lambda s: all(x.startswith('DECOY') for x in s)))
    return aux.filter(df, **kwargs)