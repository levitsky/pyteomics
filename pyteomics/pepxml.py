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

Target-decoy approach
---------------------

  :py:func:`filter` - filter PSMs from a chain of pepXML files to a specific FDR
  using TDA.

  :py:func:`filter.chain` - chain a series of filters applied independently to
  several files.

  :py:func:`filter.chain.from_iterable` - chain a series of filters applied
  independently to an iterable of files.

  :py:func:`filter_df` - filter pepXML files and return a :py:class:`pandas.DataFrame`.

  :py:func:`fdr` - estimate the false discovery rate of a PSM set using the
  target-decoy approach.

  :py:func:`qvalues` - get an array of scores and local FDR values for a PSM
  set using the target-decoy approach.

  :py:func:`is_decoy` - determine whether a PSM is decoy or not.

Miscellaneous
-------------

  :py:func:`roc_curve` - get a receiver-operator curve (min PeptideProphet
  probability in a sample vs. false discovery rate) of PeptideProphet analysis.

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
from . import xml, auxiliary as aux, _schema_defaults


class PepXML(xml.MultiProcessingXML, xml.IndexSavingXML):
    """
    Parser class for pepXML files.
    """
    file_format = 'pepXML'
    _root_element = 'msms_pipeline_analysis'
    _default_schema = _schema_defaults._pepxml_schema_defaults
    _default_version = '1.15'
    _default_iter_tag = 'spectrum_query'
    _indexed_tags = {'spectrum_query'}
    _indexed_tag_keys = {'spectrum_query': 'spectrum'}
    _default_id_attr = 'spectrum'
    _structures_to_flatten = {'search_score_summary', 'modification_info'}
    # attributes which contain unconverted values
    _convert_items = {'float': {'calc_neutral_pep_mass', 'massdiff', 'probability', 'variable', 'static'},
                      'int': {'start_scan', 'end_scan', 'index', 'num_matched_peptides'},
                      'bool': {'is_rejected'},
                      'floatarray': {'all_ntt_prob'}}.items()

    def _get_info_smart(self, element, **kwargs):
        """Extract the info in a smart way depending on the element type"""
        try:
            name = kwargs.pop('ename')
        except KeyError:
            name = xml._local_name(element)
        rec = kwargs.pop('recursive', None)
        if name == self._root_element:
            info = self._get_info(element, ename=name, recursive=(rec if rec is not None else False), **kwargs)
        else:
            info = self._get_info(element, ename=name, recursive=(rec if rec is not None else True), **kwargs)

        def safe_float(s):
            try:
                return float(s)
            except ValueError:
                if s.startswith('+-0'):
                    return 0
                return s

        converters = {'float': safe_float, 'int': int,
                      'bool': lambda x: x.lower() in {'1', 'true'},
                      'floatarray': lambda x: list(map(float, x[1:-1].split(',')))}
        for k, v in dict(info).items():
            for t, s in self._convert_items:
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
            info['proteins'] = [{'protein': info.pop('protein'), 'protein_descr': info.pop('protein_descr', None)}]
            for add_key in {'peptide_prev_aa', 'peptide_next_aa', 'protein_mw'}:
                if add_key in info:
                    info['proteins'][0][add_key] = info.pop(add_key)
            info['proteins'][0]['num_tol_term'] = info.pop('num_tol_term', 0)
            if 'alternative_protein' in info:
                info['proteins'].extend(info['alternative_protein'])
                del info['alternative_protein']
        if 'peptide' in info and 'modified_peptide' not in info:
            info['modified_peptide'] = info['peptide']
        if 'peptide' in info:
            info['modifications'] = info.pop('mod_aminoacid_mass', [])
            if 'mod_nterm_mass' in info:
                info['modifications'].insert(0, {'position': 0, 'mass': float(info.pop('mod_nterm_mass'))})
            if 'mod_cterm_mass' in info:
                info['modifications'].append({'position': 1 + len(info['peptide']), 'mass': float(info.pop('mod_cterm_mass'))})
        if 'modified_peptide' in info and info['modified_peptide'] == info.get(
                'peptide'):
            if not info.get('modifications'):
                info['modifications'] = []
            else:
                mp = info['modified_peptide']
                for mod in sorted(info['modifications'], key=lambda m: m['position'], reverse=True):
                    if mod['position'] not in {0, 1+len(info['peptide'])}:
                        p = mod['position']
                        mp = mp[:p] + '[{}]'.format(int(mod['mass'])) + mp[p:]
                info['modified_peptide'] = mp
        if 'search_hit' in info:
            info['search_hit'].sort(key=lambda x: x['hit_rank'])
        return info

    def search_hits(self):
        """
        Iterate over search hits rather than spectrum queries.
        """
        for sq in self:
            sh = sq.pop('search_hit', [])
            for item in sh:
                item.update(sq)
                yield item


def read(*args, **kwargs):
    """Parse `source` and iterate through peptide-spectrum matches.

    Parameters
    ----------
    source : str or file
        A path to a target pepXML file or the file object itself.

    read_schema : bool, optional
        If :py:const:`True`, attempt to extract information from the XML schema
        mentioned in the pepXML header. Otherwise, use default parameters.
        Not recommended without Internet connection or
        if you don't like to get the related warnings.

    iterative : bool, optional
        Defines whether iterative parsing should be used. It helps reduce
        memory usage at almost the same parsing speed. Default is
        :py:const:`True`.

    use_index : bool, optional
        Defines whether an index of byte offsets needs to be created for
        elements listed in `indexed_tags`.
        This is useful for random access to spectum queries.
        Default is :py:const:`True`.

    indexed_tags : container of bytes, optional
        If `use_index` is :py:const:`True`, elements listed in this parameter
        will be indexed. Empty set by default.

    Returns
    -------
    out : PepXML
       An iterator over dicts with PSM properties.
    """

    return PepXML(*args, **kwargs)


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

    recursive : bool, keyword only, optional
        If :py:const:`False`, subelements will not be processed when
        extracting info from elements. Default is :py:const:`True`.

    iterative : bool, keyword only, optional
        Specifies whether iterative XML parsing should be used. Iterative
        parsing significantly reduces memory usage and may be just a little
        slower. When `retrieve_refs` is :py:const:`True`, however, it is
        highly recommended to disable iterative parsing if possible.
        Default value is :py:const:`True`.

    read_schema : bool, keyword only, optional
        If :py:const:`True`, attempt to extract information from the XML schema
        mentioned in the mzIdentML header. Otherwise, use default parameters.
        Not recommended without Internet connection or
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
        A list of ROC points.
    """

    parser = etree.XMLParser(remove_comments=True, ns_clean=True)
    tree = etree.parse(source, parser=parser)

    roc_curve = []
    for roc_error_data in tree.xpath(
        "/*[local-name()='msms_pipeline_analysis'] \
        //*[local-name()='analysis_summary' and @analysis='peptideprophet'] \
        //*[local-name()='peptideprophet_summary'] \
        //*[local-name()='roc_error_data']"):
        for element in roc_error_data.xpath("*[local-name()='roc_data_point' or local-name()='error_point']"):
            data_point = dict(element.attrib)
            for key in data_point:
                data_point[key] = float(data_point[key])
            data_point["charge"] = roc_error_data.attrib["charge"]
            data_point["tag"] = etree.QName(element).localname
            roc_curve.append(data_point)

    return roc_curve


# chain = aux._make_chain(read, 'read')
chain = aux.ChainBase._make_chain(read)


def _is_decoy_prefix(psm, prefix='DECOY_'):
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
               for protein in psm['search_hit'][0]['proteins'])


def _is_decoy_suffix(psm, suffix='_DECOY'):
    return all(protein['protein'].endswith(suffix)
               for protein in psm['search_hit'][0]['proteins'])


is_decoy = _is_decoy_prefix
fdr = aux._make_fdr(_is_decoy_prefix, _is_decoy_suffix)
_key = lambda x: min(sh['search_score']['expect'] for sh in x['search_hit'])
qvalues = aux._make_qvalues(chain, _is_decoy_prefix, _is_decoy_suffix, _key)
filter = aux._make_filter(chain, _is_decoy_prefix, _is_decoy_suffix, _key, qvalues)
filter.chain = aux._make_chain(filter, 'filter', True)


def DataFrame(*args, **kwargs):
    """Read pepXML output files into a :py:class:`pandas.DataFrame`.

    Requires :py:mod:`pandas`.

    Parameters
    ----------
    *args
        pepXML file names or objects. Passed to :py:func:`chain`.

    **kwargs
        Passed to :py:func:`chain`.

    by : str, keyword only, optional
        Can be :py:const:`"spectrum_query"` (default) or :py:const:`"search_hit"`.
        One row in the resulting dataframe corresponds to one element of the given type.
        If :py:const:`"spectrum_query"` is set, only the top search hit is shown in the dataframe.

    sep : str or None, keyword only, optional
        Some values related to PSMs (such as protein information) are variable-length
        lists. If `sep` is a :py:class:`str`, they will be packed into single string using
        this delimiter. If `sep` is :py:const:`None`, they are kept as lists. Default is
        :py:const:`None`.

    recursive : bool, keyword only, optional
        If :py:const:`False`, subelements will not be processed when
        extracting info from elements. Default is :py:const:`True`.

    iterative : bool, keyword only, optional
        Specifies whether iterative XML parsing should be used. Iterative
        parsing significantly reduces memory usage and may be just a little
        slower. When `retrieve_refs` is :py:const:`True`, however, it is
        highly recommended to disable iterative parsing if possible.
        Default value is :py:const:`True`.

    read_schema : bool, keyword only, optional
        If :py:const:`True`, attempt to extract information from the XML schema
        mentioned in the mzIdentML header. Otherwise, use default parameters.
        Not recommended without Internet connection or
        if you don't like to get the related warnings.

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
    kwargs.setdefault('use_index', False)
    by = kwargs.pop('by', 'spectrum_query')

    def search_hit_info(sh):
        info = {}
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
        formatted_mods = ['{0[mass]:.3f}@{0[position]}'.format(x) for x in mods]
        if sep is not None:
            info['modifications'] = sep.join(formatted_mods)
        else:
            info['modifications'] = formatted_mods
        for k, v in sh.items():
            if isinstance(v, (str, int, float)):
                info[k] = v
        if 'analysis_result' in sh:
            for ar in sh['analysis_result']:
                if ar['analysis'] == 'peptideprophet':
                    try:
                        info.update(ar['peptideprophet_result']['parameter'])
                    except KeyError:
                        pass
                    info['peptideprophet_probability'] = ar['peptideprophet_result']['probability']
                    info['peptideprophet_ntt_prob'] = ar['peptideprophet_result']['all_ntt_prob']
                elif ar['analysis'] == 'interprophet':
                    info.update(ar['interprophet_result']['parameter'])
                    info['interprophet_probability'] = ar['interprophet_result']['probability']
                    info['interprophet_ntt_prob'] = ar['interprophet_result']['all_ntt_prob']
        return info

    def iter_spectrum_query():
        with chain(*args, **kwargs) as f:
            for item in f:
                info = {}
                for k, v in item.items():
                    if isinstance(v, (str, int, float)):
                        info[k] = v
                if 'search_hit' in item:
                    sh = item['search_hit'][0]
                    info.update(search_hit_info(sh))
                yield info

    def iter_search_hit():
        for source in args:
            with PepXML(source, **kwargs) as f:
                for sh in f.search_hits():
                    yield search_hit_info(sh)

    items = {'spectrum_query': iter_spectrum_query, 'search_hit': iter_search_hit}[by]
    return pd.DataFrame(items(), **pd_kwargs)


def filter_df(*args, **kwargs):
    """Read pepXML files or DataFrames and return a :py:class:`DataFrame` with filtered PSMs.
    Positional arguments can be pepXML files or DataFrames. Keyword parameter `fdr` is also required.
    Other parameters are optional.

    Requires :py:mod:`pandas`.

    Parameters
    ----------
    positional args
        pepXML file names, file objects, or DataFrames. Passed to :py:func:`DataFrame`.
    fdr : float, keyword only, 0 <= fdr <= 1
        Desired FDR level.
    key : str / iterable / callable, keyword only, optional
        PSM score. Default is 'expect'.
    is_decoy : str / iterable / callable, keyword only, optional
        Default is to check if all strings in the "protein" column start with `'DECOY_'`.
    sep : str or None, keyword only, optional
        Some values related to PSMs (such as protein information) are variable-length
        lists. If `sep` is a :py:class:`str`, they will be packed into single string using
        this delimiter. If `sep` is :py:const:`None`, they are kept as lists. Default is
        :py:const:`None`.
    reverse : bool, keyword only, optional
        If :py:const:`True`, then PSMs are sorted in descending order,
        i.e. the value of the key function is higher for better PSMs.
        Default is :py:const:`False`.
    decoy_prefix : str, optional
        If the default `is_decoy` function works for you, this parameter specifies which
        protein name prefix to use to detect decoy matches. If you provide your own
        `is_decoy`, or if you specify `decoy_suffix`, this parameter has no effect.
        Default is `"DECOY_"`.
    decoy_suffix : str, optional
        If the default `is_decoy` function works for you, this parameter specifies which
        protein name suffix to use to detect decoy matches. If you provide your own
        `is_decoy`, this parameter has no effect. Mutually exclusive with `decoy_prefix`.
    remove_decoy : bool, keyword only, optional
        Defines whether decoy matches should be removed from the output.
        Default is :py:const:`True`.

        .. note:: If set to :py:const:`False`, then by default the decoy
           PSMs will be taken into account when estimating FDR. Refer to the
           documentation of :py:func:`fdr` for math; basically, if
           `remove_decoy` is :py:const:`True`, then formula 1 is used
           to control output FDR, otherwise it's formula 2. This can be
           changed by overriding the `formula` argument.

    formula : int, keyword only, optional
        Can be either 1 or 2, defines which formula should be used for FDR
        estimation. Default is 1 if `remove_decoy` is :py:const:`True`,
        else 2 (see :py:func:`fdr` for definitions).
    ratio : float, keyword only, optional
        The size ratio between the decoy and target databases. Default is
        1. In theory, the "size" of the database is the number of
        theoretical peptides eligible for assignment to spectra that are
        produced by *in silico* cleavage of that database.
    correction : int or float, keyword only, optional
        Possible values are 0, 1 and 2, or floating point numbers between 0 and 1.

        0 (default): no correction;

        1: enable "+1" correction. This accounts for the probability that a false
        positive scores better than the first excluded decoy PSM;

        2: this also corrects that probability for finite size of the sample,
        so the correction will be slightly less than "+1".

        If a floating point number
        is given, then instead of the expectation value for the number of false PSMs,
        the confidence value is used. The value of `correction` is then interpreted as
        desired confidence level. E.g., if correction=0.95, then the calculated q-values
        do not exceed the "real" q-values with 95% probability.

        See `this paper <http://dx.doi.org/10.1021/acs.jproteome.6b00144>`_ for further explanation.

    pep : callable / array-like / iterable / str, keyword only, optional
        If callable, a function used to determine the posterior error probability (PEP).
        Should accept exactly one argument (PSM) and return a float.
        If array-like, should contain float values for all given PSMs.
        If string, it is used as a field name (PSMs must be in a record array
        or a :py:class:`DataFrame`).

        .. note:: If this parameter is given, then PEP values will be used to calculate
           q-values. Otherwise, decoy PSMs will be used instead. This option conflicts with:
           `is_decoy`, `remove_decoy`, `formula`, `ratio`, `correction`.
           `key` can still be provided. Without `key`, PSMs will be sorted by PEP.

    q_label : str, optional
        Field name for q-value in the output. Default is ``'q'``.

    score_label : str, optional
        Field name for score in the output. Default is ``'score'``.

    decoy_label : str, optional
        Field name for the decoy flag in the output. Default is ``'is decoy'``.

    pep_label : str, optional
        Field name for PEP in the output. Default is ``'PEP'``.

    Returns
    -------
    out : pandas.DataFrame
    """
    import pandas as pd
    sep = kwargs.get('sep')
    kwargs.setdefault('key', 'expect')
    if all(isinstance(arg, pd.DataFrame) for arg in args):
        if len(args) > 1:
            df = pd.concat(args)
        else:
            df = args[0]
    else:
        read_kw = {k: kwargs.pop(k) for k in ['iterative', 'read_schema', 'sep', 'pd_kwargs'] if k in kwargs}
        df = DataFrame(*args, **read_kw)
    if 'is_decoy' not in kwargs:
        if sep is not None:
            if 'decoy_suffix' in kwargs:
                kwargs['is_decoy'] = df['protein'].str.split(';').apply(
                    lambda s: all(x.endswith(kwargs['decoy_suffix']) for x in s))
            else:
                kwargs['is_decoy'] = df['protein'].str.split(';').apply(
                    lambda s: all(x.startswith(kwargs.get('decoy_prefix', 'DECOY_')) for x in s))
        else:
            if 'decoy_suffix' in kwargs:
                kwargs['is_decoy'] = df['protein'].apply(
                    lambda s: all(x.endswith(kwargs['decoy_suffix']) for x in s))
            else:
                kwargs['is_decoy'] = df['protein'].apply(
                    lambda s: all(x.startswith(kwargs.get('decoy_prefix', 'DECOY_')) for x in s))
    return aux.filter(df, **kwargs)
