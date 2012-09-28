"""
pepxml - reader for peptide-spectrum matches in pep.XML format
==============================================================

Summary
-------

`pep.XML <http://tools.proteomecenter.org/wiki/index.php?title=Formats:pepXML>`_
was the first widely accepted format for proteomics search engines' output. 
Even though it is to be replaced by a community standard
`mzIdentML <http://www.psidev.info/index.php?q=node/454>`_, it is still used
commonly.

This module provides minimalistic infrastructure for access to data stored in
pep.XML files. The most important function is :py:func:`read`, which 
reads peptide-spectum matches and related information and saves them into 
human-readable dicts. The rest of data can be obtained via :py:func:`get_node` 
function. This function relies on the terminology of the underlying 
`lxml library <http://lxml.de/>`_.

Data access
-----------

  :py:func:`read` - iterate through peptide-spectrum matches in a pep.XML 
  file. Data from a single PSM are converted to a human-readable dict. 

  :py:func:`get_node` - get arbitrary nodes of pep.XML file by their xpath.

  :py:func:`roc_curve` - get a receiver-operator curve (min peptideprophet
  probability is a sample vs. false discovery rate) of peptideprophet analysis.

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
from . import mzid
from . import auxiliary as aux
try:
    reduce # Python 2.7
except NameError: # Python 3.x
    from functools import reduce

# spectrum_query attributes which contain float values.
float_keys = {
    'calc_neutral_pep_mass',
    'precursor_neutral_mass',
    'massdiff',
    'start_scan',
    'end_scan',
    'assumed_charge',
    'num_tot_proteins',
    'num_missed_cleavages',
    'retention_time_sec',
    'tot_num_ions',
    'num_matched_ions',
    'index',
    'hit_rank',
    'is_rejected'
    }

# Default namespace of pepXML.
xmlns = 'http://regis-web.systemsbiology.net/pepXML'

def _peptide_length(psm):
    return len(psm['peptide'])

def _insert_default_ns(xpath, ns_key = 'd'):
    """Inserts the key for the default namespace before each node name.
    Does not modify a nodename if it already has a namespace.
    
    Parameters
    ----------
    xpath : str
        An original XPath
    ns_key : str
        A key for the default namespace.
        If empty string, `xpath` will be returned unchanged.

    Returns
    -------
    out : str
        A modified XPath.
    """
    if not ns_key: return xpath
    return '/'.join(
        [(ns_key + ':' + node if (node and node.count(':') == 0) else node)
         for node in xpath.split('/')])

def get_node(source, xpath, namespaces={'d':xmlns}):
    """Retrieves arbitrary nodes from a pepXML file by their xpath.
    Each node in the xpath is assigned to the default namespace 
    'http://regis-web.systemsbiology.net/pepXML' unless specified else.

    Parameters
    ----------
    source : str or file
        A path to a target pepXML file or the file object itself.
    xpath : str
        An XPath to target nodes. 
    namespaces : dict, optional
        A dictionary of namespaces. The default namespace key is 'd'.
        If the XML document does not have a specified namespace,
        supply an empty dictionary.
    
    Returns
    -------
    out : list of lxml.Element 
        List of target nodes.

    Examples
    --------
    >> get_node('/msms_pipeline_analysis/msms_run_summary[0]'
                '/search_summary/aminoacid_modification')

    """
    parser = etree.XMLParser(remove_comments=True, ns_clean=True)
    tree = etree.parse(source, parser=parser)

    if not namespaces: ns_key = ''
    else: ns_key = 'd'

    xpath_w_namespace = _insert_default_ns(xpath, ns_key)
    kwargs = {}
    if namespaces: kwargs['namespaces'] = namespaces
       
    return tree.xpath(xpath_w_namespace, **kwargs)

    
def _psm_from_query(query):
    """Analyze a spectrum query Element object and generate a dictionary with 
    its properties.

    Parameters
    ----------
    element : lxml.Element
        A parent element with a spectrum query.

    Returns
    -------
    out : dict
    """

    # A psm stores the properties of the spectrum query...
    psm = dict(query.attrib)
    # ... and the info about all search hits (in a list).
    search_hit_elements = query.xpath(
        "*[local-name()='search_result']"
        "/*[local-name()='search_hit']")
    if not search_hit_elements:
        return {}

    # Convert str values to float.
    for key in float_keys:
            if key in psm:
                psm[key] = float(psm[key])

    psm['search_hits'] = []
    for search_hit in search_hit_elements:
        # form a dictionary with search hit info, then add it to the list
        search_hit_info = {}
        search_hit_info.update(search_hit.attrib)

        # Use non-modified sequence of a peptide as a modified one if
        # there are no modifications.
        search_hit_info['modified_peptide'] = search_hit_info['peptide']

        # Convert str values to float in search hit.
        for key in float_keys:
            if key in search_hit_info:
                search_hit_info[key] = float(search_hit_info[key])    

        # Store alternative proteins as a list of dictionaries.
        proteins = []
        proteins.append(
            {"protein":         search_hit_info.pop("protein"),
             "protein_descr":   search_hit_info.pop("protein_descr", ""),
             "num_tol_term":    float(search_hit_info.pop("num_tol_term"))
                                if "num_tol_term" in search_hit_info else None,
             "peptide_prev_aa": search_hit_info.pop("peptide_prev_aa"),
             "peptide_next_aa": search_hit_info.pop("peptide_next_aa")})

        for subelement in search_hit.xpath("*[local-name()='search_score']"):
            search_hit_info[subelement.attrib['name']] = float(
                    subelement.attrib['value'])
        for subelement in search_hit.xpath("*[local-name()='analysis_result']"
            "/*[local-name()='peptideprophet_result']"):
            search_hit_info['peptideprophet'] = float(
                    subelement.attrib['probability'])
        for subelement in search_hit.xpath(
                "*[local-name()='alternative_protein']"):
            proteins.append(subelement.attrib)

        # Store a list of modifications.
        modifications = []
        for subelement in search_hit.xpath(
                "*[local-name()='modification_info']"):
            search_hit_info['modified_peptide'] = subelement.attrib.get(
                'modified_peptide', '')
            if 'mod_nterm_mass' in subelement.attrib:
                modifications.append(
                    {'position' : 0,
                     'mass': float(subelement.attrib['mod_nterm_mass'])})
            if 'mod_cterm_mass' in subelement.attrib:
                modifications.append(
                    {'position' : _peptide_length(search_hit_info) + 1,
                     'mass': float(subelement.attrib['mod_cterm_mass'])})
            for mod_element in subelement.xpath(
                    "*[local-name()='mod_aminoacid_mass']"):
                modification = dict(mod_element.attrib)
                modification['position'] = int(modification['position'])
                modification['mass'] = float(modification['mass'])
                modifications.append(modification)
        
        if modifications and ((not search_hit_info['modified_peptide']) or 
                search_hit_info['modified_peptide'] == 
                search_hit_info['peptide']):
            seq = list(search_hit_info['peptide'])
            splitseq = []
            indices = sorted([x['position'] for x in modifications])
            reduce(lambda x, y: splitseq.append(seq[x:y]) or y,
                    indices + [len(seq)], 0)
            mods = ['[%d]' % int(round(x['mass'])) for x in sorted(
                modifications, key=lambda y: y['position'])]
            modseq = ''.join(splitseq[0])
            for i in range(len(mods)):
                modseq += mods[i] + ''.join(splitseq[i+1])
            search_hit_info['modified_peptide'] = modseq

        search_hit_info["proteins"] = proteins
        search_hit_info["modifications"] = modifications

        psm['search_hits'].append(search_hit_info)

    psm['search_hits'].sort(key = lambda x: x['hit_rank'])
    return psm
    
def read(source):
    """Parse ``source`` and iterate through peptide-spectrum matches.

    Parameters
    ----------
    source : str or file
        A path to a target pepXML file or the file object itself.

    Returns
    -------
    out : iterator
       An iterator over the dicts with PSM properties.
    """

    for _, tag in etree.iterparse(source):
        if tag.tag.endswith('spectrum_query'):
            yield _psm_from_query(tag)
            tag.clear()

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
    
_version_info_env = {'format': 'pepXML', 'element': 'msms_pipeline_analysis'}
version_info = aux._make_version_info(_version_info_env)

_schema_defaults = {'ints': 
    {('distribution_point', 'obs_5_distr'),
     ('distribution_point', 'obs_2_distr'),
     ('enzymatic_search_constraint', 'max_num_internal_cleavages'),
     ('asapratio_lc_heavypeak', 'right_valley'),
     ('libra_summary', 'output_type'),
     ('distribution_point', 'obs_7_distr'),
     ('data_filter', 'number'),
     ('search_hit', 'num_tol_term'),
     ('asapratio_lc_lightpeak', 'right_valley'),
     ('specificity', 'min_spacing'),
     ('database_refresh_timestamp', 'min_num_enz_term'),
     ('enzymatic_search_constraint', 'min_number_termini'),
     ('distribution_point', 'obs_3_distr'),
     ('search_database', 'size_in_db_entries'),
     ('alternative_protein', 'num_tol_term'),
     ('search_hit', 'tot_num_ions'),
     ('mixture_model', 'precursor_ion_charge'),
     ('search_hit', 'num_matched_ions'),
     ('dataset_derivation', 'generation_no'),
     ('spectrum_query', 'assumed_charge'),
     ('distribution_point', 'obs_4_distr'),
     ('libra_summary', 'normalization'),
     ('asapratio_lc_heavypeak', 'left_valley'),
     ('distribution_point', 'obs_6_distr'),
     ('asapratio_lc_lightpeak', 'left_valley'),
     ('search_database', 'size_of_residues'),
     ('asapratio_peptide_data', 'cidIndex'),
     ('mixture_model', 'num_iterations'),
     ('mod_aminoacid_mass', 'position'),
     ('asapratio_summary', 'area_flag'),
     ('mixture_model', 'tot_num_spectra'),
     ('distribution_point', 'obs_1_distr'),
     ('asapratio_contribution', 'charge'),
     ('libra_summary', 'centroiding_preference')},
    'floats':
    {('asapratio_contribution', 'error'),
     ('asapratio_lc_heavypeak', 'area_error'),
     ('modification_info', 'mod_nterm_mass'),
     ('distribution_point', 'model_4_neg_distr'),
     ('distribution_point', 'model_5_pos_distr'),
     ('spectrum_query', 'precursor_neutral_mass'),
     ('asapratio_lc_heavypeak', 'time_width'),
     ('xpressratio_summary', 'masstol'),
     ('affected_channel', 'correction'),
     ('distribution_point', 'model_7_neg_distr'),
     ('error_point', 'error'),
     ('intensity', 'target_mass'),
     ('roc_data_point', 'sensitivity'),
     ('distribution_point', 'model_4_pos_distr'),
     ('distribution_point', 'model_2_neg_distr'),
     ('distribution_point', 'model_3_pos_distr'),
     ('mixture_model', 'prior_probability'),
     ('roc_data_point', 'error'),
     ('intensity', 'normalized'),
     ('modification_info', 'mod_cterm_mass'),
     ('asapratio_lc_lightpeak', 'area_error'),
     ('distribution_point', 'fvalue'),
     ('distribution_point', 'model_1_neg_distr'),
     ('peptideprophet_summary', 'min_prob'),
     ('asapratio_result', 'mean'),
     ('point', 'pos_dens'),
     ('fragment_masses', 'mz'),
     ('mod_aminoacid_mass', 'mass'),
     ('distribution_point', 'model_6_neg_distr'),
     ('asapratio_lc_lightpeak', 'time_width'),
     ('asapratio_result', 'heavy2light_error'),
     ('peptideprophet_result', 'probability'),
     ('error_point', 'min_prob'),
     ('peptideprophet_summary', 'est_tot_num_correct'),
     ('roc_data_point', 'min_prob'),
     ('asapratio_result', 'heavy2light_mean'),
     ('distribution_point', 'model_5_neg_distr'),
     ('mixturemodel', 'neg_bandwidth'),
     ('asapratio_result', 'error'),
     ('xpressratio_result', 'light_mass'),
     ('point', 'neg_dens'),
     ('asapratio_lc_lightpeak', 'area'),
     ('distribution_point', 'model_1_pos_distr'),
     ('xpressratio_result', 'mass_tol'),
     ('mixturemodel', 'pos_bandwidth'),
     ('xpressratio_result', 'light_area'),
     ('asapratio_peptide_data', 'heavy_mass'),
     ('distribution_point', 'model_2_pos_distr'),
     ('search_hit', 'calc_neutral_pep_mass'),
     ('intensity', 'absolute'),
     ('asapratio_peptide_data', 'light_mass'),
     ('distribution_point', 'model_3_neg_distr'),
     ('aminoacid_modification', 'mass'),
     ('asapratio_lc_heavypeak', 'time'),
     ('asapratio_lc_lightpeak', 'time'),
     ('asapratio_lc_lightpeak', 'background'),
     ('mixture_model', 'est_tot_correct'),
     ('point', 'value'),
     ('asapratio_lc_heavypeak', 'background'),
     ('terminal_modification', 'mass'),
     ('fragment_masses', 'offset'),
     ('xpressratio_result', 'heavy_mass'),
     ('search_hit', 'protein_mw'),
     ('libra_summary', 'mass_tolerance'),
     ('spectrum_query', 'retention_time_sec'),
     ('distribution_point', 'model_7_pos_distr'),
     ('asapratio_lc_heavypeak', 'area'),
     ('alternative_protein', 'protein_mw'),
     ('asapratio_contribution', 'ratio'),
     ('xpressratio_result', 'heavy_area'),
     ('distribution_point', 'model_6_pos_distr')},
    'bools':
    {('sample_enzyme', 'independent'),
     ('intensity', 'reject'),
     ('libra_result', 'is_rejected')},
    'intlists': set(),
    'floatlists': set(),
    'charlists': set(),
    'lists': {'point', 'aminoacid_modification', 'msms_run_summary',
            'mixturemodel', 'search_hit', 'mixturemodel_distribution',
            'sequence_search_constraint', 'specificity', 'alternative_protein',
            'analysis_result', 'data_filter', 'fragment_masses', 'error_point',
            'parameter', 'spectrum_query', 'search_result', 'affected_channel',
            'analysis_summary', 'roc_data_point', 'distribution_point',
            'search_summary', 'mod_aminoacid_mass', 'search_score', 'intensity',
            'analysis_timestamp', 'mixture_model', 'terminal_modification',
            'contributing_channel', 'inputfile'}}

_schema_env = {'format': 'pepXML', 'version_info': version_info,
        'default_version': '1.15', 'defaults': _schema_defaults}
_schema_info = aux._make_schema_info(_schema_env)
