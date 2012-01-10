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
pep.XML files. The most important function is :py:func:`iter_psm`, which 
reads peptide-spectum matches and related information and saves them into 
human-readable dicts. The rest of data can be obtained via :py:func:`get_node` 
function. This functions relies on the terminology of the underlying 
`lxml library <http://lxml.de/>`_.

Data access
-----------

  :py:func:`iter_psm` - iterate through peptide-spectrum matches in a pep.XML 
  file. Data from a single PSM are converted to a human-readable dict. 

  :py:func:`get_node` - get arbitrary nodes of pep.XML file by their xpath.

  :py:func:`roc_curve` - get a receiver-operator curve (min peptideprophet
  probability is a sample vs. false discovery rate) of peptideprophet analysis.

-------------------------------------------------------------------------------
"""
# Licensed under the MIT license:
# http://www.opensource.org/licenses/mit-license.php 

from lxml import etree

# A list of the spectrum_query attributes which contain float values.
float_keys = [
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
    ]

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

    Returns
    -------
    out : str
        A modified XPath.
    """
    return '/'.join(
        [(ns_key + ':' + node if (node and node.count(':') == 0) else node)
         for node in xpath.split('/')])

def get_node(source, xpath, namespaces={'d':xmlns}):
    """Retrieves arbitrary nodes from a pepxml file by their xpath.
    Each node in the xpath is assigned to the default namespace 
    'http://regis-web.systemsbiology.net/pepXML' unless specified else.

    Parameters
    ----------
    source : str or file
        A path or an URL to a target pepXML file or the file object itself.
    xpath : str
        An XPath to target nodes. 
    namespaces : dict, optional
        A dictionary of namespaces. The default namespace key is 'd'.
    
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

    xpath_w_namespace = _insert_default_ns(xpath)

    return tree.xpath(xpath_w_namespace, namespaces=namespaces)
    
def _psm_from_query(query, namespaces={'d':xmlns}):
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
    # ... and the best hit from peptide database.
    search_hit_elements = query.xpath(
        'd:search_result/d:search_hit[@hit_rank=\'1\']',
        namespaces=namespaces)
    if not search_hit_elements:
        return {}
    psm.update(search_hit_elements[0].attrib)

    # Use non-modified sequence of a peptide as a modified one if
    # there are no modifications.
    psm['modified_peptide'] = psm['peptide']

    # Convert str values to float.
    for key in float_keys:
        if key in psm:
            psm[key] = float(psm[key])    

    # Store alternative proteins as a list of dictionaries.
    proteins = []
    proteins.append(
        {"protein":         psm.pop("protein"),
         "protein_descr":   psm.pop("protein_descr", ""),
         "num_tol_term":    float(psm.pop("num_tol_term")),
         "peptide_prev_aa": psm.pop("peptide_prev_aa"),
         "peptide_next_aa": psm.pop("peptide_next_aa")})

    # Store a list of modifications.
    modifications = []
    for subelement in query.xpath('d:search_result/d:search_hit/'
                                  'd:search_score', 
                                  namespaces=namespaces):
        psm[subelement.attrib['name']] = float(subelement.attrib['value'])
    for subelement in query.xpath('d:search_result/d:search_hit/'
                                  'd:analysis_result/d:peptideprophet_result',
                                  namespaces=namespaces):
        psm['peptideprophet'] = float(subelement.attrib['probability'])
    for subelement in query.xpath('d:search_result/d:search_hit/'
                                  'd:alternative_protein',
                                  namespaces=namespaces):
        proteins.append(subelement.attrib)
    for subelement in query.xpath('d:search_result/d:search_hit/'
                                  'd:modification_info',
                                  namespaces=namespaces):
        psm['modified_peptide'] = subelement.attrib.get('modified_peptide',
                                                        '')
        if 'mod_nterm_mass' in subelement.attrib:
            modifications.append(
                {'position' : 0,
                 'mass': float(subelement.attrib['mod_nterm_mass'])})
        if 'mod_cterm_mass' in subelement.attrib:
            modifications.append(
                {'position' : _peptide_length(psm) + 1,
                 'mass': float(subelement.attrib['mod_cterm_mass'])})
        for mod_element in subelement.xpath('d:mod_aminoacid_mass',
                                            namespaces=namespaces):
            modification = dict(mod_element.attrib)
            modification['position'] = int(modification['position'])
            modification['mass'] = float(modification['mass'])
            modifications.append(modification)                

    psm["proteins"] = proteins
    psm["modifications"] = modifications

    return psm
    
def iter_psm(source):
    """Parse ``source`` and iterate through peptide-spectrum matches.

    Parameters
    ----------
    source : str or file
        A path or an URL to a target pepXML file or the file object itself.

    Returns
    -------
    out : iterator
       An iterator over the dicts with PSM properties.
    """

    parser = etree.XMLParser(remove_comments=True, ns_clean=True) 
    tree = etree.parse(source, parser=parser)

    for spectrum_query in tree.xpath(
        '/d:msms_pipeline_analysis/d:msms_run_summary/d:spectrum_query',
        namespaces = {'d': xmlns}):
        
        yield _psm_from_query(spectrum_query)

def roc_curve(source):
    """Parse source and return a ROC curve for peptideprophet analysis.

    Parameters
    ----------
    source : str or file
        A path or an URL to a target pepXML file or the file object itself.

    Returns
    -------
    out : list
        A list of ROC points, sorted by ascending min prob.
    """

    parser = etree.XMLParser(remove_comments=True, ns_clean=True) 
    tree = etree.parse(source, parser=parser)

    roc_curve = []
    for roc_element in tree.xpath(
        '/d:msms_pipeline_analysis'
        '/d:analysis_summary[@analysis=\'peptideprophet\']'
        '/d:peptideprophet_summary/d:roc_data_point',
        namespaces = {'d': xmlns}):
        
        roc_data_point = dict(roc_element.attrib)
        for key in roc_data_point:
            roc_data_point[key] = float(roc_data_point[key])
        roc_curve.append(roc_data_point)

    return sorted(roc_curve, key=lambda x: x['min_prob'])
    
