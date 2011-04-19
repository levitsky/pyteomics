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

def peptide_length(psm):
    return len(psm['peptide'])

def get_node(source, xpath, namespaces={'d':xmlns}):
    """Retrieves arbitrary nodes from a pepxml file by their xpath.

    Arguments: 

    source -- any of the following:    
              - a file name/path
              - a file object
              - a file-like object
              - a URL using the HTTP or FTP protocol

    xpath -- a string with XPath to required objects. Usually, pepXML
    has a default namespace
    'http://regis-web.systemsbiology.net/pepXML' which should be
    prepended to every nodename.

    namespaces -- a dictionary of namespaces. The default pepxml
    namespace has a key 'd'.
    
    Example:
    pepxml.get_node('/d:msms_pipeline_analysis/d:msms_run_summary[0]'
                    '/d:search_summary/d:aminoacid_modification')

    Returns:
    a list of lxml.Element objects.
    """

    parser = etree.XMLParser(remove_comments=True, ns_clean=True)
    tree = etree.parse(source, parser=parser)
    
    return tree.xpath(xpath, namespaces=namespaces)
    
def psm_from_query(query):
    """Analyze an Element object with a spectrum query and generate a
    dictionary with a peptide spectrum match.    
    """

    # A psm stores the properties of the spectrum query...
    psm = dict(query.attrib)
    # ... and the best hit from peptide database.
    search_hit_elements = query.xpath(
        'd:search_result/d:search_hit[@hit_rank=\'1\']',
        namespaces={'d':xmlns})
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
                                  namespaces={'d':xmlns}):
        psm[subelement.attrib['name']] = float(subelement.attrib['value'])
    for subelement in query.xpath('d:search_result/d:search_hit/'
                                  'd:analysis_result/d:peptideprophet_result',
                                  namespaces={'d':xmlns}):
        psm['peptideprophet'] = float(subelement.attrib['probability'])
    for subelement in query.xpath('d:search_result/d:search_hit/'
                                  'd:altenative_protein',
                                  namespaces={'d':xmlns}):
        proteins.append(subelement.attrib)
    for subelement in query.xpath('d:search_result/d:search_hit/'
                                  'd:modification_info',
                                  namespaces={'d':xmlns}):
        psm['modified_peptide'] = subelement.attrib.get('modified_peptide',
                                                        '')
        if 'mod_nterm_mass' in subelement.attrib:
            modifications.append(
                {'position' : 0,
                 'mass': float(subelement.attrib['mod_nterm_mass'])})
        if 'mod_cterm_mass' in subelement.attrib:
            modifications.append(
                {'position' : peptide_length(psm) + 1,
                 'mass': float(subelement.attrib['mod_cterm_mass'])})
        for mod_element in subelement.xpath('d:mod_aminoacid_mass',
                                            namespaces={'d':xmlns}):
            modification = dict(mod_element.attrib)
            modification['position'] = int(modification['position'])
            modification['mass'] = float(modification['mass'])
            modifications.append(modification)                

    psm["proteins"] = proteins
    psm["modifications"] = modifications

    return psm
    
def iter_psm(source):
    """Parse source and iterate through a list of peptide-spectrum
    matches from the ``source``.

    Arguments:
    source -- any of the following:    
              - a file name/path
              - a file object
              - a file-like object
              - a URL using the HTTP or FTP protocol

    Returns: a list of PSM dicts.
    """

    parser = etree.XMLParser(remove_comments=True, ns_clean=True) 
    tree = etree.parse(source, parser=parser)

    for spectrum_query in tree.xpath(
        '/d:msms_pipeline_analysis/d:msms_run_summary/d:spectrum_query',
        namespaces = {'d': xmlns}):
        
        yield psm_from_query(spectrum_query)

def roc_curve(source):
    """Parse source and return a ROC curve for peptideprophet analysis.

    Arguments:
    source -- any of the following:    
              - a file name/path
              - a file object
              - a file-like object
              - a URL using the HTTP or FTP protocol

    Returns: a list of ROC points, sorted by ascending min prob.
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
    
if __name__ == "__main__":
    pass
