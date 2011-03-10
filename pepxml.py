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

def name(element):
    return str(element.tag).split('}')[1]

def psm_from_query(query):
    """Analyze an Element object with spectrum query and generate a
    dictionary with a peptide spectrum match.    
    """

    psm = dict(query.attrib)         #attributes of <spectrum_query>
    psm.update(query[0][0].attrib)   #attributes of <search_hit>

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
    for child in query[0][0]:
        if name(child) == 'search_score':
            psm[child.attrib['name']] = float(child.attrib['value'])
        elif name(child) == 'alternative_protein':
            proteins.append(child.attrib)
        elif name(child) == 'modification_info':
            if len(child):
                modifications.append(child[0].attrib)
        elif name(child) == 'analysis_result':
            if child.attrib.get('analysis', '') == 'peptideprophet':
                psm['peptideprophet'] = float(child[0].attrib['probability'])
        else:
            print "Unexpected child named %s" % name(child)

    psm["proteins"] = proteins
    psm["modifications"] = modifications

    # Generate the modified sequence of a peptide.
    shift = 0
    psm["non_mod_seq"] = psm["peptide"]
    if len(psm["modifications"]):
        for mod in psm["modifications"]:
            psm["peptide"] = (
                psm["peptide"][:int(mod["position"]) + shift] + '[' 
                + mod["mass"].split('.')[0] + ']' 
                + psm["peptide"][int(mod["position"]) + shift:])
            shift += len(mod["mass"]) + 2
    return psm
    
def psm_list(source):
    """Parse source and return a list of peptide-spectrum matches from
    the ``source``.

    The ``source`` can be any of the following:    
    - a file name/path
    - a file object
    - a file-like object
    - a URL using the HTTP or FTP protocol
    """

    parser = etree.XMLParser(remove_comments=True, ns_clean=True) 
    tree = etree.parse(source, parser=parser)

    PSMs = []
    for elem in tree.iter():
        if name(elem) == "spectrum_query":
            PSMs.append(psm_from_query(elem))
    return PSMs
    
if __name__ == "__main__":
    pass
