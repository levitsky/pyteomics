# Licensed under the MIT license:
# http://www.opensource.org/licenses/mit-license.php 

import array
import zlib
import base64

from lxml import etree

# A list of the spectrum attributes which contain float values.
float_keys = [
    'base peak m/z',
    'base peak intensity',
    'total ion current',
    'lowest observed m/z',
    'highest observed m/z',
    'scan start time',
    'ion injection time',
    'scan window lower limit',
    'scan window upper limit',
    'ms level'
    ]

# Default namespace of mzML.
xmlns = 'http://psi.hupo.org/ms/mzml'

def get_node(source, xpath, namespaces={'d':xmlns}):
    """Retrieves arbitrary nodes from an mzml file by their xpath.

    Arguments: 

    source -- any of the following:    
              - a file name/path
              - a file object
              - a file-like object
              - a URL using the HTTP or FTP protocol

    xpath -- a string with XPath to required objects. Usually, mzml
    has a default namespace 'http://psi.hupo.org/ms/mzml' which should
    be prepended to every nodename.

    namespaces -- a dictionary of namespaces. The default pepxml
    namespace has a key 'd'.
    
    Example:
    mzml.get_node('/d:indexedmzML/d:mzML/d:run/d:spectrumList/d:spectrum')

    Returns:
    a list of lxml.Element objects.
    """

    parser = etree.XMLParser(remove_comments=True, ns_clean=True)
    tree = etree.parse(source, parser=parser)
    
    return tree.xpath(xpath, namespaces=namespaces)


def fill_params(output_dict, element, xpath, namespaces={'d':xmlns}):
    """Obtain subelements of the given element with xpath and read their
    children cvParam and userParam to the given dictionary.

    Keyword arguments:
    output_dict -- a dict to fill with params;
    element     -- a parent element;
    xpath       -- an xpath of subelements relative to the parent element;
    namespaces  -- a dict of namespace abbreviations used in the xpath.
    """

    for subelement in element.xpath(xpath, namespaces=namespaces):
        for param_element in subelement.xpath('d:cvParam | d:userParam',
                                              namespaces=namespaces):
            output_dict[param_element.attrib['name']] = (
                param_element.attrib['value'])

def decode_base64_data_array(source, type_code, is_compressed):
    """Read a base64-coded binary array.
    
    Arguments:
    source -- a string with a base64-coded binary array.

    type_code -- a type code. Allowable values:
              'l' - signed long, 4 bytes
              'L' - unsigned long, 4 bytes
              'f' - float, 4 bytes
              'd' - double, 4 bytes

    is_compressed -- if True then the byte array will be decompressed
                     with zlib

    Return: a python array
    """
    if type_code not in ['l', 'L', 'f', 'd']:
        print 'Unexpected type of a variable'
        return []
    
    decoded_source = base64.decodestring(source)
    if is_compressed:
        decoded_source = zlib.decompress(decoded_source)
    output = array.array(type_code)
    output.fromstring(decoded_source)
    return output

def spectrum_from_element(element):
    """Analyze an Element object with a spectrum and generate a
    dictionary with its description.
    """

    spectrum = {}
    spectrum.update(element.attrib)    
    fill_params(spectrum, element, 'self::*')

    # Read the scan list.
    fill_params(spectrum, element, 'd:scanList')
    spectrum['scanList'] = []
    for scan_elem in element.xpath('d:scanList/d:scan',
                              namespaces={'d':xmlns}):
        scan = {}
        fill_params(scan, scan_elem, 'self::*')
        scan['scanWindowList'] = []
        for scan_window_elem in scan_elem.xpath('d:scanWindowList/'
                                                'd:scanWindow',
                                                namespaces={'d':xmlns}):
            scan_window = {}
            fill_params(scan_window, scan_window_elem, 'self::*') 
            scan['scanWindowList'].append(scan_window)
        spectrum['scanList'].append(scan)                    

    # Read the list of precursors.
    spectrum['precursorList'] = []
    for precursor_elem in element.xpath('d:precursorList/d:precursor',
                                        namespaces={'d':xmlns}):
        precursor = {}
        precursor.update(precursor_elem.attrib)
        fill_params(precursor, precursor_elem, 'self::*')
        fill_params(precursor, precursor_elem, 'd:isolationWindow')
        fill_params(precursor, precursor_elem, 'd:activation')

        # Read the list of selected ions for given precursor.
        precursor['selectedIonList'] = []
        for selectedIon_elem in precursor_elem.xpath('d:selectedIonList/'
                                                     'd:selectedIon',
                                                     namespaces={'d':xmlns}):
            selectedIon = {}
            fill_params(selectedIon, selectedIon_elem, 'self::*')
            precursor['selectedIonList'].append(selectedIon)        
        spectrum['precursorList'].append(precursor)

    # Read binary arrays with m/zs and intensities.
    for array_element in element.xpath(
        'd:binaryDataArrayList/d:binaryDataArray',
        namespaces={'d':xmlns}): 
        
        # Define the contents of the array.
        array_type = None
        for cvParam in array_element.xpath(
            'd:cvParam[@name=\'m/z array\'] | '
            'd:cvParam[@name=\'intensity array\']',
            namespaces={'d':xmlns}):
            
            array_type = cvParam.attrib['name']

        # Define the type of an array element.
        type_code = None
        for cvParam in array_element.xpath(
            'd:cvParam[@name=\'64-bit float\'] | '
            'd:cvParam[@name=\'32-bit float\']',
            namespaces={'d':xmlns}):
            
            type_code = {'32-bit float': 'f',
                         '64-bit float': 'd'
                         }[cvParam.attrib['name']]
            
        # Decode, decompress and read the array.
        is_compressed = bool(
            array_element.xpath('d:cvParam[@name=\'zlib compression\']',
                                namespaces={'d':xmlns}))
        
        if array_element.attrib['encodedLength'] != '0':
            spectrum[array_type] = decode_base64_data_array(
                array_element.xpath('d:binary/text()',
                                    namespaces={'d':xmlns})[0],
                type_code,
                is_compressed)
        else:
            spectrum[array_type] = []

    # Convert selected values to float.
    for key in spectrum:
        if key in float_keys:
            spectrum[key] = float(spectrum[key])

    return spectrum

def iter_spectrum(source):
    """Parse source and iterate through a list of peptide-spectrum
    matches from the ``source``.

    Arguments:
    source -- any of the following:    
              - a file name/path
              - a file object
              - a file-like object
              - a URL using the HTTP or FTP protocol

    Returns: a generator which yields PSM dicts one by one.
    """

    parser = etree.XMLParser(remove_comments=True, ns_clean=True)
    tree = etree.parse(source, parser=parser)

    for spectrum_element in tree.xpath(
        '/d:indexedmzML/d:mzML/d:run/d:spectrumList/d:spectrum',
        namespaces = {'d': xmlns}):
        
        yield spectrum_from_element(spectrum_element)
