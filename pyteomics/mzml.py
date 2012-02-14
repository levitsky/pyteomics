"""
mzml - reader for mass spectrometry data in mzML format
=======================================================

Summary
-------

mzML is a standard rich XML-format for raw mass spectrometry data storage.
Please refer to http://www.psidev.info/index.php?q=node/257 for the detailed 
specification of the format and the structure of mzML files.

This module provides minimalistic infrastructure for access to data stored in
mzML files. The most important function is :py:func:`iter_spectrum`, which 
reads spectra and related information as saves them into human-readable dicts.
The rest of data can be obtained via a combination of :py:func:`get_node` and
:py:func:`read_params` functions. These functions rely on the terminology of 
the underlying `lxml library <http://lxml.de/>`_. 

Data access
-----------

  :py:func:`iter_spectrum` - iterate through spectra in mzML file. Data from a
  single spectrum are converted to a human-readable dict. Spectra themselves are 
  stored under 'm/z array' and 'intensity array' keys.

  :py:func:`get_node` - get arbitrary nodes of mzML file by their xpath.

  :py:func:`read_params` - read children cvParams and userParams into a dict.

-------------------------------------------------------------------------------

"""
# Licensed under the MIT license:
# http://www.opensource.org/licenses/mit-license.php 

import numpy
import zlib

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

# Default mzML namespace.
xmlns = 'http://psi.hupo.org/ms/mzml'

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
    """Retrieves arbitrary nodes from an mzml file by their xpath.
    Each node in the xpath is assigned to the default namespace 
    'http://psi.hupo.org/ms/mzml' unless specified else.

    Parameters
    ----------
    source : str or file
        A path or an URL to a target mzML file or the file object itself.
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
    >> get_node('/indexedmzML/mzML/run/spectrumList/spectrum')

    """
    parser = etree.XMLParser(remove_comments=True, ns_clean=True)
    tree = etree.parse(source, parser=parser)

    xpath_w_namespace = _insert_default_ns(xpath)

    return tree.xpath(xpath_w_namespace, namespaces=namespaces)

def _decode_base64_data_array(source, dtype, is_compressed):
    """Read a base64-encoded binary array.
    
    Parameters
    ----------
    source : str 
        A binary array encoded with base64.
    dtype : str
        The type of the array in numpy dtype notation.
    is_compressed : bool
        If True then the array will be decompressed with zlib.

    Returns
    -------
    out : numpy.array
    """

    decoded_source = source.decode('base64')
    if is_compressed:
        decoded_source = zlib.decompress(decoded_source)
    output = numpy.frombuffer(decoded_source, dtype=dtype)
    return output

def read_params(element, xpath, namespaces):
    """
    Obtain children nodes of a given node by their xpath and read their
    children cvParam and userParam.

    Parameters
    ----------
    element : lxml.Element
        A parent element.
    xpath : str
        An XPath of children nodes relative to the parent element.
    namespaces : dict
        A dictionary of namespaces. The default namespace key is 'd'.

    Returns
    -------
    out : dict
    """

    output = {}
    for subelement in element.xpath(xpath, namespaces=namespaces):
        for param_element in subelement.xpath('d:cvParam | d:userParam',
                                              namespaces=namespaces):
            output[param_element.attrib['name']] = (
                param_element.attrib['value'])

    return output

def _spectrum_from_element(element, namespaces):
    """Analyze a spectrum Element object and generate a dictionary with its 
    properties.

    Parameters
    ----------
    element : lxml.Element
        A parent element with a spectrum.
    namespaces : dict
        A dictionary of namespaces. The default namespace key is 'd'.

    Returns
    -------
    out : dict
    """

    spectrum = {}
    spectrum.update(element.attrib)    
    spectrum.update(read_params(element, 'self::*', namespaces=namespaces))

    # Read the scan list.
    spectrum.update(
        read_params(element, 'd:scanList', namespaces=namespaces))
    spectrum['scanList'] = []

    for scan_elem in element.xpath('d:scanList/d:scan', namespaces=namespaces):
        scan = read_params(scan_elem, 'self::*', namespaces=namespaces)
        scan['scanWindowList'] = []
        for scan_window_elem in scan_elem.xpath('d:scanWindowList/'
                                                'd:scanWindow',
                                                namespaces=namespaces):
            scan_window = read_params(
                scan_window_elem, 'self::*', namespaces=namespaces) 
            scan['scanWindowList'].append(scan_window)
        spectrum['scanList'].append(scan)                    

    # Read the list of precursors.
    spectrum['precursorList'] = []
    for precursor_elem in element.xpath('d:precursorList/d:precursor',
                                        namespaces=namespaces):
        precursor = {}
        precursor.update(precursor_elem.attrib)
        precursor.update(
            read_params(precursor_elem, 'self::*', namespaces=namespaces))
        precursor.update(
            read_params(precursor_elem, 'd:isolationWindow',
                         namespaces=namespaces))
        precursor.update(
            read_params(precursor_elem, 'd:activation', namespaces=namespaces))

        # Read the list of selected ions for given precursor.
        precursor['selectedIonList'] = []
        for selectedIon_elem in precursor_elem.xpath('d:selectedIonList/'
                                                     'd:selectedIon',
                                                     namespaces=namespaces):
            selectedIon = read_params(
                selectedIon_elem, 'self::*', namespaces=namespaces)
            precursor['selectedIonList'].append(selectedIon)        
        spectrum['precursorList'].append(precursor)

    # Read binary arrays with m/zs and intensities.
    for array_element in element.xpath(
        'd:binaryDataArrayList/d:binaryDataArray', namespaces=namespaces): 
        
        # Define the contents of the array.
        array_type = None
        for cvParam in array_element.xpath(
            'd:cvParam[@name=\'m/z array\'] | '
            'd:cvParam[@name=\'intensity array\']',
            namespaces=namespaces):
            
            array_type = cvParam.attrib['name']

        # Define the type of an array element.
        type_code = None
        for cvParam in array_element.xpath(
            'd:cvParam[@name=\'64-bit float\'] | '
            'd:cvParam[@name=\'32-bit float\']',
            namespaces=namespaces):
            
            type_code = {'32-bit float': 'f',
                         '64-bit float': 'd'
                         }[cvParam.attrib['name']]
            
        # Decode, decompress and read the array.
        is_compressed = bool(
            array_element.xpath('d:cvParam[@name=\'zlib compression\']',
                                namespaces=namespaces))
        
        if array_element.attrib['encodedLength'] != '0':
            spectrum[array_type] = _decode_base64_data_array(
                array_element.xpath('d:binary/text()', namespaces=namespaces)[0],
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
    """Parse ``source`` and iterate through spectra.

    Parameters
    ----------
    source : str or file
        A path or an URL to a target mzML file or the file object itself.

    Returns
    -------
    out : iterator
       An iterator over the dicts with spectra properties.
    """

    parser = etree.XMLParser(remove_comments=True, ns_clean=True)
    tree = etree.parse(source, parser=parser)
    namespaces = {'d': xmlns}

    for spectrum_element in tree.xpath(
        _insert_default_ns('/indexedmzML/mzML/run/spectrumList/spectrum'),
        namespaces=namespaces):
        
        yield _spectrum_from_element(spectrum_element, namespaces=namespaces)
