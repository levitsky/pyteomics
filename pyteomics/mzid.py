"""
mzid - reader for peptide identification data in mzIdentML format
======================================================================

Summary
-------

`mzIdentML <http://www.psidev.info/mzidentml>`_  is one of the standards
developed by the Proteomics Informatics working group of the HUPO Proteomics
Standard Initiative.

This module provides a minimalistic way to extract information from mzIdentML
files. The main idea is the same as in :py:mod:`pyteomics.pepxml`: the top-level
function :py:func:`read` allows iterating over entries in
`<SpectrumIdentificationResult>` tags, i.e. groups of identifications
for a certain spectrum. Note that each entry can contain more than one PSM
(peptide-spectrum match). They are accessible with `"search_hits"` key, just
like in :py:mod:`pyteomics.pepxml`.

Data access
-----------

  :py:func:`read` - iterate through peptide-spectrum matches in a pep.XML 
  file. Data from a single PSM group are converted to a human-readable dict. 

  :py:func:`get_by_id` - get an element by its ID and extract the data from it.

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
from functools import wraps
from warnings import warn
from collections import defaultdict
try: # Python 2.7
    from urllib import urlopen
except ImportError: # Python 3.x
    from urllib.request import urlopen
from .auxiliary import PyteomicsError

def _keepstate(func):
    """Decorator to help keep the position in open files passed as arguments"""
    @wraps(func)
    def wrapped(source, *args, **kwargs):
        if hasattr(source, 'seek') and hasattr(source, 'tell'):
            pos = source.tell()
            source.seek(0)
            res = func(source, *args, **kwargs)
            source.seek(pos)
            return res
        else:
            return func(source, *args, **kwargs)
    return wrapped

def _local_name(element):
    if element.tag.startswith('{'):
        return element.tag.rsplit('}', 1)[1]
    else:
        return element.tag

def _get_info(source, element, recursive=False):
    """Extract info from element's attributes, possibly recursive.
    <cvParam> and <userParam> elements are treated in a special way."""
    name = _local_name(element)
    if name in ('cvParam', 'userParam'):
        if 'value' in element.attrib:
            return {element.attrib['name']: element.attrib['value']}
        else:
            return {'name': element.attrib['name']}
    info = dict(element.attrib)
    # process subelements
    if recursive:
        for child in element.iterchildren():
            cname = _local_name(child)
            if cname in ('cvParam', 'userParam'):
                info.update(_get_info(source, child))
            else:
                if cname not in _schema_info(source, 'lists'):
                    info[cname] = _get_info(source, child, True)
                else:
                    if cname not in info:
                        info[cname] = []
                    info[cname].append(_get_info(source, child, True))
    # process element text
    if element.text and element.text.strip():
        stext = element.text.strip()
        if stext:
            if info:
                info[name] = stext
            else:
                return stext
    # convert types
    for k, v in info.items():
        if k in _schema_info(source, 'ints'):
            info[k] = int(v)
        elif k in _schema_info(source, 'floats'):
            info[k] = float(v)
        
    return info

def _get_info_smart(source, element):
    """Extract the info in a smart way depending on the element type"""
    name = _local_name(element)
    if name == 'MzIdentML':
        return _get_info(source, element, False)
    else:
        return _get_info(source, element, True)

@_keepstate
def get_by_id(source, elem_id):
    """Parse ``source`` and return the element with `id` attribute equal to
    ``elem_id``. Returns :py:const:`None` if no such element is found.

    Parameters
    ----------
    source : str or file
        A path to a target mzIdentML file of the file object itself.
    
    elem_id : str
        The value of the `id` attribute to match.

    Returns
    -------
    out : :py:class:`lxml.etree.Element` or :py:const:`None`
    """

    found = False
    for event, elem in etree.iterparse(source, events=('start', 'end')):
        if event == 'start':
            if elem.attrib.get('id') == elem_id:
                found = True
        else:
            if elem.attrib.get('id') == elem_id:
                return _get_info_smart(source, elem)
            if not found:
                elem.clear()
    return None

def read(source):
    """Parse ``source`` and iterate through peptide-spectrum matches.

    Parameters
    ----------
    source : str or file
        A path to a target mzIdentML file or the file object itself.

    Returns
    -------
    out : iterator
       An iterator over the dicts with PSM properties.
    """

    return _itertag(source, 'SpectrumIdentificationResult')

@_keepstate
def _itertag(source, localname):
    """Parse ``source`` and yield info on elements with specified local name.
    Case-insensitive. Namespace-aware."""
    found = False
    for ev, elem in etree.iterparse(source, events=('start', 'end')):
        if ev == 'start':
            if _local_name(elem).lower() == localname.lower():
                found = True
        else:
            if _local_name(elem).lower() == localname.lower():
                yield _get_info_smart(source, elem)
                found = False
                if not found:
                    elem.clear()

@_keepstate
def mzid_version_info(source):
    for _, elem in etree.iterparse(source, events=('start',)):
        if _local_name(elem) == 'MzIdentML':
            return elem.attrib.get('version'), elem.attrib.get((
                '{{{}}}'.format(elem.nsmap['xsi'])
                if 'xsi' in elem.nsmap else '') + 'schemaLocation')

_schema_info_cache = defaultdict(dict)
def _schema_info(source, key):
    '''Stores defaults for version 1.1.0, tries to retrieve the schema for
    other versions. Keys are: 'floats', 'ints', 'lists'.'''

    if key in _schema_info_cache[source]:
        return _schema_info_cache[source][key]
    
    version, schema = mzid_version_info(source)
    defaults = {'ints': ('charge', 'chargeState', 'rank', 'location',
                'missedCleavages', 'start', 'end', 'length', 'year'),
            'floats': ('massDelta', 'experimentalMassToCharge',
                  'calculatedMassToCharge', 'calculatedPI', 'avgMassDelta',
                  'monoisotopicMassDelta', 'mass'),
            'lists': ('Residue', 'AnalysisSoftware', 'SpectrumIdentificationList',
                'SourceFile', 'SpectrumIdentificationProtocol',
                'ProteinDetectionHypothesis', 'SpectraData', 'Enzyme',
                'Modification', 'MassTable', 'DBSequence',
                'InputSpectra', 'cv', 'IonType', 'SearchDatabaseRef',
                'Peptide', 'SearchDatabase', 'ContactRole', 'cvParam',
                'ProteinAmbiguityGroup', 'SubSample', 'SpectrumIdentificationItem',
                'TranslationTable', 'AmbiguousResidue', 'SearchModification',
                'SubstitutionModification', 'PeptideEvidenceRef',
                'PeptideEvidence', 'SpecificityRules',
                'SpectrumIdentificationResult', 'Filter', 'FragmentArray',
                'InputSpectrumIdentifications', 'BibliographicReference',
                'SpectrumIdentification', 'Sample', 'Affiliation',
                'PeptideHypothesis',
                'Measure', 'SpectrumIdentificationItemRef')}
    if version == '1.1.0':
        ret = defaults[key]
    else:
        try:
            if not schema:
                schema_url = ''
                raise PyteomicsError(
                        'Schema information not found in {}.'.format(source))
            schema_url = schema.split()[-1]
            schema_file = urlopen(schema_url)
            if key == 'ints':
                ret = set(elem['name'] for elem in _itertag(
                    schema_file, 'attribute') if elem.get('type') == 'xsd:int')
            elif key == 'floats':
                ret = set(elem['name'] for elem in _itertag(
                    schema_file, 'attribute') if elem.get('type') in (
                        'xsd:float', 'xsd:double'))
            elif key == 'lists':
                ret = set(elem['name'] for elem in _itertag(
                    schema_file, 'element') if elem.get('maxOccurs', '1') != '1')
            else:
                raise PyteomicsError('Unknown key ' + key)

        except Exception, e:
            warn("Unknown MzIdentML version `{}`. Attempt to use schema "
                    "information from <{}> failed. Reason:\n{!r}: {}\n"
                    "Falling back to defaults for 1.1.0".format(
                        version, schema_url, type(e), e.message))
            ret = defaults[key]
    _schema_info_cache[source][key] = ret
    return ret
