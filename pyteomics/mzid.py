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
import numpy

from . import auxiliary as aux

# 'keys' should contain keys whose value is a dict
def _make_get_info(env):
    def _get_info(source, element, recursive=False, retrieve_refs=False):
        """Extract info from element's attributes, possibly recursive.
        <cvParam> and <userParam> elements are treated in a special way."""
        name = aux._local_name(element)
        kwargs = dict(recursive=recursive, retrieve_refs=retrieve_refs)
        if name in {'cvParam', 'userParam'}:
            if 'value' in element.attrib:
                try:
                    value = float(element.attrib['value'])
                except ValueError:
                    value = element.attrib['value']
                return {element.attrib['name']: value}
            else:
                return {'name': element.attrib['name']}

        info = dict(element.attrib)
        # process subelements
        if recursive:
            for child in element.iterchildren():
                cname = aux._local_name(child)
                if cname in {'cvParam', 'userParam'}:
                    info.update(_get_info(source, child))
                else:
                    if cname not in env['schema_info'](source, 'lists'):
                        info[cname] = _get_info_smart(source, child, **kwargs)
                    else:
                        if cname not in info:
                            info[cname] = []
                        info[cname].append(_get_info_smart(source, child, **kwargs))
        # process element text
        if element.text and element.text.strip():
            stext = element.text.strip()
            if stext:
                if info:
                    info[name] = stext
                else:
                    return stext
        # convert types
        def str_to_bool(s):
            if s.lower() in {'true', '1'}: return True
            if s.lower() in {'false', '0'}: return False
            raise aux.PyteomicsError('Cannot convert string to bool: ' + s)

        converters = {'ints': int, 'floats': float, 'bools': str_to_bool,
                'intlists': lambda x: numpy.fromstring(x, dtype=int, sep=' '),
                'floatlists': lambda x: numpy.fromstring(x, sep=' '),
                'charlists': list}
        for k, v in info.items():
            for t, a in converters.items():
                if (aux._local_name(element), k) in env['schema_info'](source, t):
                    info[k] = a(v)
        # resolve refs
        # loop is needed to resolve refs pulled from other refs
        if retrieve_refs:
            while True:
                refs = False
                for k, v in dict(info).items():
                    if k.endswith('_ref'):
                        refs = True
                        info.update(get_by_id(source, v))
                        del info[k]
                        del info['id']
                if not refs:
                    break
        # flatten the excessive nesting
        for k, v in dict(info).items():
            if k in env['keys']:
                info.update(v)
                del info[k]
        # another simplification
        for k, v in dict(info).items():
            if isinstance(v, dict) and 'name' in v and len(v) == 1:
                info[k] = v['name']
        if len(info) == 2 and 'name' in info and (
                'value' in info or 'values' in info):
            name = info.pop('name')
            info = {name: info.popitem()[1]}
        return info
    return _get_info

def _get_info_smart(source, element, **kw):
    """Extract the info in a smart way depending on the element type"""
    name = aux._local_name(element)
    kwargs = dict(kw)
    rec = kwargs.pop('recursive', None)
    if name == 'MzIdentML':
        return _get_info(source, element, rec if rec is not None else False,
                **kwargs)
    else:
        return _get_info(source, element, rec if rec is not None else True,
                **kwargs)

@aux._keepstate
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
    for event, elem in etree.iterparse(source, events=('start', 'end'),
            remove_comments=True):
        if event == 'start':
            if elem.attrib.get('id') == elem_id:
                found = True
        else:
            if elem.attrib.get('id') == elem_id:
                return _get_info_smart(source, elem)
            if not found:
                elem.clear()
    return None

def _make_itertag(env):
    @aux._keepstate
    def _itertag(source, localname, **kwargs):
        """Parse ``source`` and yield info on elements with specified local name.
        Case-insensitive. Namespace-aware."""
        found = False
        for ev, elem in etree.iterparse(source, events=('start', 'end'),
                remove_comments=True):
            if ev == 'start':
                if aux._local_name(elem).lower() == localname.lower():
                    found = True
            else:
                if aux._local_name(elem).lower() == localname.lower():
                    yield env['get_info_smart'](source, elem, **kwargs)
                    found = False
                    if not found:
                        elem.clear()
    return _itertag

_schema_defaults = {'ints': {('DBSequence', 'length'),
                     ('IonType', 'charge'),
                     ('BibliographicReference', 'year'),
                     ('SubstitutionModification', 'location'),
                     ('PeptideEvidence', 'end'),
                     ('Enzyme', 'missedCleavages'),
                     ('PeptideEvidence', 'start'),
                     ('Modification', 'location'),
                     ('SpectrumIdentificationItem', 'rank'),
                     ('SpectrumIdentificationItem', 'chargeState')},
            'floats': {('SubstitutionModification', 'monoisotopicMassDelta'),
                     ('SpectrumIdentificationItem', 'experimentalMassToCharge'),
                     ('Residue', 'mass'),
                     ('SpectrumIdentificationItem', 'calculatedPI'),
                     ('Modification', 'avgMassDelta'),
                     ('SearchModification', 'massDelta'),
                     ('Modification', 'monoisotopicMassDelta'),
                     ('SubstitutionModification', 'avgMassDelta'),
                     ('SpectrumIdentificationItem', 'calculatedMassToCharge')},
            'bools': {('PeptideEvidence', 'isDecoy'),
                     ('SearchModification', 'fixedMod'),
                     ('Enzymes', 'independent'),
                     ('Enzyme', 'semiSpecific'),
                     ('SpectrumIdentificationItem', 'passThreshold'),
                     ('ProteinDetectionHypothesis', 'passThreshold')},
            'lists': {'SourceFile', 'SpectrumIdentificationProtocol',
                    'ProteinDetectionHypothesis', 'SpectraData', 'Enzyme',
                    'Modification', 'MassTable', 'DBSequence',
                    'InputSpectra', 'cv', 'IonType', 'SearchDatabaseRef',
                    'Peptide', 'SearchDatabase', 'ContactRole', 'cvParam',
                    'ProteinAmbiguityGroup', 'SubSample',
                    'SpectrumIdentificationItem', 'TranslationTable',
                    'AmbiguousResidue', 'SearchModification',
                    'SubstitutionModification', 'PeptideEvidenceRef',
                    'PeptideEvidence', 'SpecificityRules',
                    'SpectrumIdentificationResult', 'Filter', 'FragmentArray',
                    'InputSpectrumIdentifications', 'BibliographicReference',
                    'SpectrumIdentification', 'Sample', 'Affiliation',
                    'PeptideHypothesis',
                    'Measure', 'SpectrumIdentificationItemRef'},
            'intlists': {('IonType', 'index'), ('MassTable', 'msLevel')},
            'floatlists': {('FragmentArray', 'values')},
            'charlists': {('Modification', 'residues'),
                    ('SearchModification', 'residues')}}

_version_info_env = {'format': 'mzIdentML', 'element': 'MzIdentML'}
version_info = aux._make_version_info(_version_info_env)

_schema_env = {'format': 'MzIdentML', 'version_info': version_info,
        'default_version': '1.1.0', 'defaults': _schema_defaults}
_schema_info = aux._make_schema_info(_schema_env)

_get_info_env = {'keys':  {'Fragmentation',},
        'schema_info': _schema_info}
_get_info = _make_get_info(_get_info_env)

_itertag_env = {'get_info_smart': _get_info_smart}
_itertag = _make_itertag(_itertag_env)

def read(source, **kwargs):
    """Parse ``source`` and iterate through peptide-spectrum matches.

    Parameters
    ----------
    source : str or file
        A path to a target mzIdentML file or the file object itself.
    recursive : bool, optional
        If :py:const:`False`, subelements will not be processed when
        extracting info from elements. Default is :py:const:`True`.
    retrieve_refs : bool, optional
        If :py:const:`True`, additional information from references will be
        automatically added to the results. The file processing time will
        increase. Default is :py:const:`False`.

    Returns
    -------
    out : iterator
       An iterator over the dicts with PSM properties.
    """
    
    return _itertag(source, 'SpectrumIdentificationResult', **kwargs)
