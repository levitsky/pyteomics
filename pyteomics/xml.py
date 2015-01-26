import re
import warnings
warnings.formatwarning = lambda msg, *args: str(msg) + '\n'
import socket
from functools import wraps
from traceback import format_exc
import operator as op
import ast
import numpy as np
from lxml import etree
from .auxiliary import _file_obj, PyteomicsError
try: # Python 2.7
    from urllib2 import urlopen, URLError
except ImportError: # Python 3.x
    from urllib.request import urlopen, URLError


def _local_name(element):
    """Strip namespace from the XML element's name"""
    if element.tag and element.tag[0] == '{':
        return element.tag.rsplit('}', 1)[1]
    return element.tag


def _keepstate(func):
    """Decorator to help keep the position in open files passed as
    positional arguments to functions"""
    @wraps(func)
    def wrapped(self, *args, **kwargs):
        position = self.tell()
        self.seek(0)
        res = func(self, *args, **kwargs)
        self.seek(position)
        return res
    return wrapped


class XMLValueConverter(object):
    @classmethod
    def str_to_bool(cls, s):
        if s.lower() in {'true', '1'}:
            return True
        if s.lower() in {'false', '0'}:
            return False
        raise PyteomicsError('Cannot convert string to bool: ' + s)

    @classmethod
    def str_to_num(cls, s, numtype):
        return numtype(s) if s else None

    @classmethod
    def to(cls, t):
        def convert_from(s):
            return cls.str_to_num(s, t)
        return convert_from

    @classmethod
    def converters(cls):
        return {
            'ints': cls.to(int), 'floats': cls.to(float), 'bools': cls.str_to_bool,
            'intlists': lambda x: np.fromstring(x.replace('\n', ' '), dtype=int, sep=' '),
            'floatlists': lambda x: np.fromstring(x.replace('\n', ' '), sep=' '),
            'charlists': list
        }


class XMLParserBase(object):
    # Configurable data
    file_format = "XML"
    root_element = None
    version_info_element = None
    default_schema = {}
    default_version = 0
    default_iter_tag = None
    structures_to_flatten = []

    # Configurable plugin logic
    converters = XMLValueConverter.converters()

    # Must be implemented by subclasses
    def get_info_smart(self, element, **kwargs):
        raise NotImplementedError

    def __init__(self, source, read_schema=True, **kwargs):
        self.source = _file_obj(source, 'rb')

        if kwargs.pop('build_tree', False):
            self.build_tree()
        else:
            self.tree = None

        # For handling
        self.version_info = self.get_version_info(self.version_info_element)
        self.schema_info = self.get_schema_info(read_schema)
        self.source.file.seek(0)

        self._iter = self.iterfind(self.default_iter_tag, **kwargs)

    def __iter__(self):
        return self._iter

    def __next__(self):
        try:
            return next(self._iter)
        except StopIteration:
            self.__exit__(None, None, None)
            raise

    next = __next__  # Python 2 support

    # context manager support
    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        self.source.__exit__(*args, **kwargs)

    def __getattr__(self, attr):
        return getattr(self.source, attr)

    @_keepstate
    def get_version_info(self, element):
        """
        Provide version information about the {0} file.

        Parameters:
        -----------
        source : str or file
            file object or path to file

        Returns:
        --------
        out : tuple
            A (version, schema URL) tuple, both elements are strings or None.
        """
        vinfo = None
        for _, elem in etree.iterparse(
                self.source, events=('start',), remove_comments=True):
            if _local_name(elem) == element:
                return (elem.attrib.get('version'),
                        elem.attrib.get(('{{{}}}'.format(elem.nsmap['xsi'])
                            if 'xsi' in elem.nsmap else '') + 'schemaLocation'))

    @_keepstate
    def get_schema_info(self, read_schema=True):
        """Stores defaults for the schema, tries to retrieve the schema for
        other versions. Keys are: 'floats', 'ints', 'bools', 'lists',
        'intlists', 'floatlists', 'charlists'."""
        if not read_schema:
            return self.default_schema

        version, schema = self.version_info
        if version == self.default_version:
            return self.default_schema

        ret = {}
        try:
            if not schema:
                schema_url = ''
                raise PyteomicsError(
                        'Schema information not found in {}.'.format(
                            self.source.name))
            schema_url = schema.split()[-1]
            if not (schema_url.startswith('http://') or
                    schema_url.startswith('file://')):
                schema_url = 'file://' + schema_url
            schema_file = urlopen(schema_url)
            p = etree.XMLParser(remove_comments=True)
            schema_tree = etree.parse(schema_file, parser=p)
            types = {'ints': {'int', 'long', 'nonNegativeInteger', 'positiveInt',
                              'integer', 'unsignedInt'},
                     'floats': {'float', 'double'},
                     'bools': {'boolean'},
                     'intlists': {'listOfIntegers'},
                     'floatlists': {'listOfFloats'},
                     'charlists': {'listOfChars', 'listOfCharsOrAny'}}
            for k, val in types.items():
                tuples = set()
                for elem in schema_tree.iter():
                    if _local_name(elem) == 'attribute' and elem.attrib.get(
                            'type', '').split(':')[-1] in val:
                        anc = elem.getparent()
                        while not (
                                (_local_name(anc) == 'complexType'
                                    and 'name' in anc.attrib)
                                or _local_name(anc) == 'element'):
                            anc = anc.getparent()
                            if anc is None:
                                break
                        else:
                            if _local_name(anc) == 'complexType':
                                elnames = [x.attrib['name'] for x in
                                           schema_tree.iter()
                                           if x.attrib.get('type', ''
                                               ).split(':')[-1] == anc.attrib['name']]
                            else:
                                elnames = (anc.attrib['name'],)
                            for elname in elnames:
                                tuples.add(
                                    (elname, elem.attrib['name']))
                ret[k] = tuples
            ret['lists'] = set(elem.attrib['name'] for elem in schema_tree.xpath(
                '//*[local-name()="element"]') if 'name' in elem.attrib and
                elem.attrib.get('maxOccurs', '1') != '1')
        except Exception as e:
            if isinstance(e, (URLError, socket.error, socket.timeout)):
                warnings.warn("Can't get the {0} schema for version "
                "`{1}` from <{2}> at the moment.\n"
                "Using defaults for {3}.\n"
                "You can disable reading the schema by specifying "
                "`read_schema=False`.".format(self.file_format, version,
                    schema_url, self.default_version))
            else:
                warnings.warn("Unknown {0} version `{1}`. "
                    "Attempt to use schema\n"
                    "information from <{2}> failed.\n"
                    "Exception information:\n{3}\n"
                    "Falling back to defaults for {0[default_version]}\n"
                    "NOTE: This is just a warning, probably from a badly-"
                    "generated XML file.\nYou will still most probably get "
                    "decent results.\nLook here for suppressing warnings:\n"
                    "http://docs.python.org/library/warnings.html#"
                    "temporarily-suppressing-warnings\n"
                    "You can also disable reading the schema by specifying "
                    "`read_schema=False`.\n"
                    "If you think this shouldn't have happened, please "
                    "report this to\n"
                    "http://hg.theorchromo.ru/pyteomics/issues\n"
                    "".format(self.file_format, version, schema_url,
                        format_exc()))
            ret = self.default_schemau
        return ret

    def handle_param(self, element, **kwargs):
        '''Unpacks cvParam and userParam tags into key-value pairs'''
        if 'value' in element.attrib:
            try:
                value = float(element.attrib['value'])
            except ValueError:
                value = element.attrib['value']
            return {element.attrib['name']: value}
        else:
            return {'name': element.attrib['name']}

    def get_info(self, element, **kwargs):
        """Extract info from element's attributes, possibly recursive.
        <cvParam> and <userParam> elements are treated in a special way."""
        name = _local_name(element)
        schema_info = self.schema_info
        if name in {'cvParam', 'userParam'}:
            return self.handle_param(element)

        info = dict(element.attrib)
        # process subelements
        if kwargs.get('recursive'):
            for child in element.iterchildren():
                cname = _local_name(child)
                if cname in {'cvParam', 'userParam'}:
                    newinfo = self.handle_param(child, **kwargs)
                    if not ('name' in info and 'name' in newinfo):
                        info.update(newinfo)
                    else:
                        if not isinstance(info['name'], list):
                            info['name'] = [info['name']]
                        info['name'].append(newinfo.pop('name'))
                else:
                    if cname not in schema_info['lists']:
                        info[cname] = self.get_info_smart(child, **kwargs)
                    else:
                        info.setdefault(cname, []).append(
                                self.get_info_smart(child, **kwargs))

        # process element text
        if element.text and element.text.strip():
            stext = element.text.strip()
            if stext:
                if info:
                    info[name] = stext
                else:
                    return stext

        # convert types
        converters = self.converters
        for k, v in info.items():
            for t, a in converters.items():
                if (_local_name(element), k) in schema_info[t]:
                    info[k] = a(v)

        # resolve refs
        if kwargs.get('retrieve_refs'):
            self._retrieve_refs(info, **kwargs)

        # flatten the excessive nesting
        for k, v in dict(info).items():
            if k in self.structures_to_flatten:
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

    @_keepstate
    def build_tree(self):
        '''Build and store the ElementTree instance for `self.source`'''
        p = etree.XMLParser(remove_comments=True)
        self.tree = etree.parse(self.source, parser=p)

    def clear_tree(self):
        self.tree = None

    @_keepstate
    def iterfind(self, path, **kwargs):
        """Parse `source` and yield info on elements with specified local name
        or by specified "XPath". Only local names separated with slashes are
        accepted. An asterisk (`*`) means any element.
        You can specify a single condition in the end, such as:
        "/path/to/element[some_value>1.5]"
        Note: you can do much more powerful filtering using plain Python.
        The path can be absolute or "free". Please don't specify
        namespaces."""
        try:
            path, _, cond = re.match(pattern_path, path).groups()
        except AttributeError:
            raise PyteomicsError('Invalid path: ' + path)
        if path.startswith('//') or not path.startswith('/'):
            absolute = False
            if path.startswith('//'):
                path = path[2:]
                if path.startswith('/') or '//' in path:
                    raise PyteomicsError("Too many /'s in a row.")
        else:
            absolute = True
            path = path[1:]
        nodes = path.rstrip('/').split('/')
        if not self.tree:
            localname = nodes[0].lower()
            found = False
            for ev, elem in etree.iterparse(self, events=('start', 'end'),
                    remove_comments=True):
                name_lc = _local_name(elem).lower()
                if ev == 'start':
                    if name_lc == localname or localname == '*':
                        found += True
                else:
                    if name_lc == localname or localname == '*':
                        if (absolute and elem.getparent() is None) or not absolute:
                            for child in get_rel_path(elem, nodes[1:]):
                                info = self.get_info_smart(child, **kwargs)
                                if cond is None or satisfied(info, cond):
                                    yield info
                        if not localname == '*':
                            found -= 1
                    if not found:
                        elem.clear()
        else:
            xpath = ('/' if absolute else '//') + '/'.join(
                    '*[local-name()="{}"]'.format(node) for node in nodes)
            for elem in self.tree.xpath(xpath):
                info = self.get_info_smart(elem, tree=self.tree, **kwargs)
                if cond is None or satisfied(info, cond):
                    yield info

# XPath emulator tools
pattern_path = re.compile('([\w/*]*)(\[(\w+[<>=]{1,2}[^\]]+)\])?')
pattern_cond = re.compile('^\s*(\w+)\s*([<>=]{,2})\s*([^\]]+)$')

def get_rel_path(element, names):
    if not names:
        yield element
    else:
        for child in element.iterchildren():
            if _local_name(child).lower() == names[0].lower(
                    ) or names[0] == '*':
                if len(names) == 1:
                    yield child
                else:
                    for gchild in get_rel_path(child, names[1:]):
                        yield gchild

def satisfied(d, cond):
    func = {'<': 'lt', '<=': 'le', '=': 'eq', '==': 'eq',
            '!=': 'ne', '<>': 'ne', '>': 'gt', '>=': 'ge'}
    try:
        lhs, sign, rhs = re.match(pattern_cond, cond).groups()
        if lhs in d:
            return getattr(op, func[sign])(
                    d[lhs], ast.literal_eval(rhs))
        return False
    except (AttributeError, KeyError, ValueError):
        raise PyteomicsError('Invalid condition: ' + cond)

_mzid_schema_defaults = {'ints': {('DBSequence', 'length'),
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


