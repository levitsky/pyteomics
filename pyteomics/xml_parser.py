import re
import warnings

from functools import wraps
from traceback import format_exc

import numpy as np
from lxml import etree

from pyteomics.auxiliary import _file_obj, PyteomicsError, urlopen, URLError, _mzid_schema_defaults


def _local_name(element):
    """Strip namespace from the XML element's name"""
    if element.tag and element.tag[0] == '{':
        return element.tag.rsplit('}', 1)[1]
    return element.tag


def _oo_keepstate(func):
    """Decorator to help keep the position in open files passed as
    positional arguments to functions"""
    @wraps(func)
    def wrapped(self, *args, **kwargs):
        self._position = self.tell()
        self.seek(0)

        try:
            res = func(self, *args, **kwargs)
        except Exception as e:
            self.source.file.seek(self._position)
            raise e

        self.seek(self._position)
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

    def __init__(self, source, read_schema=True, **kwargs):
        self.source = _file_obj(source, 'rb')

        # Set by @_oo_keepstate
        self._position = None

        self.id_dict = {}
        self.tree = None

        # For handling
        self.version_info = self.get_version_info(self.version_info_element)
        self.schema_info = self.get_schema_info(read_schema)
        self.source.file.seek(0)

        self.retrieve_refs = kwargs.pop("retrieve_refs", False)
        self.iterative = kwargs.pop("iterative", not self.retrieve_refs)

    def seek(self, i):
        self.source.file.seek(i)

    def tell(self):
        return self.source.file.tell()

    @_oo_keepstate
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
        for _, elem in etree.iterparse(self.source, events=('start',), remove_comments=True):
            if _local_name(elem) == element:
                vinfo = (elem.attrib.get('version'), elem.attrib.get(('{{{}}}'.format(
                    elem.nsmap['xsi']) if 'xsi' in elem.nsmap else '') + 'schemaLocation'))
                break
        return vinfo

    @_oo_keepstate
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

    def _retrieve_refs(self, info, tree=None):
        '''Retrieves and embeds the data for each attribute in `info` that
        ends in _ref. Removes the id attribute from `info`'''
        for k, v in dict(info).items():
            if k.endswith('_ref'):
                info.update(self.get_by_id(v, retrieve_refs=True, tree=tree))
                del info[k]
                info.pop('id', None)

    @_oo_keepstate
    def get_info(self, element, recursive=False, retrieve_refs=False, **kw):
        """Extract info from element's attributes, possibly recursive.
        <cvParam> and <userParam> elements are treated in a special way."""
        name = _local_name(element)
        kwargs = dict(recursive=recursive, retrieve_refs=retrieve_refs)
        kwargs.update(kw)
        schema_info = self.schema_info
        if name in {'cvParam', 'userParam'}:
            return self.handle_param(element)

        info = dict(element.attrib)
        # process subelements
        if recursive:
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
        if retrieve_refs:
            self._retrieve_refs(info, kw.get("tree"))

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

    def get_info_smart(self, element, **kwargs):
        """Extract the info in a smart way depending on the element type"""
        name = _local_name(element)
        kwargs = dict(kwargs)
        rec = kwargs.pop("recursive", None)

        # Try not to recursively unpack the root element
        # unless the user really wants to.
        if name == self.root_element:
            return self.get_info(element, rec if rec is not None else False,
                                 **kwargs)
        else:
            return self.get_info(element, rec if rec is not None else True,
                                 **kwargs)

    @_oo_keepstate
    def _build_id_cache(self):
        '''Constructs a cache for each element in the document, indexed by id
        attribute'''
        self.id_dict = None
        stack = 0
        id_dict = {}
        for event, elem in etree.iterparse(self.source, events=('start', 'end'),
                remove_comments=True):
            if event == 'start':
                if 'id' in elem.attrib:
                    stack += 1
            else:
                if 'id' in elem.attrib:
                    stack -= 1
                    id_dict[elem.attrib['id']] = elem
                elif stack == 0:
                    elem.clear()
        self.id_dict = id_dict

    @_oo_keepstate
    def get_by_id(self, elem_id, tree=None, **kwargs):
        """Parse ``source`` and return the element with `id` attribute equal to
        `elem_id`. Returns :py:const:`None` if no such element is found.

        Parameters
        ----------
        source : str or file
            A path to a target mzIdentML file of the file object itself.

        elem_id : str
            The value of the `id` attribute to match.

        Returns
        -------
        out : :py:class:`dict` or :py:const:`None`
        """
        # Carried tree over from functional implementation. A better name to describe what it
        # means? If None, don't use id_dict
        if self.tree is None and tree is None:
            found = False
            for event, elem in etree.iterparse(self.source,
                    events=('start', 'end'), remove_comments=True):
                if event == 'start':
                    if elem.attrib.get('id') == elem_id:
                        found = True
                else:
                    if elem.attrib.get('id') == elem_id:
                        return self.get_info_smart(elem, **kwargs)
                    if not found:
                        elem.clear()
            return None
        # Otherwise do build and use the id_dict to cache elements
        else:
            if len(self.id_dict) == 0:
                self._build_id_cache()
            return self.get_info_smart(self.id_dict[elem_id], **kwargs)

    @_oo_keepstate
    def _parse_tree(self):
        '''Build and store the ElementTree instance for `self.source`'''
        p = etree.XMLParser(remove_comments=True)
        self.tree = etree.parse(self.source, parser=p)

    @_oo_keepstate
    def iterfind(self, path, iterative=True, **kwargs):
        """Parse `source` and yield info on elements with specified local name
        or by specified "XPath". Only local names separated with slashes are
        accepted. An asterisk (`*`) means any element.
        You can specify a single condition in the end, such as:
        "/path/to/element[some_value>1.5]"
        Note: you can do much more powerful filtering using plain Python.
        The path can be absolute or "free". Please don't specify
        namespaces."""
        if self.tree is None:
            self._parse_tree()
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
        if iterative:
            localname = nodes[0].lower()
            found = False
            for ev, elem in etree.iterparse(self.source, events=('start', 'end'),
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

    def __iter__(self):
        for i in self.iterfind(self.default_iter_tag,
                               iterative=self.iterative,
                               retrieve_refs=self.retrieve_refs):
            yield i


class MzIdentMLParser(XMLParserBase):
    file_formatformat = "mzIdentML"
    root_element = "MzIdentML"
    version_info_element = "MzIdentML"
    default_schema = _mzid_schema_defaults
    default_version = "1.1.0"
    default_iter_tag = "SpectrumIdentificationResult"
    structures_to_flatten = {'Fragmentation'}


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
