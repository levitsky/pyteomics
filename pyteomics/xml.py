"""
xml - utilities for XML parsing
===============================

This module is not intended for end users. It implements the abstract classes
for all XML parsers, :py:class:`XML` and :py:class:`IndexedXML`, and some utility functions.

Dependencies
------------

This module requres :py:mod:`lxml` and :py:mod:`numpy`.

--------------------------------------------------------------------------------
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

import re
import socket
from traceback import format_exc
import warnings
from collections import OrderedDict, namedtuple
from itertools import islice
from lxml import etree
import numpy as np

from .auxiliary import FileReader, PyteomicsError, basestring, _file_obj, HierarchicalOffsetIndex
from .auxiliary import unitint, unitfloat, unitstr, cvstr
from .auxiliary import _keepstate_method as _keepstate
from .auxiliary import TaskMappingMixin, IndexedReaderMixin, IndexSavingMixin

try:  # Python 2.7
    from urllib2 import urlopen, URLError
except ImportError:  # Python 3.x
    from urllib.request import urlopen, URLError


def _local_name(element):
    """Strip namespace from the XML element's name"""
    tag = element.tag
    if tag and tag[0] == '{':
        return tag.rpartition('}')[2]
    return tag


def xsd_parser(schema_url):
    """Parse an XSD file from the specified URL into a schema dictionary
    that can be used by :class:`XML` parsers to automatically cast data to
    the appropriate type.

    Parameters
    ----------
    schema_url : str
        The URL to retrieve the schema from

    Returns
    -------
    dict
    """
    ret = {}
    if not (schema_url.startswith('http://') or
            schema_url.startswith('https://') or
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
                anc_name = _local_name(anc)
                while not (
                        (anc_name == 'complexType' and 'name' in anc.attrib) or anc_name == 'element'):
                    anc = anc.getparent()
                    anc_name = _local_name(anc)
                    if anc is None:
                        break
                else:
                    if anc_name == 'complexType':
                        elnames = [x.attrib['name'] for x in
                                   schema_tree.iter()
                                   if x.attrib.get('type', '').split(':')[-1] == anc.attrib['name']]
                    else:
                        elnames = (anc.attrib['name'],)
                    for elname in elnames:
                        tuples.add(
                            (elname, elem.attrib['name']))
        ret[k] = tuples
    ret['lists'] = set(elem.attrib['name'] for elem in schema_tree.xpath(
        '//*[local-name()="element"]') if 'name' in elem.attrib and
        elem.attrib.get('maxOccurs', '1') != '1')
    return ret


class XMLValueConverter(object):
    # Adapted from http://stackoverflow.com/questions/2764269/parsing-an-xsduration-datatype-into-a-python-datetime-timedelta-object
    _duration_parser = re.compile(
        (r'(?P<sign>-?)P(?:(?P<years>\d+\.?\d*)Y)?(?:(?P<months>\d+\.?\d*)M)?(?:(?P<days>\d+\.?\d*)D)?(?:T(?:(?P<hours>\d+\.?\d*)H)?(?:(?P<minutes>\d+\.?\d*)M)?(?:(?P<seconds>\d+\.?\d*)S)?)?'))

    @classmethod
    def duration_str_to_float(cls, s):
        # Not a duration, so pass along
        if not s.startswith('P'):
            try:
                return unitfloat(s, 'duration')
            except ValueError:
                return unitstr(s, 'duration')
        match = cls._duration_parser.search(s)
        if match:
            matchdict = match.groupdict()
            hours = float(matchdict.get('hours', 0) or 0)
            minutes = float(matchdict.get('minutes', 0) or 0)
            seconds = float(matchdict.get('seconds', 0) or 0)
            minutes += hours * 60.
            minutes += (seconds / 60.)
            return unitfloat(minutes, 'minute')
        else:
            return unitstr(s, 'duration')

    @classmethod
    def str_to_bool(cls, s):
        if s.lower() in {'true', '1', 'y'}:
            return True
        if s.lower() in {'false', '0', 'n'}:
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
            'ints': cls.to(unitint), 'floats': cls.to(unitfloat), 'bools': cls.str_to_bool,
            'intlists': lambda x: np.fromstring(x.replace('\n', ' '), dtype=int, sep=' '),
            'floatlists': lambda x: np.fromstring(x.replace('\n', ' '), sep=' '),
            'charlists': list,
            'duration': cls.duration_str_to_float
        }


class _XMLParam(namedtuple("XMLParam", ("name", "value", "type"))):
    '''A holder for semantic parameters used in several common XML formats

    Attributes
    ----------
    name: :class:`~.cvstr`
        The name of the attribute, carrying the accession and unit information
    value: :class:`~.unitfloat`, :class:`~.unitint` or :class:`~.unitstr`
        The value of the parameter
    type: :class:`str`
        The parameter's local XML tag name.
    '''
    __slots__ = ()

    def is_empty(self):
        value = self.value
        return value == "" or value is None


class XML(FileReader):
    """Base class for all format-specific XML parsers. The instances can be used
    as context managers and as iterators.
    """
    # Configurable data
    file_format = 'XML'
    _root_element = None
    _default_schema = {}
    _read_schema = False
    _default_version = 0
    _default_iter_tag = None
    _default_iter_path = None
    _structures_to_flatten = []
    _schema_location_param = 'schemaLocation'
    _default_id_attr = 'id'
    _huge_tree = False
    _retrieve_refs_enabled = None  # only some subclasses implement this
    _iterative = True

    # Configurable plugin logic
    _converters = XMLValueConverter.converters()
    _element_handlers = {}

    # Must be implemented by subclasses
    def _get_info_smart(self, element, **kwargs):
        raise NotImplementedError

    def __init__(self, source, read_schema=None, iterative=None, build_id_cache=False, **kwargs):
        """Create an XML parser object.

        Parameters
        ----------
        source : str or file
            File name or file-like object corresponding to an XML file.
        read_schema : bool, optional
            Defines whether schema file referenced in the file header
            should be used to extract information about value conversion.
            Default is :py:const:`False`.
        iterative : bool, optional
            Defines whether an :py:class:`ElementTree` object should be
            constructed and stored on the instance or if iterative parsing
            should be used instead. Iterative parsing keeps the memory usage
            low for large XML files. Default is :py:const:`True`.
        build_id_cache : bool, optional
            Defines whether a dictionary mapping IDs to XML tree elements
            should be built and stored on the instance. It is used in
            :py:meth:`XML.get_by_id`, e.g. when using
            :py:class:`pyteomics.mzid.MzIdentML` with ``retrieve_refs=True``.
        huge_tree : bool, optional
            This option is passed to the `lxml` parser and defines whether
            security checks for XML tree depth and node size should be disabled.
            Default is :py:const:`False`.
            Enable this option for trusted files to avoid XMLSyntaxError exceptions
            (e.g. `XMLSyntaxError: xmlSAX2Characters: huge text node`).
        """

        super(XML, self).__init__(source, mode='rb', parser_func=self.iterfind, pass_file=False,
                args=(self._default_iter_path or self._default_iter_tag,), kwargs=kwargs)
        if iterative is None:
            iterative = self._iterative
        if iterative:
            self._tree = None
        else:
            self.build_tree()
        if build_id_cache:
            self.build_id_cache()
        else:
            self._id_dict = None

        self.version_info = self._get_version_info()
        if read_schema is not None:
            self._read_schema = read_schema
        self.schema_info = self._get_schema_info(read_schema)

        self._converters_items = self._converters.items()
        self._huge_tree = kwargs.get('huge_tree', self._huge_tree)
        self._retrieve_refs_enabled = kwargs.get('retrieve_refs')

    def __reduce_ex__(self, protocol):
        return self.__class__, (
            self._source_init, self._read_schema, self._tree is None,
            False,
        ), self.__getstate__()

    def __getstate__(self):
        state = super(XML, self).__getstate__()
        state['_huge_tree'] = self._huge_tree
        state['_retrieve_refs_enabled'] = self._retrieve_refs_enabled
        state['_id_dict'] = self._id_dict
        return state

    def __setstate__(self, state):
        super(XML, self).__setstate__(state)
        self._huge_tree = state['_huge_tree']
        self._retrieve_refs_enabled = state['_retrieve_refs_enabled']
        self._id_dict = state['_id_dict']

    @_keepstate
    def _get_version_info(self):
        """
        Provide version information about the XML file.

        Returns
        -------
        out : tuple
            A (version, schema URL) tuple, both elements are strings or None.
        """
        for _, elem in etree.iterparse(
                self._source, events=('start',), remove_comments=True, huge_tree=self._huge_tree):
            if _local_name(elem) == self._root_element:
                return (elem.attrib.get('version'),
                        elem.attrib.get(('{{{}}}'.format(elem.nsmap['xsi'])
                            if 'xsi' in elem.nsmap else '') + self._schema_location_param))

    @_keepstate
    def _get_schema_info(self, read_schema=True):
        """Stores defaults for the schema, tries to retrieve the schema for
        other versions. Keys are: 'floats', 'ints', 'bools', 'lists',
        'intlists', 'floatlists', 'charlists'."""
        if not read_schema:
            return self._default_schema

        version, schema = self.version_info
        if version == self._default_version:
            return self._default_schema

        ret = {}
        try:
            if not schema:
                schema_url = ''
                raise PyteomicsError(
                        'Schema information not found in {}.'.format(self.name))
            schema_url = schema.split()[-1]
            ret = xsd_parser(schema_url)
        except Exception as e:
            if isinstance(e, (URLError, socket.error, socket.timeout)):
                warnings.warn("Can't get the {0.file_format} schema for version "
                "`{1}` from <{2}> at the moment.\n"
                "Using defaults for {0._default_version}.\n"
                "You can disable reading the schema by specifying "
                "`read_schema=False`.".format(self, version, schema_url))
            else:
                warnings.warn("Unknown {0.file_format} version `{1}`.\n"
                    "Attempt to use schema "
                    "information from <{2}> failed.\n"
                    "Exception information:\n{3}\n"
                    "Falling back to defaults for {0._default_version}\n"
                    "NOTE: This is just a warning, probably from a badly-"
                    "generated XML file.\nYou will still most probably get "
                    "decent results.\nLook here for suppressing warnings:\n"
                    "http://docs.python.org/library/warnings.html#"
                    "temporarily-suppressing-warnings\n"
                    "You can also disable reading the schema by specifying "
                    "`read_schema=False`.\n"
                    "If you think this shouldn't have happened, please "
                    "report this to\n"
                    "http://github.com/levitsky/pyteomics/issues\n"
                    "".format(self, version, schema_url, format_exc()))
            ret = self._default_schema
        return ret

    def _handle_param(self, element, **kwargs):
        """Unpacks cvParam and userParam tags into key-value pairs"""
        types = {'int': unitint, 'float': unitfloat, 'string': unitstr}
        attribs = element.attrib
        unit_info = None
        unit_accesssion = None
        if 'unitCvRef' in attribs or 'unitName' in attribs:
            unit_accesssion = attribs.get('unitAccession')
            unit_name = attribs.get('unitName', unit_accesssion)
            unit_info = unit_name
        accession = attribs.get('accession')
        value = attribs.get('value', '')
        try:
            if attribs.get('type') in types:
                value = types[attribs['type']](value, unit_info)
            else:
                value = unitfloat(value, unit_info)
        except ValueError:
            value = unitstr(value, unit_info)

        # return {cvstr(attribs['name'], accession, unit_accesssion): value}
        return _XMLParam(cvstr(attribs['name'], accession, unit_accesssion), value, _local_name(element))

    def _handle_referenceable_param_group(self, param_group_ref, **kwargs):
        raise NotImplementedError()
        return []

    def _find_immediate_params(self, element, **kwargs):
        return element.xpath(
            './*[local-name()="cvParam" or local-name()="userParam" or local-name()="UserParam" or local-name()="referenceableParamGroupRef"]')

    def _insert_param(self, info_dict, param):
        key = param.name
        if key in info_dict:
            if isinstance(info_dict[key], list):
                info_dict[key].append(param.value)
            else:
                info_dict[key] = [info_dict[key], param.value]
        else:
            info_dict[key] = param.value

    def _promote_empty_parameter_to_name(self, info, params):
        empty_values = []
        not_empty_values = []
        for param in params:
            if param.is_empty():
                empty_values.append(param)
            else:
                not_empty_values.append(param)

        if len(empty_values) == 1 and 'name' not in info:
            info['name'] = empty_values[0].name
            return info, not_empty_values
        return info, params

    def _get_info(self, element, **kwargs):
        """Extract info from element's attributes, possibly recursive.
        <cvParam> and <userParam> elements are treated in a special way."""
        try:
            name = kwargs.pop('ename')
        except KeyError:
            name = _local_name(element)
        schema_info = self.schema_info
        if name in {'cvParam', 'userParam', 'UserParam'}:
            return self._handle_param(element, **kwargs)
        elif name == "referenceableParamGroupRef":
            return self._handle_referenceable_param_group(element, **kwargs)

        info = dict(element.attrib)
        # process subelements
        params = []
        if kwargs.get('recursive'):
            for child in element.iterchildren():
                cname = _local_name(child)
                if cname in {'cvParam', 'userParam', 'UserParam'}:
                    newinfo = self._handle_param(child, **kwargs)
                    params.append(newinfo)
                elif cname == "referenceableParamGroupRef":
                    params.extend(self._handle_referenceable_param_group(child, **kwargs))
                else:
                    if cname not in schema_info['lists']:
                        info[cname] = self._get_info_smart(child, ename=cname, **kwargs)
                    else:
                        info.setdefault(cname, []).append(
                            self._get_info_smart(child, ename=cname, **kwargs))
        else:
            # handle the case where we do not want to unpack all children, but
            # *Param tags are considered part of the current entity, semantically
            for child in self._find_immediate_params(element, **kwargs):
                param_or_group = self._handle_param(child, **kwargs)
                if isinstance(param_or_group, list):
                    params.extend(param_or_group)
                else:
                    params.append(param_or_group)

        handler = self._element_handlers.get(name)
        if handler is not None:
            info, params = handler(self, info, params)

        for param in params:
            self._insert_param(info, param)

        # process element text
        if element.text:
            stext = element.text.strip()
            if stext:
                if info:
                    info[name] = stext
                else:
                    return stext

        # convert types
        try:
            for k, v in info.items():
                for t, a in self._converters_items:
                    if t in schema_info and (name, k) in schema_info[t]:
                        info[k] = a(v)
        except ValueError as e:
            message = 'Error when converting types: {}'.format(e.args)
            if not self._read_schema:
                message += '\nTry reading the file with read_schema=True'
            raise PyteomicsError(message)

        # resolve refs
        if kwargs.get('retrieve_refs', self._retrieve_refs_enabled):
            self._retrieve_refs(info, **kwargs)

        # flatten the excessive nesting
        for k, v in dict(info).items():
            if k in self._structures_to_flatten:
                if isinstance(v, list):
                    for vi in v:
                        info.update(vi)
                else:
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
        """Build and store the :py:class:`ElementTree` instance
        for the underlying file"""
        p = etree.XMLParser(remove_comments=True, huge_tree=True)
        self._tree = etree.parse(self._source, parser=p)

    def clear_tree(self):
        """Remove the saved :py:class:`ElementTree`."""
        self._tree = None

    def _retrieve_refs(self, info, **kwargs):
        """Retrieves and embeds the data for each attribute in `info` that
        ends in _ref. Removes the id attribute from `info`.

        This implementation is a stub and must be implemented for each specific
        subclass. It is only called if :attr:`retrieve_refs` """
        raise NotImplementedError(
            ("_retrieve_refs is not implemented for {}. "
             "Do not use `retrieve_refs=True`.").format(
                self.__class__.__name__))

    def iterfind(self, path, **kwargs):
        """Parse the XML and yield info on elements with specified local
        name or by specified "XPath".

        Parameters
        ----------
        path : str
            Element name or XPath-like expression. The path is very close to
            full XPath syntax, but local names should be used for all elements in the path.
            They will be substituted with local-name() checks, up to the (first) predicate.
            The path can be absolute or "free". Please don't specify namespaces.
        **kwargs : passed to :py:meth:`self._get_info_smart`.

        Returns
        -------
        out : iterator
        """
        return Iterfind(self, path, **kwargs)

    @_keepstate
    def _iterfind_impl(self, path, **kwargs):
        """Parse the XML and yield info on elements with specified local
        name or by specified "XPath".

        Parameters
        ----------
        path : str
            Element name or XPath-like expression. The path is very close to
            full XPath syntax, but local names should be used for all elements in the path.
            They will be substituted with local-name() checks, up to the (first) predicate.
            The path can be absolute or "free". Please don't specify namespaces.
        **kwargs : passed to :py:meth:`self._get_info_smart`.

        Returns
        -------
        out : iterator
        """
        try:
            path, tail = re.match(pattern_path, path).groups()
        except AttributeError:
            raise PyteomicsError('Invalid path: ' + path)
        if path[:2] == '//' or path[0] != '/':
            absolute = False
            if path[:2] == '//':
                path = path[2:]
                if path[0] == '/' or '//' in path:
                    raise PyteomicsError("Too many /'s in a row.")
        else:
            absolute = True
            path = path[1:]
        nodes = path.rstrip('/').split('/')
        if not nodes:
            raise PyteomicsError('Invalid path: ' + path)

        if not self._tree:
            if tail:
                if tail[0] == '[':
                    tail = '(.)' + tail
                else:
                    raise PyteomicsError('Cannot parse path tail: ' + tail)
                xpath = etree.XPath(tail)
            localname = nodes[0]
            found = False
            for ev, elem in etree.iterparse(self, events=('start', 'end'), remove_comments=True, huge_tree=self._huge_tree):
                name_lc = _local_name(elem)
                if ev == 'start':
                    if name_lc == localname or localname == '*':
                        found += 1
                else:
                    if name_lc == localname or localname == '*':
                        if (absolute and elem.getparent() is None) or not absolute:
                            for child in get_rel_path(elem, nodes[1:]):
                                if tail:
                                    for elem in xpath(child):
                                        info = self._get_info_smart(elem, **kwargs)
                                        yield info
                                else:
                                    info = self._get_info_smart(child, **kwargs)
                                    yield info
                        if not localname == '*':
                            found -= 1
                    if not found:
                        elem.clear()
        else:
            xpath = ('/' if absolute else '//') + '/'.join(
                    '*[local-name()="{}"]'.format(node) if node != '*' else '*' for node in nodes ) + tail
            for elem in self._tree.xpath(xpath):
                info = self._get_info_smart(elem, **kwargs)
                yield info

    @_keepstate
    def build_id_cache(self):
        """Construct a cache for each element in the document, indexed by id
        attribute"""
        stack = 0
        id_dict = {}
        for event, elem in etree.iterparse(self._source, events=('start', 'end'),
                remove_comments=True, huge_tree=self._huge_tree):
            if event == 'start':
                if 'id' in elem.attrib:
                    stack += 1
            else:
                if 'id' in elem.attrib:
                    stack -= 1
                    id_dict[elem.attrib['id']] = elem
                elif stack == 0:
                    elem.clear()
        self._id_dict = id_dict

    def clear_id_cache(self):
        """Clear the element ID cache"""
        self._id_dict = {}

    def _find_by_id_no_reset(self, elem_id, id_key=None):
        """
        An almost exact copy of :meth:`get_by_id` with the difference that it does
        not reset the file reader's position before iterative parsing.

        Parameters
        ----------
        elem_id : str
            The element id to query for

        Returns
        -------
        lxml.Element
        """
        found = False
        if id_key is None:
            id_key = self._default_id_attr
        for event, elem in etree.iterparse(
                self._source, events=('start', 'end'), remove_comments=True, huge_tree=self._huge_tree):
            if event == 'start':
                if elem.attrib.get(id_key) == elem_id:
                    found = True
            else:
                if elem.attrib.get(id_key) == elem_id:
                    return elem
                if not found:
                    elem.clear()
        raise KeyError(elem_id)

    @_keepstate
    def get_by_id(self, elem_id, **kwargs):
        """Parse the file and return the element with `id` attribute equal
        to `elem_id`. Returns :py:const:`None` if no such element is found.

        Parameters
        ----------
        elem_id : str
            The value of the `id` attribute to match.

        Returns
        -------
        out : :py:class:`dict` or :py:const:`None`
        """
        if not self._id_dict:
            elem = self._find_by_id_no_reset(elem_id)
        else:
            elem = self._id_dict[elem_id]
        return self._get_info_smart(elem, **kwargs)


# XPath emulator tools
pattern_path = re.compile(r'([\w/*]*)(.*)')


def get_rel_path(element, names):
    if not names:
        yield element
    else:
        for child in element.iterchildren():
            if names[0] == '*' or _local_name(child) == names[0]:
                if len(names) == 1:
                    yield child
                else:
                    for gchild in get_rel_path(child, names[1:]):
                        yield gchild


def xpath(tree, path, ns=None):
    """Return the results of XPath query with added namespaces.
    Assumes the ns declaration is on the root element or absent.

    Parameters
    ----------

    tree : ElementTree
    path : str
    ns   : str or None, optional
    """
    if hasattr(tree, 'getroot'):
        root = tree.getroot()
    else:
        root = tree
        while root.getparent() is not None:
            root = root.getparent()
    ns = root.nsmap.get(ns)

    def repl(m):
        s = m.group(1)
        if not ns: return s
        if not s: return 'd:'
        return '/d:'
    new_path = re.sub(r'(\/|^)(?![\*\/])', repl, path)
    n_s = ({'d': ns} if ns else None)
    return tree.xpath(new_path, namespaces=n_s)


def _make_version_info(cls):
    def version_info(source):
        return cls(source).version_info
    version_info.__doc__ = """
    Provide version information about the {0.file_format} file.

    .. note:: This function is provided for backward compatibility only.
        It simply creates an :py:class:`{0.__name__}` instance
        and returns its :py:data:`!version_info` attribute.

    Parameters
    ----------
    source : str or file
        File name or file-like object.

    Returns
    -------
    out : tuple
        A (version, schema URL) tuple, both elements are strings or None.
    """.format(cls)
    return version_info


class ByteCountingXMLScanner(_file_obj):
    """
    Carry out the construction of a byte offset index for `source` XML file
    for each type of tag in :attr:`indexed_tags`.

    Inheris from :py:class:`pyteomics.auxiliary._file_obj` to support the object-oriented
    :py:func:`_keep_state` interface.
    """
    entities = {
        'quot': '"',
        'amp': '&',
        'apos': "'",
        'lt': '<',
        'gt': '>',
    }

    xml_entity_pattern = re.compile(r"&({});".format('|'.join(entities.keys())))

    def __init__(self, source, indexed_tags, block_size=1000000):
        """
        Parameters
        ----------
        indexed_tags : iterable of bytes
            The XML tags (without namespaces) to build indices for.
        block_size : int, optional
            The size of the each chunk or "block" of the file to hold in memory as a
            partitioned string at any given time. Defaults to `1000000`.
        """
        super(ByteCountingXMLScanner, self).__init__(source, 'rb')
        self.indexed_tags = ensure_bytes(indexed_tags)
        self.block_size = block_size

    def _chunk_iterator(self):
        """
        Read a file in large blocks and chunk up each block into parts
        resembling XML tags, yielding each chunk.

        Assumes the file is opened in binary mode.
        """
        f = self.file
        read_size = self.block_size
        delim = b'<'
        buff = f.read(read_size)
        started_with_delim = buff.startswith(delim)
        parts = buff.split(delim)
        tail = parts[-1]
        front = parts[:-1]
        i = 0
        for part in front:
            i += 1
            if part == b"":
                continue
            if i == 1:
                if started_with_delim:
                    yield delim + part
                else:
                    yield part
            else:
                yield delim + part
        running = True
        while running:
            buff = f.read(read_size)
            if not buff:
                running = False
                buff = tail
            else:
                buff = tail + buff
            parts = buff.split(delim)
            tail = parts[-1]
            front = parts[:-1]
            for part in front:
                yield delim + part

    def _generate_offsets(self):
        """
        Iterate over the lines of an XML file where each line contains exactly one tag,
        tracking the byte count for each line. When a line contains a tag whose name matches
        a name in :attr:`indexed_tags`, yield the byte offset, the tag type, and it's attributes.

        Yields
        ------
        offset : int
            The byte offset of a matched tag's opening line
        tag_type : bytes
            The type of tag matched
        attr_dict : dict
            The attributes on the matched tag
        """
        i = 0
        packed = b"|".join(self.indexed_tags)
        pattern = re.compile((r"^\s*<(%s)\s" % packed.decode()).encode())
        attrs = re.compile(br"(\S+)=[\"']([^\"']*)[\"']")
        for line in self._chunk_iterator():
            match = pattern.match(line)
            if match:
                yield i, match.group(1), dict(attrs.findall(line))
            i += len(line)

    def _entity_sub_cb(self, match):
        ent = match.group(1)
        return self.entities[ent]

    def replace_entities(self, key):
        '''Replace XML entities in a string with their character representation

        Uses the minimal mapping of XML entities pre-defined for all XML documents and
        does not attempt to deal with external DTD defined entities. This mapping is found
        in :attr:`entities`.

        Parameters
        ----------
        key : str
            The string to substitute

        Returns
        -------
        str
        '''
        return self.xml_entity_pattern.sub(self._entity_sub_cb, key)

    @_keepstate
    def build_byte_index(self, lookup_id_key_mapping=None):
        """
        Builds a byte offset index for one or more types of tags.

        Parameters
        ----------
        lookup_id_key_mapping : Mapping, optional
            A mapping from tag name to the attribute to look up the identity
            for each entity of that type to be extracted. Defaults to 'id' for
            each type of tag.

        Returns
        -------
        defaultdict(dict)
            Mapping from tag type to dict from identifier to byte offset
        """
        if lookup_id_key_mapping is None:
            lookup_id_key_mapping = {}
        lookup_id_key_mapping = {ensure_bytes_single(key): ensure_bytes_single(value)
            for key, value in lookup_id_key_mapping.items()}

        for name in self.indexed_tags:
            bname = ensure_bytes_single(name)
            lookup_id_key_mapping.setdefault(bname, 'id')
            lookup_id_key_mapping[bname] = ensure_bytes_single(lookup_id_key_mapping[bname])

        indices = HierarchicalOffsetIndex()
        g = self._generate_offsets()
        for offset, offset_type, attrs in g:
            k = attrs[lookup_id_key_mapping[offset_type]].decode('utf-8')
            if '&' in k:
                k = self.replace_entities(k)
            indices[offset_type.decode('utf-8')][k] = offset
        return indices

    @classmethod
    def scan(cls, source, indexed_tags):
        inst = cls(source, indexed_tags)
        return inst.build_byte_index()


class TagSpecificXMLByteIndex(object):
    """
    Encapsulates the construction and querying of a byte offset index
    for a set of XML tags.

    This type mimics an immutable Mapping.

    Attributes
    ----------
    indexed_tags : iterable of bytes
        The tag names to index, not including a namespace
    offsets : defaultdict(OrderedDict(str, int))
        The hierarchy of byte offsets organized ``{"tag_type": {"id": byte_offset}}``
    indexed_tag_keys: dict(str, str)
        A mapping from tag name to unique identifier attribute

    Parameters
    ----------
    index_tags: iterable of bytes
        The tag names to include in the index

    """
    _default_indexed_tags = []
    _default_keys = {}
    _scanner_class = ByteCountingXMLScanner

    def __init__(self, source, indexed_tags=None, keys=None):
        if keys is None:
            keys = self._default_keys.copy()
        if indexed_tags is None:
            indexed_tags = self._default_indexed_tags
        self.indexed_tags = indexed_tags
        self.indexed_tag_keys = keys
        self.source = source
        self.offsets = HierarchicalOffsetIndex()
        self.build_index()

    def __getstate__(self):
        state = {}
        state['indexed_tags'] = self.indexed_tags
        state['indexed_tag_keys'] = self.indexed_tag_keys
        state['offsets'] = self.offsets
        return state

    def __setstate__(self, state):
        self.indexed_tags = state['indexed_tags']
        self.indexed_tag_keys = state['indexed_tag_keys']
        self.offsets = state['offsets']

    def __getitem__(self, key):
        return self.offsets[key]

    def build_index(self):
        """
        Perform the byte offset index building for py:attr:`source`.

        Returns
        -------
        offsets: defaultdict
            The hierarchical offset, stored in offsets
        """
        scanner = self._scanner_class(self.source, self.indexed_tags)
        self.offsets = scanner.build_byte_index(self.indexed_tag_keys)
        return self.offsets

    def items(self):
        return self.offsets.items()

    def keys(self):
        return self.offsets.keys()

    def __iter__(self):
        return iter(self.keys())

    def __len__(self):
        return sum(len(group) for key, group in self.items())

    @classmethod
    def build(cls, source, indexed_tags=None, keys=None):
        indexer = cls(source, indexed_tags, keys)
        return indexer.offsets


def ensure_bytes_single(string):
    if isinstance(string, bytes):
        return string
    try:
        return string.encode('utf-8')
    except (AttributeError, UnicodeEncodeError):
        raise PyteomicsError('{!r} could not be encoded'.format(string))


def ensure_bytes(strings):
    if isinstance(strings, basestring):
        strings = [strings]
    return [ensure_bytes_single(string) for string in strings]


def _flatten_map(hierarchical_map):
    all_records = []
    for key, records in hierarchical_map.items():
        all_records.extend(records.items())

    all_records.sort(key=lambda x: x[1])
    return OrderedDict(all_records)


class IndexedXML(IndexedReaderMixin, XML):
    """Subclass of :py:class:`XML` which uses an index of byte offsets for some
    elements for quick random access.
    """
    _indexed_tags = set()
    _indexed_tag_keys = {}
    _use_index = True

    def __init__(self, source, read_schema=False, iterative=True, build_id_cache=False,
                 use_index=None, *args, **kwargs):
        """Create an indexed XML parser object.

        Parameters
        ----------
        source : str or file
            File name or file-like object corresponding to an XML file.
        read_schema : bool, optional
            Defines whether schema file referenced in the file header
            should be used to extract information about value conversion.
            Default is :py:const:`False`.
        iterative : bool, optional
            Defines whether an :py:class:`ElementTree` object should be
            constructed and stored on the instance or if iterative parsing
            should be used instead. Iterative parsing keeps the memory usage
            low for large XML files. Default is :py:const:`True`.
        use_index : bool, optional
            Defines whether an index of byte offsets needs to be created for
            elements listed in `indexed_tags`.
            This is useful for random access to spectra in mzML or elements of mzIdentML files,
            or for iterative parsing of mzIdentML with ``retrieve_refs=True``.
            If :py:const:`True`, `build_id_cache` is ignored.
            If :py:const:`False`, the object acts exactly like :py:class:`XML`.
            Default is :py:const:`True`.
        indexed_tags : container of bytes, optional
            If `use_index` is :py:const:`True`, elements listed in this parameter
            will be indexed. Empty set by default.
        """
        tags = kwargs.get('indexed_tags')
        tag_index_keys = kwargs.get('indexed_tag_keys')

        if tags is not None:
            self._indexed_tags = tags
        if tag_index_keys is not None:
            self._indexed_tag_keys = tag_index_keys

        if use_index is not None:
            self._use_index = use_index

        if use_index:
            build_id_cache = False
            if self._default_iter_path and self._default_iter_path != self._default_iter_tag:
                warnings.warn('_default_iter_path differs from _default_iter_tag and index is enabled. '
                    '_default_iter_tag will be used in the index, mind the consequences.')
        super(IndexedXML, self).__init__(source, read_schema, iterative, build_id_cache, *args, **kwargs)

        self._offset_index = self.build_byte_index()

    @property
    def default_index(self):
        return self._offset_index[self._default_iter_tag]

    def __reduce_ex__(self, protocol):
        reconstructor, args, state = XML.__reduce_ex__(self, protocol)
        args = args + (False, )
        return reconstructor, args, state

    def __getstate__(self):
        state = super(IndexedXML, self).__getstate__()
        state['_indexed_tags'] = self._indexed_tags
        state['_indexed_tag_keys'] = self._indexed_tag_keys
        state['_use_index'] = self._use_index
        state['_offset_index'] = self._offset_index
        return state

    def __setstate__(self, state):
        super(IndexedXML, self).__setstate__(state)
        self._indexed_tags = state['_indexed_tags']
        self._indexed_tag_keys = state['_indexed_tag_keys']
        self._use_index = state['_use_index']
        self._offset_index = state['_offset_index']

    @_keepstate
    def build_byte_index(self):
        """
        Build up an index of offsets for elements.

        Returns
        -------

        out : TagSpecificXMLByteIndex
        """
        if not self._indexed_tags or not self._use_index:
            return
        return TagSpecificXMLByteIndex.build(
            self._source, self._indexed_tags, self._indexed_tag_keys)

    @_keepstate
    def _find_by_id_reset(self, elem_id, id_key=None):
        return self._find_by_id_no_reset(elem_id, id_key=id_key)

    @_keepstate
    def get_by_id(self, elem_id, id_key=None, element_type=None, **kwargs):
        """
        Retrieve the requested entity by its id. If the entity
        is a spectrum described in the offset index, it will be retrieved
        by immediately seeking to the starting position of the entry, otherwise
        falling back to parsing from the start of the file.

        Parameters
        ----------
        elem_id : str
            The id value of the entity to retrieve.
        id_key : str, optional
            The name of the XML attribute to use for lookup.
            Defaults to :py:attr:`self._default_id_attr`.

        Returns
        -------
        dict
        """
        try:
            index = self._offset_index
            if element_type is None:
                offset, element_type = index.find_no_type(elem_id)
            else:
                offset = index.find(elem_id, element_type)
            self._source.seek(offset)
            if id_key is None:
                id_key = self._indexed_tag_keys.get(element_type)
            elem = self._find_by_id_no_reset(elem_id, id_key=id_key)
        except (KeyError, AttributeError, etree.LxmlError):
            elem = self._find_by_id_reset(elem_id, id_key=id_key)
        data = self._get_info_smart(elem, **kwargs)
        return data

    def __contains__(self, key):
        return key in self._offset_index[self._default_iter_tag]

    def __len__(self):
        return len(self._offset_index[self._default_iter_tag])

    def iterfind(self, path, **kwargs):
        """Parse the XML and yield info on elements with specified local
        name or by specified "XPath".

        Parameters
        ----------
        path : str
            Element name or XPath-like expression. The path is very close to
            full XPath syntax, but local names should be used for all elements in the path.
            They will be substituted with local-name() checks, up to the (first) predicate.
            The path can be absolute or "free". Please don't specify namespaces.
        **kwargs : passed to :py:meth:`self._get_info_smart`.

        Returns
        -------
        out : iterator
        """
        if path in self._indexed_tags and self._use_index:
            return IndexedIterfind(self, path, **kwargs)
        return Iterfind(self, path, **kwargs)

    def __init_subclass__(cls, **kwargs):  # only works on Python 3.x
        super(IndexedXML, cls).__init_subclass__(**kwargs)
        if hasattr(cls, '_build_index'):
            warnings.warn("The method `_build_index` has been renamed to `build_byte_index`.")

class MultiProcessingXML(IndexedXML, TaskMappingMixin):
    """XML reader that feeds indexes to external processes
    for parallel parsing and analysis of XML entries."""

    def _task_map_iterator(self):
        """Returns the :class:`Iteratable` to use when dealing work items onto the input IPC
        queue used by :meth:`map`

        Returns
        -------
        :class:`Iteratable`
        """
        return iter(self._offset_index[self._default_iter_tag])


class IndexSavingXML(IndexSavingMixin, IndexedXML):
    """An extension to the IndexedXML type which
    adds facilities to read and write the byte offset
    index externally.
    """
    _index_class = HierarchicalOffsetIndex

    def _read_byte_offsets(self):
        """Read the byte offset index JSON file at :attr:`_byte_offset_filename`
        and populate :attr:`_offset_index`
        """
        with open(self._byte_offset_filename, 'r') as f:
            index = self._index_class.load(f)
            if index.schema_version is None:
                raise TypeError("Legacy Offset Index!")
            return index


class Iterfind(object):
    def __init__(self, parser, tag_name, **kwargs):
        self.parser = parser
        self.tag_name = tag_name
        self.config = kwargs
        self._iterator = None

    def __repr__(self):
        template = "{self.__class__.__name__}({self.tag_name!r}{config})"
        if self.config:
            config = ", " + repr(self.config)
        else:
            config = ''
        return template.format(self=self, config=config)

    def __iter__(self):
        return self

    def _make_iterator(self):
        return self.parser._iterfind_impl(self.tag_name, **self.config)

    def __next__(self):
        if self._iterator is None:
            self._iterator = self._make_iterator()
        return next(self._iterator)

    def next(self):
        return self.__next__()

    @property
    def is_indexed(self):
        return False

    def reset(self):
        self._iterator = None
        self.parser.reset()

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        self.reset()

    def map(self, *args,**kwargs):
        raise NotImplementedError("This query isn't indexed, it cannot be mapped with multiprocessing")

    def _get_by_index(self, idx):
        self.reset()
        value = next(islice(self, idx, idx + 1))
        return value

    def _get_by_slice(self, slc):
        self.reset()
        value = list(islice(self, slc.start, slc.stop, slc.step))
        return value

    def __getitem__(self, i):
        if isinstance(i, slice):
            return self._get_by_slice(i)
        return self._get_by_index(i)


class IndexedIterfind(TaskMappingMixin, Iterfind):

    def __init__(self, parser, tag_name, **kwargs):
        TaskMappingMixin.__init__(self, **kwargs)
        Iterfind.__init__(self, parser, tag_name, **kwargs)

    def _task_map_iterator(self):
        """Returns the :class:`Iteratable` to use when dealing work items onto the input IPC
        queue used by :meth:`map`

        Returns
        -------
        :class:`Iteratable`
        """
        return iter(self._index)

    @property
    def _offset_index(self):
        return self._index

    @property
    def _index(self):
        return self.parser.index[self.tag_name]

    def _get_reader_for_worker_spec(self):
        return self.parser

    def _yield_from_index(self):
        for key in self._task_map_iterator():
            yield self.parser.get_by_id(key, **self.config)

    def _make_iterator(self):
        if self.is_indexed:
            return self._yield_from_index()
        warnings.warn("Non-indexed iterator created from %r" % (self, ))
        return super(IndexedIterfind, self)._make_iterator()

    @property
    def is_indexed(self):
        if hasattr(self.parser, 'index'):
            if self.parser.index is not None:
                index = self.parser.index
                if isinstance(index, HierarchicalOffsetIndex):
                    return bool(self.tag_name in index and index[self.tag_name])
        return False

    def _get_by_index(self, idx):
        index = self._index
        key = index.from_index(idx)
        return self.parser.get_by_id(key)

    def _get_by_slice(self, slc):
        index = self._index
        keys = index.from_slice(slc)
        return self.parser.get_by_ids(keys)

    def __len__(self):
        index = self._index
        return len(index)
