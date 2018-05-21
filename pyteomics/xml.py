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
import warnings
warnings.formatwarning = lambda msg, *args, **kw: str(msg) + '\n'
import socket
from traceback import format_exc
import operator as op
import ast
import os
import json
import numpy as np
from lxml import etree
from collections import OrderedDict, defaultdict
from .auxiliary import FileReader, PyteomicsError, basestring, _file_obj
from .auxiliary import unitint, unitfloat, unitstr, cvstr
from .auxiliary import _keepstate_method as _keepstate
from .auxiliary import BinaryDataArrayTransformer

try: # Python 2.7
    from urllib2 import urlopen, URLError
except ImportError: # Python 3.x
    from urllib.request import urlopen, URLError


def _local_name(element):
    """Strip namespace from the XML element's name"""
    if element.tag and element.tag[0] == '{':
        return element.tag.rpartition('}')[2]
    return element.tag


class XMLValueConverter(object):
    # Adapted from http://stackoverflow.com/questions/2764269/parsing-an-xsduration-datatype-into-a-python-datetime-timedelta-object
    _duration_parser = re.compile(
        (r'(?P<sign>-?)P(?:(?P<years>\d+\.?\d*)Y)?(?:(?P<months>\d+\.?\d*)M)?(?:(?P<days>\d+\.?\d*)D)?(?:T(?:(?P<hours>\d+\.?\d*)H)?(?:(?P<minutes>\d+\.?\d*)M)?(?:(?P<seconds>\d+\.?\d*)S)?)?'))

    @classmethod
    def duration_str_to_float(cls, s):
        # Not a duration, so pass along unchanged
        if not s.startswith("P"):
            return unitstr(s, "duration")
        match = cls._duration_parser.search(s)
        if match:
            matchdict = match.groupdict()
            hours = float(matchdict.get('hours', 0) or 0)
            minutes = float(matchdict.get('minutes', 0) or 0)
            seconds = float(matchdict.get('seconds', 0) or 0)
            minutes += hours * 60.
            minutes += (seconds / 60.)
            return unitfloat(minutes, "minute")
        else:
            return unitstr(s, "duration")

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


class XML(FileReader):
    """Base class for all format-specific XML parsers. The instances can be used
    as context managers and as iterators.
    """
    # Configurable data
    file_format = 'XML'
    _root_element = None
    _default_schema = {}
    _default_version = 0
    _default_iter_tag = None
    _structures_to_flatten = []
    _schema_location_param = 'schemaLocation'
    _default_id_attr = 'id'
    _huge_tree = False
    _skip_empty_cvparam_values = False
    _retrieve_refs_enabled = None # only some subclasses implement this

    # Configurable plugin logic
    _converters = XMLValueConverter.converters()

    # Must be implemented by subclasses
    def _get_info_smart(self, element, **kwargs):
        raise NotImplementedError

    def __init__(self, source, read_schema=False,
                 iterative=True, build_id_cache=False, **kwargs):
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
        skip_empty_cvparam_values : bool, optional
            .. warning ::
                This parameter affects the format of the produced dictionaries.

            By default, when parsing cvParam elements, "value" attributes with empty values are not
            treated differently from others. When this parameter is set to :py:const:`True`,
            these empty values are flattened. You can enable this to obtain the same output structure
            regardless of the presence of an empty "value". Default is :py:const:`False`.
        """

        super(XML, self).__init__(source, 'rb', self.iterfind, False,
                (self._default_iter_tag,), kwargs)

        if iterative:
            self._tree = None
        else:
            self.build_tree()
        if build_id_cache:
            self.build_id_cache()
        else:
            self._id_dict = None

        self.version_info = self._get_version_info()
        self._read_schema = read_schema
        self.schema_info = self._get_schema_info(read_schema)

        self._converters_items = self._converters.items()
        self._huge_tree = kwargs.get('huge_tree', self._huge_tree)
        self._skip_empty_cvparam_values = kwargs.get('skip_empty_cvparam_values', False)
        self._retrieve_refs_enabled = kwargs.get('retrieve_refs')

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
                        anc_name = _local_name(anc)
                        while not (
                                (anc_name == 'complexType'
                                    and 'name' in anc.attrib)
                                or anc_name == 'element'):
                            anc = anc.getparent()
                            anc_name = _local_name(anc)
                            if anc is None:
                                break
                        else:
                            if anc_name == 'complexType':
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
                    "http://hg.theorchromo.ru/pyteomics/issues\n"
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
            unit_name = attribs.get("unitName", unit_accesssion)
            unit_info = unit_name
        accession = attribs.get("accession")
        if 'value' in attribs and (not self._skip_empty_cvparam_values or
            attribs['value'] != ''):
            try:
                if attribs.get('type') in types:
                    value = types[attribs['type']](attribs['value'], unit_info)
                else:
                    value = unitfloat(attribs['value'], unit_info)
            except ValueError:
                value = unitstr(attribs['value'], unit_info)
            return {cvstr(attribs['name'], accession, unit_accesssion): value}
        else:
            return {'name': cvstr(attribs['name'], accession, unit_accesssion)}

    def _get_info(self, element, **kwargs):
        """Extract info from element's attributes, possibly recursive.
        <cvParam> and <userParam> elements are treated in a special way."""
        try:
            name = kwargs.pop('ename')
        except KeyError:
            name = _local_name(element)
        schema_info = self.schema_info
        if name in {'cvParam', 'userParam'}:
            return self._handle_param(element, **kwargs)

        info = dict(element.attrib)
        # process subelements
        if kwargs.get('recursive'):
            for child in element.iterchildren():
                cname = _local_name(child)
                if cname in {'cvParam', 'userParam', 'UserParam'}:
                    newinfo = self._handle_param(child, **kwargs)
                    if not ('name' in info and 'name' in newinfo):
                        for key in set(info) & set(newinfo):
                            if isinstance(info[key], list):
                                info[key].append(newinfo.pop(key))
                            else:
                                info[key] = [info[key], newinfo.pop(key)]
                        info.update(newinfo)
                    else:
                        if not isinstance(info['name'], list):
                            info['name'] = [info['name']]
                        info['name'].append(newinfo.pop('name'))
                else:
                    if cname not in schema_info['lists']:
                        info[cname] = self._get_info_smart(child, ename=cname, **kwargs)
                    else:
                        info.setdefault(cname, []).append(
                                self._get_info_smart(child, ename=cname, **kwargs))

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

    @_keepstate
    def iterfind(self, path, **kwargs):
        """Parse the XML and yield info on elements with specified local
        name or by specified "XPath".

        Parameters
        ----------
        path : str
            Element name or XPath-like expression. Only local names separated
            with slashes are accepted. An asterisk (`*`) means any element.
            You can specify a single condition in the end, such as:
            ``"/path/to/element[some_value>1.5]"``
            Note: you can do much more powerful filtering using plain Python.
            The path can be absolute or "free". Please don't specify
            namespaces.
        **kwargs : passed to :py:meth:`self._get_info_smart`.

        Returns
        -------
        out : iterator
        """
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
        if not self._tree:
            localname = nodes[0].lower()
            found = False
            for ev, elem in etree.iterparse(self, events=('start', 'end'),
                    remove_comments=True, huge_tree=self._huge_tree):
                name_lc = _local_name(elem).lower()
                if ev == 'start':
                    if name_lc == localname or localname == '*':
                        found += True
                else:
                    if name_lc == localname or localname == '*':
                        if (absolute and elem.getparent() is None) or not absolute:
                            for child in get_rel_path(elem, nodes[1:]):
                                info = self._get_info_smart(child, **kwargs)
                                if cond is None or satisfied(info, cond):
                                    yield info
                        if not localname == '*':
                            found -= 1
                    if not found:
                        elem.clear()
        else:
            xpath = ('/' if absolute else '//') + '/'.join(
                    '*[local-name()="{}"]'.format(node) for node in nodes)
            for elem in self._tree.xpath(xpath):
                info = self._get_info_smart(elem, **kwargs)
                if cond is None or satisfied(info, cond):
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
        elem = None
        if not self._id_dict:
            elem = self._find_by_id_no_reset(elem_id)
        elif elem_id in self._id_dict:
            elem = self._id_dict[elem_id]
        if elem is not None:
            return self._get_info_smart(elem, **kwargs)

# XPath emulator tools
pattern_path = re.compile(r'([\w/*]*)(\[(\w+[<>=]{1,2}[^\]]+)\])?')
pattern_cond = re.compile(r'^\s*(\w+)\s*([<>=]{,2})\s*([^\]]+)$')

def get_rel_path(element, names):
    if not names:
        yield element
    else:
        for child in element.iterchildren():
            if names[0] == '*' or _local_name(child).lower() == names[0].lower():
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
    version_info.__doc__ =  """
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


class ByteEncodingOrderedDict(OrderedDict):

    def __getitem__(self, key):
        try:
            return super(ByteEncodingOrderedDict, self).__getitem__(key)
        except KeyError:
            key = ensure_bytes_single(key)
            return super(ByteEncodingOrderedDict, self).__getitem__(key)

    def __setitem__(self, key, value):
        key = ensure_bytes_single(key)
        return super(ByteEncodingOrderedDict, self).__setitem__(key, value)


class ByteCountingXMLScanner(_file_obj):
    """
    Carry out the construction of a byte offset index for `source` XML file
    for each type of tag in :attr:`indexed_tags`.

    Inheris from :py:class:`pyteomics.auxiliary._file_obj` to support the object-oriented
    :py:func:`_keep_state` interface.
    """
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
            if len(buff) == 0:
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
        attrs = re.compile(br"(\S+)=[\"']([^\"']+)[\"']")
        for line in self._chunk_iterator():
            match = pattern.match(line)
            if match:
                yield i, match.group(1), dict(attrs.findall(line))
            i += len(line)

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
        defaultdict(ByteEncodingOrderedDict)
            Mapping from tag type to ByteEncodingOrderedDict from identifier to byte offset
        """
        if lookup_id_key_mapping is None:
            lookup_id_key_mapping = {}

        for name in self.indexed_tags:
            lookup_id_key_mapping.setdefault(name, "id")
            lookup_id_key_mapping[name] = ensure_bytes_single(lookup_id_key_mapping[name])

        indices = defaultdict(ByteEncodingOrderedDict)
        g = self._generate_offsets()
        for offset, offset_type, attrs in g:
            indices[offset_type][attrs[lookup_id_key_mapping[offset_type]]] = offset
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
        self.offsets = defaultdict(ByteEncodingOrderedDict)
        self.build_index()

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


class FlatTagSpecificXMLByteIndex(TagSpecificXMLByteIndex):
    """
    An alternative interface on top of :py:class:`TagSpecificXMLByteIndex` that assumes
    that identifiers across different tags are globally unique, as in MzIdentML.

    Attributes
    ----------
    offsets : ByteEncodingOrderedDict
        The mapping between ids and byte offsets.
    """
    def build_index(self):
        hierarchical_index = super(FlatTagSpecificXMLByteIndex, self).build_index()
        flat_index = []

        for tag_type in hierarchical_index.values():
            flat_index.extend(tag_type.items())

        flat_index.sort(key=lambda x: x[1])
        self.offsets = ByteEncodingOrderedDict(flat_index)
        return self.offsets

    def __len__(self):
        return len(self.offsets)


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
    return ByteEncodingOrderedDict(all_records)


class IndexedXML(XML):
    """Subclass of :py:class:`XML` which uses an index of byte offsets for some
    elements for quick random access.
    """
    _indexed_tags = set()
    _indexed_tag_keys = {}

    def __init__(self, *args, **kwargs):
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
        use_index = kwargs.get('use_index', True)

        if tags is not None:
            self._indexed_tags = tags
        if tag_index_keys is not None:
            self._indexed_tag_keys = tag_index_keys

        self._use_index = use_index

        self._indexed_tags = ensure_bytes(self._indexed_tags)
        self._indexed_tag_keys = {
            ensure_bytes_single(k): ensure_bytes_single(v)
            for k, v in self._indexed_tag_keys.items()
        }

        if use_index:
            kwargs['build_id_cache'] = False
        super(IndexedXML, self).__init__(*args, **kwargs)
        self._offset_index = ByteEncodingOrderedDict()
        self._build_index()

    @_keepstate
    def _build_index(self):
        """
        Build up a `dict` of `dict` of offsets for elements. Calls :func:`find_index_list`
        on :attr:`_source` and assigns the return value to :attr:`_offset_index`
        """
        if not self._indexed_tags or not self._use_index:
            return
        self._offset_index = FlatTagSpecificXMLByteIndex(
            self._source, self._indexed_tags, self._indexed_tag_keys)

    @_keepstate
    def _find_by_id_reset(self, elem_id, id_key=None):
        return self._find_by_id_no_reset(elem_id, id_key=id_key)

    @_keepstate
    def get_by_id(self, elem_id, id_key=None, **kwargs):
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
            offset = index[elem_id]
            self._source.seek(offset)
            elem = self._find_by_id_no_reset(elem_id, id_key=id_key)
        except (KeyError, etree.LxmlError):
            elem = self._find_by_id_reset(elem_id, id_key=id_key)
        data = self._get_info_smart(elem, **kwargs)
        return data

    def __getitem__(self, elem_id):
        return self.get_by_id(elem_id)


def save_byte_index(index, fp):
    """Write the byte offset index to the provided
    file

    Parameters
    ----------
    index : ByteEncodingOrderedDict
        The byte offset index to be saved
    fp : file
        The file to write the index to

    Returns
    -------
    file
    """
    encoded_index = dict()
    for key, offset in index.items():
        encoded_index[key.decode("utf8")] = offset
    json.dump(encoded_index, fp)
    return fp


def load_byte_index(fp):
    """Read a byte offset index from a file

    Parameters
    ----------
    fp : file
        The file to read the index from

    Returns
    -------
    ByteEncodingOrderedDict
    """
    data = json.load(fp)
    index = ByteEncodingOrderedDict()
    for key, value in sorted(data.items(), key=lambda x: x[1]):
        index[key] = value
    return index


class PrebuiltOffsetIndex(FlatTagSpecificXMLByteIndex):
    """An Offset Index class which just holds offsets
    and performs no extra scanning effort.

    Attributes
    ----------
    offsets : ByteEncodingOrderedDict
    """

    def __init__(self, offsets):
        self.offsets = offsets


class IndexSavingXML(IndexedXML):
    """An extension to the IndexedXML type which
    adds facilities to read and write the byte offset
    index externally.
    """

    _save_byte_index_to_file = staticmethod(save_byte_index)
    _load_byte_index_from_file = staticmethod(load_byte_index)

    @property
    def _byte_offset_filename(self):
        path = self._source.name
        byte_offset_filename = os.path.splitext(path)[0] + '-byte-offsets.json'
        return byte_offset_filename

    def _check_has_byte_offset_file(self):
        """Check if the file at :attr:`_byte_offset_filename` exists

        Returns
        -------
        bool
            Whether the file exists
        """
        path = self._byte_offset_filename
        return os.path.exists(path)

    def _read_byte_offsets(self):
        """Read the byte offset index JSON file at :attr:`_byte_offset_filename`
        and populate :attr:`_offset_index`
        """
        with open(self._byte_offset_filename, 'r') as f:
            index = PrebuiltOffsetIndex(self._load_byte_index_from_file(f))
            self._offset_index = index

    def write_byte_offsets(self):
        """Write the byte offsets in :attr:`_offset_index` to the file
        at :attr:`_byte_offset_filename`
        """
        with open(self._byte_offset_filename, 'w') as f:
            self._save_byte_index_to_file(self._offset_index, f)

    @_keepstate
    def _build_index(self):
        """Build the byte offset index by either reading these offsets
        from the file at :attr:`_byte_offset_filename`, or falling back
        to the method used by :class:`IndexedXML` if this operation fails
        due to an IOError
        """
        try:
            self._read_byte_offsets()
        except (IOError, AttributeError):
            super(IndexSavingXML, self)._build_index()

    @classmethod
    def prebuild_byte_offset_file(cls, path):
        """Construct a new XML reader, build its byte offset index and
        write it to file

        Parameters
        ----------
        path : str
            The path to the file to parse
        """
        with cls(path, use_index=True) as inst:
            inst.write_byte_offsets()

class ArrayConversionMixin(BinaryDataArrayTransformer):
    _dtype_dict = {}
    _array_keys = ['m/z array', 'intensity array']

    def __init__(self, *args, **kwargs):
        self._dtype_dict = {None: None}
        dtype = kwargs.pop('dtype', None)
        if isinstance(dtype, dict):
            self._dtype_dict.update(dtype)
        elif dtype:
            self._dtype_dict = {k: dtype for k in self._array_keys}
            self._dtype_dict[None] = dtype
        super(ArrayConversionMixin, self).__init__(*args, **kwargs)

    def _convert_array(self, k, array):
        dtype = self._dtype_dict.get(k)
        if dtype is not None:
            return array.astype(dtype)
        return array

    def _finalize_record_conversion(self, array, record):
        key = record.key
        return self._convert_array(key, array)
