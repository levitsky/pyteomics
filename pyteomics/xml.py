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
warnings.formatwarning = lambda msg, *args: str(msg) + '\n'
import socket
from functools import wraps
from traceback import format_exc
import operator as op
import ast
import numpy as np
from lxml import etree
import re
from collections import OrderedDict, defaultdict
from .auxiliary import FileReader, PyteomicsError, basestring, _file_obj
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
    """Decorator for :py:class:`XML` methods to help keep the position
    in the underlying file.
    """
    @wraps(func)
    def wrapped(self, *args, **kwargs):
        position = self.tell()
        self.seek(0)
        try:
            return func(self, *args, **kwargs)
        finally:
            self.seek(position)
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

    # Configurable plugin logic
    _converters = XMLValueConverter.converters()

    # Must be implemented by subclasses
    def _get_info_smart(self, element, **kwargs):
        raise NotImplementedError

    def __init__(self, source, read_schema=True,
            iterative=True, build_id_cache=False, **kwargs):
        """Create an XML parser object.

        Parameters
        ----------
        source : str or file
            File name or file-like object corresponding to an XML file.
        read_schema : bool, optional
            Defines whether schema file referenced in the file header
            should be used to extract information about value conversion.
            Default is :py:const:`True`.
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
        self.schema_info = self._get_schema_info(read_schema)

    @_keepstate
    def _get_version_info(self):
        """
        Provide version information about the XML file.

        Returns
        -------
        out : tuple
            A (version, schema URL) tuple, both elements are strings or None.
        """
        vinfo = None
        for _, elem in etree.iterparse(
                self._source, events=('start',), remove_comments=True):
            if _local_name(elem) == self._root_element:
                return (elem.attrib.get('version'),
                        elem.attrib.get(('{{{}}}'.format(elem.nsmap['xsi'])
                            if 'xsi' in elem.nsmap else '') + 'schemaLocation'))

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
        if 'value' in element.attrib:
            try:
                value = float(element.attrib['value'])
            except ValueError:
                value = element.attrib['value']
            return {element.attrib['name']: value}
        else:
            return {'name': element.attrib['name']}

    def _get_info(self, element, **kwargs):
        """Extract info from element's attributes, possibly recursive.
        <cvParam> and <userParam> elements are treated in a special way."""
        name = _local_name(element)
        schema_info = self.schema_info
        if name in {'cvParam', 'userParam'}:
            return self._handle_param(element)

        info = dict(element.attrib)
        # process subelements
        if kwargs.get('recursive'):
            for child in element.iterchildren():
                cname = _local_name(child)
                if cname in {'cvParam', 'userParam'}:
                    newinfo = self._handle_param(child, **kwargs)
                    if not ('name' in info and 'name' in newinfo):
                        info.update(newinfo)
                    else:
                        if not isinstance(info['name'], list):
                            info['name'] = [info['name']]
                        info['name'].append(newinfo.pop('name'))
                else:
                    if cname not in schema_info['lists']:
                        info[cname] = self._get_info_smart(child, **kwargs)
                    else:
                        info.setdefault(cname, []).append(
                                self._get_info_smart(child, **kwargs))

        # process element text
        if element.text and element.text.strip():
            stext = element.text.strip()
            if stext:
                if info:
                    info[name] = stext
                else:
                    return stext

        # convert types
        converters = self._converters
        for k, v in info.items():
            for t, a in converters.items():
                if (_local_name(element), k) in schema_info[t]:
                    info[k] = a(v)

        # resolve refs
        if kwargs.get('retrieve_refs'):
            self._retrieve_refs(info, **kwargs)

        # flatten the excessive nesting
        for k, v in dict(info).items():
            if k in self._structures_to_flatten:
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
                    remove_comments=True, huge_tree=True):
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
        self._id_dict = id_dict

    def clear_id_cache(self):
        """Clear the element ID cache"""
        self._id_dict = {}

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
            found = False
            for event, elem in etree.iterparse(self._source,
                    events=('start', 'end'), remove_comments=True):
                if event == 'start':
                    if elem.attrib.get('id') == elem_id:
                        found = True
                else:
                    if elem.attrib.get('id') == elem_id:
                        return self._get_info_smart(elem, **kwargs)
                    if not found:
                        elem.clear()
            return None
        else:
            try:
                return self._get_info_smart(self._id_dict[elem_id], **kwargs)
            except KeyError:
                return None

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
    new_path = re.sub('(\/|^)(?![\*\/])', repl, path)
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
        parts = buff.split(delim)
        tail = parts[-1]
        front = parts[:-1]
        for part in front:
            if part == b"":
                continue
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
        pattern = re.compile(br"^\s*<(%s)\s" % packed)
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

    Parameters
    ----------
    index_tags: iterable of bytes
        The tag names to include in the index

    """
    _default_indexed_tags = []
    _scanner_class = ByteCountingXMLScanner

    def __init__(self, source, indexed_tags=None):
        if indexed_tags is None:
            indexed_tags = self._default_indexed_tags
        self.indexed_tags = indexed_tags
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
        self.offsets = self._scanner_class.scan(self.source, self.indexed_tags)
        return self.offsets

    def items(self):
        return self.offsets.items()

    def keys(self):
        return self.offsets.keys()

    def __iter__(self):
        return iter(self.keys())


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


def ensure_bytes_single(string):
    if isinstance(string, bytes):
        return string
    try:
        return string.encode('utf-8')
    except (AttributeError, UnicodeEncodeError):
        raise PyteomicsError('%{!r} could not be encoded'.format(string))


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

    def __init__(self, *args, **kwargs):
        """Create an XML parser object.

        Parameters
        ----------
        source : str or file
            File name or file-like object corresponding to an XML file.
        read_schema : bool, optional
            Defines whether schema file referenced in the file header
            should be used to extract information about value conversion.
            Default is :py:const:`True`.
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
            If :py:const:`True`, `build_id_cache`
            is ignored. If :py:const:`False`, the object acts exactly like :py:class:`XML`.
            Default is :py:const:`True`.
        indexed_tags : container of bytes, optional
            If `use_index` is :py:const:`True`, elements listed in this parameter
            will be indexed. Empty set by default.
        """
        tags = kwargs.get('indexed_tags')
        use_index = kwargs.get('use_index', True)

        if tags is not None:
            self._indexed_tags = (tags)

        self._use_index = use_index

        self._indexed_tags = ensure_bytes(self._indexed_tags)

        if use_index: kwargs['build_id_cache'] = False
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
        self._offset_index = FlatTagSpecificXMLByteIndex(self._source, self._indexed_tags)

    def _find_by_id_no_reset(self, elem_id):
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
        for event, elem in etree.iterparse(self._source, events=('start', 'end'), remove_comments=True):
            if event == 'start':
                if elem.attrib.get('id') == elem_id:
                    found = True
            else:
                if elem.attrib.get('id') == elem_id:
                    return elem
                if not found:
                    elem.clear()

    @_keepstate
    def _find_by_id_reset(self, elem_id):
        return self._find_by_id_no_reset(elem_id)

    def get_by_id(self, elem_id, **kwargs):
        """
        Retrieve the requested entity by its id. If the entity
        is a spectrum described in the offset index, it will be retrieved
        by immediately seeking to the starting position of the entry, otherwise
        falling back to parsing from the start of the file.

        Parameters
        ----------
        elem_id : str
            The id value of the entity to retrieve.

        Returns
        -------
        dict
        """
        try:
            index = self._offset_index
            offset = index[elem_id]
            self._source.seek(offset)
            elem = self._find_by_id_no_reset(elem_id)
            data = self._get_info_smart(elem, **kwargs)
            return data
        except (KeyError, etree.LxmlError):
            elem = self._find_by_id_reset(elem_id)
            data = self._get_info_smart(elem, **kwargs)
            return data

    def __getitem__(self, elem_id):
        return self.get_by_id(elem_id)



_mzid_schema_defaults = {
            'ints': {('DBSequence', 'length'),
                     ('IonType', 'charge'),
                     ('BibliographicReference', 'year'),
                     ('SubstitutionModification', 'location'),
                     ('PeptideEvidence', 'end'),
                     ('Enzyme', 'missedCleavages'),
                     ('PeptideEvidence', 'start'),
                     ('Modification', 'location'),
                     ('SpectrumIdentificationItem', 'rank'),
                     ('SpectrumIdentificationItem', 'chargeState'),
                     ('SearchDatabase', 'numDatabaseSequences')},
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

_tandem_schema_defaults = {'ints': {
        ('group', 'z'), ('aa', 'at')} | {('domain', k) for k in [
            'missed_cleavages', 'start', 'end', 'y_ions', 'b_ions',
            'a_ions', 'x_ions', 'c_ions', 'z_ions']},

            'floats': {('group', k) for k in [
                'fI', 'sumI', 'maxI', 'mh', 'expect', 'rt']} | {
                   ('domain', k) for k in [
                       'expect', 'hyperscore', 'b_score', 'y_score',
                       'a_score', 'x_score', 'c_score', 'z_score',
                       'nextscore', 'delta', 'mh']} | {
                   ('protein', 'expect'), ('protein', 'sumI'),
                   ('aa', 'modified')},

            'bools': set(),
            'lists': {'group', 'trace', 'attribute', 'protein', 'aa', 'note'},
            'floatlists': {('values', 'values')},
            'intlists': set(), 'charlists': set()}

_mzml_schema_defaults = {'ints': {
    ('spectrum', 'index'),
     ('instrumentConfigurationList', 'count'),
     ('binaryDataArray', 'encodedLength'),
     ('cvList', 'count'),
     ('binaryDataArray', 'arrayLength'),
     ('scanWindowList', 'count'),
     ('componentList', 'count'),
     ('sourceFileList', 'count'),
     ('productList', 'count'),
     ('referenceableParamGroupList', 'count'),
     ('scanList', 'count'),
     ('spectrum', 'defaultArrayLength'),
     ('dataProcessingList', 'count'),
     ('sourceFileRefList', 'count'),
     ('scanSettingsList', 'count'),
     ('selectedIonList', 'count'),
     ('chromatogram', 'defaultArrayLength'),
     ('precursorList', 'count'),
     ('chromatogram', 'index'),
     ('processingMethod', 'order'),
     ('targetList', 'count'),
     ('sampleList', 'count'),
     ('softwareList', 'count'),
     ('binaryDataArrayList', 'count'),
     ('spectrumList', 'count'),
     ('chromatogramList', 'count')},
        'floats': {},
        'bools': {},
        'lists': {'scan', 'spectrum', 'sample', 'cv', 'dataProcessing',
            'cvParam', 'source', 'userParam', 'detector', 'product',
            'referenceableParamGroupRef', 'selectedIon', 'sourceFileRef',
            'binaryDataArray', 'analyzer', 'scanSettings',
            'instrumentConfiguration', 'chromatogram', 'target',
            'processingMethod', 'precursor', 'sourceFile',
            'referenceableParamGroup', 'contact', 'scanWindow', 'software'},
        'intlists': {},
        'floatlists': {},
        'charlists': {}}

_pepxml_schema_defaults = {'ints':
    {('xpressratio_summary', 'xpress_light'),
     ('distribution_point', 'obs_5_distr'),
     ('distribution_point', 'obs_2_distr'),
     ('enzymatic_search_constraint', 'max_num_internal_cleavages'),
     ('asapratio_lc_heavypeak', 'right_valley'),
     ('libra_summary', 'output_type'),
     ('distribution_point', 'obs_7_distr'),
     ('spectrum_query', 'index'),
     ('data_filter', 'number'),
     ('roc_data_point', 'num_incorr'),
     ('search_hit', 'num_tol_term'),
     ('search_hit', 'num_missed_cleavages'),
     ('asapratio_lc_lightpeak', 'right_valley'),
     ('libra_summary', 'normalization'),
     ('specificity', 'min_spacing'),
     ('database_refresh_timestamp', 'min_num_enz_term'),
     ('enzymatic_search_constraint', 'min_number_termini'),
     ('xpressratio_result', 'light_lastscan'),
     ('distribution_point', 'obs_3_distr'),
     ('spectrum_query', 'end_scan'),
     ('analysis_result', 'id'),
     ('search_database', 'size_in_db_entries'),
     ('search_hit', 'hit_rank'),
     ('alternative_protein', 'num_tol_term'),
     ('search_hit', 'num_tot_proteins'),
     ('asapratio_summary', 'elution'),
     ('search_hit', 'tot_num_ions'),
     ('error_point', 'num_incorr'),
     ('mixture_model', 'precursor_ion_charge'),
     ('roc_data_point', 'num_corr'),
     ('search_hit', 'num_matched_ions'),
     ('dataset_derivation', 'generation_no'),
     ('xpressratio_result', 'heavy_firstscan'),
     ('xpressratio_result', 'heavy_lastscan'),
     ('error_point', 'num_corr'),
     ('spectrum_query', 'assumed_charge'),
     ('analysis_timestamp', 'id'),
     ('xpressratio_result', 'light_firstscan'),
     ('distribution_point', 'obs_4_distr'),
     ('asapratio_lc_heavypeak', 'left_valley'),
     ('fragment_masses', 'channel'),
     ('distribution_point', 'obs_6_distr'),
     ('affected_channel', 'channel'),
     ('search_result', 'search_id'),
     ('contributing_channel', 'channel'),
     ('asapratio_lc_lightpeak', 'left_valley'),
     ('asapratio_peptide_data', 'area_flag'),
     ('search_database', 'size_of_residues'),
     ('asapratio_peptide_data', 'cidIndex'),
     ('mixture_model', 'num_iterations'),
     ('mod_aminoacid_mass', 'position'),
     ('spectrum_query', 'start_scan'),
     ('asapratio_summary', 'area_flag'),
     ('mixture_model', 'tot_num_spectra'),
     ('search_summary', 'search_id'),
     ('xpressratio_timestamp', 'xpress_light'),
     ('distribution_point', 'obs_1_distr'),
     ('intensity', 'channel'),
     ('asapratio_contribution', 'charge'),
     ('libra_summary', 'centroiding_preference')},
    'floats':
    {('asapratio_contribution', 'error'),
     ('asapratio_lc_heavypeak', 'area_error'),
     ('modification_info', 'mod_nterm_mass'),
     ('distribution_point', 'model_4_neg_distr'),
     ('distribution_point', 'model_5_pos_distr'),
     ('spectrum_query', 'precursor_neutral_mass'),
     ('asapratio_lc_heavypeak', 'time_width'),
     ('xpressratio_summary', 'masstol'),
     ('affected_channel', 'correction'),
     ('distribution_point', 'model_7_neg_distr'),
     ('error_point', 'error'),
     ('intensity', 'target_mass'),
     ('roc_data_point', 'sensitivity'),
     ('distribution_point', 'model_4_pos_distr'),
     ('distribution_point', 'model_2_neg_distr'),
     ('distribution_point', 'model_3_pos_distr'),
     ('mixture_model', 'prior_probability'),
     ('roc_data_point', 'error'),
     ('intensity', 'normalized'),
     ('modification_info', 'mod_cterm_mass'),
     ('asapratio_lc_lightpeak', 'area_error'),
     ('distribution_point', 'fvalue'),
     ('distribution_point', 'model_1_neg_distr'),
     ('peptideprophet_summary', 'min_prob'),
     ('asapratio_result', 'mean'),
     ('point', 'pos_dens'),
     ('fragment_masses', 'mz'),
     ('mod_aminoacid_mass', 'mass'),
     ('distribution_point', 'model_6_neg_distr'),
     ('asapratio_lc_lightpeak', 'time_width'),
     ('asapratio_result', 'heavy2light_error'),
     ('peptideprophet_result', 'probability'),
     ('error_point', 'min_prob'),
     ('peptideprophet_summary', 'est_tot_num_correct'),
     ('roc_data_point', 'min_prob'),
     ('asapratio_result', 'heavy2light_mean'),
     ('distribution_point', 'model_5_neg_distr'),
     ('mixturemodel', 'neg_bandwidth'),
     ('asapratio_result', 'error'),
     ('xpressratio_result', 'light_mass'),
     ('point', 'neg_dens'),
     ('asapratio_lc_lightpeak', 'area'),
     ('distribution_point', 'model_1_pos_distr'),
     ('xpressratio_result', 'mass_tol'),
     ('mixturemodel', 'pos_bandwidth'),
     ('xpressratio_result', 'light_area'),
     ('asapratio_peptide_data', 'heavy_mass'),
     ('distribution_point', 'model_2_pos_distr'),
     ('search_hit', 'calc_neutral_pep_mass'),
     ('intensity', 'absolute'),
     ('asapratio_peptide_data', 'light_mass'),
     ('distribution_point', 'model_3_neg_distr'),
     ('aminoacid_modification', 'mass'),
     ('asapratio_lc_heavypeak', 'time'),
     ('asapratio_lc_lightpeak', 'time'),
     ('asapratio_lc_lightpeak', 'background'),
     ('mixture_model', 'est_tot_correct'),
     ('point', 'value'),
     ('asapratio_lc_heavypeak', 'background'),
     ('terminal_modification', 'mass'),
     ('fragment_masses', 'offset'),
     ('xpressratio_result', 'heavy_mass'),
     ('search_hit', 'protein_mw'),
     ('libra_summary', 'mass_tolerance'),
     ('spectrum_query', 'retention_time_sec'),
     ('distribution_point', 'model_7_pos_distr'),
     ('asapratio_lc_heavypeak', 'area'),
     ('alternative_protein', 'protein_mw'),
     ('asapratio_contribution', 'ratio'),
     ('xpressratio_result', 'heavy_area'),
     ('distribution_point', 'model_6_pos_distr')},
    'bools':
    {('sample_enzyme', 'independent'),
     ('intensity', 'reject'),
     ('libra_result', 'is_rejected')},
    'intlists': set(),
    'floatlists': set(),
    'charlists': set(),
    'lists': {'point', 'aminoacid_modification', 'msms_run_summary',
            'mixturemodel', 'search_hit', 'mixturemodel_distribution',
            'sequence_search_constraint', 'specificity', 'alternative_protein',
            'analysis_result', 'data_filter', 'fragment_masses', 'error_point',
            'parameter', 'spectrum_query', 'search_result', 'affected_channel',
            'analysis_summary', 'roc_data_point', 'distribution_point',
            'search_summary', 'mod_aminoacid_mass', 'search_score', 'intensity',
            'analysis_timestamp', 'mixture_model', 'terminal_modification',
            'contributing_channel', 'inputfile'}}

