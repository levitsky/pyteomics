"""
auxiliary - common functions and objects
========================================

Math
----

  :py:func:`linear_regression` - a wrapper for numpy linear regression

Project infrastructure
----------------------

  :py:class:`PyteomicsError` - a pyteomics-specific exception

Helpers
-------

  :py:class:`Charge` - a subclass of :py:class:`int` for charge states

  :py:class:`ChargeList` - a subclass of :py:class:`list` for lists of charges.

  :py:func:`print_tree` - display the structure of a complex nested
  :py:class:`dict`.

  :py:func:`memoize` - makes a
  `memoization <http://stackoverflow.com/a/1988826/1258041>`_
  `function decorator <http://stackoverflow.com/a/1594484/1258041>`_

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

from __future__ import print_function
import numpy as np
from functools import wraps
from lxml import etree
from traceback import format_exc
import re
import operator
try: # Python 2.7
    from urllib2 import urlopen, URLError
except ImportError: # Python 3.x
    from urllib.request import urlopen, URLError
import sys
from contextlib import contextmanager
import socket
import warnings
warnings.formatwarning = lambda msg, *args: str(msg) + '\n'

class PyteomicsError(Exception):
    """Exception raised for errors in Pyteomics library.

    Attributes
    ----------
    message : str
        Error message.
    """

    def __init__(self, msg):
        self.message = msg

    def __str__(self):
        return "Pyteomics error, message: %s" % (repr(self.message),)

def linear_regression(x, y, a=None, b=None):
    """Calculate coefficients of a linear regression y = a * x + b.

    Parameters
    ----------
    x, y : array_like of float
    a : float, optional
        If specified then the slope coefficient is fixed and equals a.
    b : float, optional
        If specified then the free term is fixed and equals b.

    Returns
    -------
    out : 4-tuple of float
        The structure is (a, b, r, stderr), where
        a -- slope coefficient,
        b -- free term,
        r -- Peason correlation coefficient,
        stderr -- standard deviation.
    """

    if not isinstance(x, np.ndarray):
        x = np.array(x)
    if not isinstance(y, np.ndarray):
        y = np.array(y)

    if (a is not None and b is None):
        b = (y - a * x).mean()
    elif (a is not None and b is not None):
        pass
    else:
        a, b = np.polyfit(x, y, 1)

    r = np.corrcoef(x, y)[0, 1]
    stderr = (y - a * x - b).std()

    return a, b, r, stderr

class Charge(int):
    """A subclass of :py:class:`int`. Can be constructed from strings in "N+"
    or "N-" format, and the string representation of a :py:class:`Charge` is
    also in that format.
    """
    def __new__(cls, *args, **kwargs):
        try:
            return super(Charge, cls).__new__(cls, *args)
        except ValueError as e:
            if isinstance(args[0], str):
                try:
                    num, sign = re.match(r'^(\d+)(\+|-)$', args[0]).groups()
                    return super(Charge, cls).__new__(cls,
                        sign + num, *args[1:], **kwargs)
                except:
                    pass
            raise PyteomicsError(*e.args)

    def __str__(self):
        return str(abs(self)) + '+-'[self<0]

class ChargeList(list):
    """Just a list of :py:class:`Charge`s. When printed, looks like an
    enumeration of the list contents. Can also be constructed from such
    strings (e.g. "2+, 3+ and 4+").
    """
    def __init__(self, *args, **kwargs):
        if args and isinstance(args[0], str):
            self.extend(map(Charge,
                re.split(r'(?:,\s*)|(?:\s*and\s*)', args[0])))
        else:
            super(ChargeList, self).__init__(*args, **kwargs)

    def __str__(self):
        if len(self) > 1:
            return ', '.join(map(str, self[:-1])) + ' and {}'.format(self[-1])
        elif self:
            return str(self[0])
        return super(ChargeList, self).__str__()

def print_tree(d, indent_str=' -> ', indent_count=1):
    """Read a nested dict (with strings as keys) and print its structure.
    """
    def structure(d):
        out = {}
        for k, v in d.items():
            if isinstance(v, dict):
                out[k] = structure(v)
            elif isinstance(v, list) and v and isinstance(v[0], dict):
                out['{} [list]'.format(k)] = structure(v[0])
            else:
                out[k] = None
        return out

    def _print(d, level=0):
        for k, v in d.items():
            print('{}{}'.format(indent_str * indent_count * level, k))
            if v is not None:
                _print(v, level+1)
    _print(structure(d))

def memoize(maxsize=1000):
    """Make a memoization decorator. A negative value of `maxsize` means
    no size limit."""
    def deco(f):
        """Memoization decorator. Items of `kwargs` must be hashable."""
        memo = {}
        @wraps(f)
        def func(*args, **kwargs):
            key = (args, frozenset(kwargs.items()))
            if key not in memo:
                if len(memo) == maxsize:
                    memo.popitem()
                memo[key] = f(*args, **kwargs)
            return memo[key]
        return func
    return deco


### Public API ends here ###

### Next section: File reading helpers
def _keepstate(func):
    """Decorator to help keep the position in open file passed as first argument
    to functions"""
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

class _file_obj(object):
    """Check if `f` is a file name and open the file in `mode`.
    A context manager."""
    def __init__(self, f, mode):
        if f is None:
            self.file = {'r': sys.stdin, 'a': sys.stdout, 'w': sys.stdout
                    }[mode[0]]
            self.none = True
        elif isinstance(f, str):
            self.file = open(f, mode)
        else:
            self.file = f
        self.close_file = (self.file is not f)
    def __enter__(self):
        return self
    def __exit__(self, *args, **kwargs):
        if (not self.close_file) or hasattr(self, 'none'):
            return  # do nothing
        # clean up
        exit = getattr(self.file, '__exit__', None)
        if exit is not None:
            return exit(*args, **kwargs)
        else:
            exit = getattr(self.file, 'close', None)
            if exit is not None:
                exit()
    def __getattr__(self, attr):
        return getattr(self.file, attr)
    def __iter__(self):
        return iter(self.file)

def _file_reader(mode='r'):
    # a lot of the code below is borrowed from
    # http://stackoverflow.com/a/14095585/1258041
    def decorator(func):
        """A decorator implementing the context manager protocol for functions
        that read files.

        Note: 'close' must be in kwargs! Otherwise it won't be respected."""
        class CManager(object):
            def __init__(self, source, *args, **kwargs):
                self.file = _file_obj(source, mode)
                try:
                    self.reader = func(self.file, *args, **kwargs)
                except:  # clean up on any error
                    self.__exit__(*sys.exc_info())
                    raise

            # context manager support
            def __enter__(self):
                return self

            def __exit__(self, *args, **kwargs):
                self.file.__exit__(*args, **kwargs)

            # iterator support
            def __iter__(self):
                return self

            def __next__(self):
                try:
                    return next(self.reader)
                except StopIteration:
                    self.__exit__(None, None, None)
                    raise

            next = __next__  # Python 2 support

            # delegate everything else to file object
            def __getattr__(self, attr):
                return getattr(self.file, attr)

        @wraps(func)
        def helper(*args, **kwargs):
            return CManager(*args, **kwargs)
        return helper
    return decorator

def _make_chain(reader):
    def _iter(files, kwargs):
        for f in files:
            with reader(f, **kwargs) as r:
                for item in r:
                    yield item

    def chain(*files, **kwargs):
        return _iter(files, kwargs)

    def from_iterable(files, **kwargs):
        return _iter(files, kwargs)

    @contextmanager
    def _chain(*files, **kwargs):
        """Chain :py:func:`read` for several files.
        Positional arguments should be file names or file objects.
        Keyword arguments are passed to the :py:func:`read` function.
        """
        yield chain(*files, **kwargs)

    @contextmanager
    def _from_iterable(files, **kwargs):
        """Chain :py:func:`read` for several files.
        Keyword arguments are passed to the :py:func:`read` function.

        Parameters
        ----------
        files : iterable
            Iterable of file names or file objects.
        """
        yield from_iterable(files, **kwargs)


    _chain.from_iterable = _from_iterable
    return _chain

### End of file helpers section ###

def _parse_charge(s, list_only=False):
    if not list_only:
        try:
            return Charge(s)
        except PyteomicsError:
            pass
    return ChargeList(s)


### XML-related stuff below ###

def _local_name(element):
    """Strip namespace from the XML element's name"""
    if element.tag.startswith('{'):
        return element.tag.rsplit('}', 1)[1]
    return element.tag

def _make_version_info(env):
    @_keepstate
    def version_info(source):
        with _file_obj(source, 'rb') as s:
            s.seek(0)
            for _, elem in etree.iterparse(s, events=('start',),
                    remove_comments=True):
                if _local_name(elem) == env['element']:
                    vinfo = elem.attrib.get('version'), elem.attrib.get((
                        '{{{}}}'.format(elem.nsmap['xsi'])
                        if 'xsi' in elem.nsmap else '') + 'schemaLocation')
                    break
            return vinfo
    version_info.__doc__ = """
        Provide version information about the {0} file.

        Parameters:
        -----------
        source : str or file
            {0} file object or path to file

        Returns:
        --------
        out : tuple
            A (version, schema URL) tuple, both elements are strings or None.
        """.format(env['format'])
    return version_info

def _make_schema_info(env):
    @memoize(100)
    def _schema_info(source, **kwargs):
        read_schema = kwargs.get('read_schema', True)
        if not read_schema: return env['defaults']
        version, schema = env['version_info'](source)
        if version == env['default_version']:
            return env['defaults']
        ret = {}
        try:
            if not schema:
                schema_url = ''
                raise PyteomicsError(
                        'Schema information not found in {}.'.format(
                            getattr(source, 'name', source)))
            schema_url = schema.split()[-1]
            if not (schema_url.startswith('http://') or
                    schema_url.startswith('file://')):
                schema_url = 'file://' + schema_url
            schema_file = urlopen(schema_url)
            p = etree.XMLParser(remove_comments=True)
            schema_tree = etree.parse(schema_file, parser=p)
            types = {'ints': {'int', 'long', 'nonNegativeInteger',
                        'positiveInt', 'integer', 'unsignedInt'},
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
                                    schema_tree.iter() if x.attrib.get(
                                        'type', '').split(':')[-1] ==
                                    anc.attrib['name']]
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
                warnings.warn("Can't get the {0[format]} schema for version "
                "`{1}` from <{2}> at the moment.\n"
                "Using defaults for {0[default_version]}.\n"
                "You can disable reading the schema by specifying "
                "`read_schema=False`.".format(env, version,
                    schema_url))
            else:
                warnings.warn("Unknown {0[format]} version `{1}`. "
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
                    "".format(env, version, schema_url,
                        format_exc()))
            ret = env['defaults']
        return ret
    _schema_info.__doc__ = """
        Stores defaults for version {}, tries to retrieve the schema for
        other versions. Keys are: 'floats', 'ints', 'bools', 'lists',
        'intlists', 'floatlists', 'charlists'.""".format(env['default_version'])
    return _schema_info

def _make_get_info(env):
    def _get_info(source, element, recursive=False, retrieve_refs=False, **kw):
        """Extract info from element's attributes, possibly recursive.
        <cvParam> and <userParam> elements are treated in a special way."""
        name = _local_name(element)
        kwargs = dict(recursive=recursive, retrieve_refs=retrieve_refs)
        kwargs.update(kw)
        schema_info = env['schema_info'](source, **kwargs)
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
                cname = _local_name(child)
                if cname in {'cvParam', 'userParam'}:
                    newinfo = _get_info(source, child, **kw)
                    if not ('name' in info and 'name' in newinfo):
                        info.update(newinfo)
                    else:
                        if not isinstance(info['name'], list):
                            info['name'] = [info['name']]
                        info['name'].append(newinfo.pop('name'))
                else:
                    if cname not in schema_info['lists']:
                        info[cname] = env['get_info_smart'](
                                source, child, **kwargs)
                    else:
                        info.setdefault(cname, []).append(
                                env['get_info_smart'](source, child, **kwargs))
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
            raise PyteomicsError('Cannot convert string to bool: ' + s)

        def str_to_num(s, numtype):
            return numtype(s) if s else None

        to = lambda t: lambda s: str_to_num(s, t)

        converters = {'ints': to(int), 'floats': to(float), 'bools': str_to_bool,
                'intlists': lambda x: np.fromstring(x.replace('\n', ' '),
                    dtype=int, sep=' '),
                'floatlists': lambda x: np.fromstring(x.replace('\n', ' '),
                    sep=' '),
                'charlists': list}
        for k, v in info.items():
            for t, a in converters.items():
                if (_local_name(element), k) in schema_info[t]:
                    info[k] = a(v)
        # resolve refs
        # loop is needed to resolve refs pulled from other refs
        if retrieve_refs:
            while True:
                refs = False
                for k, v in dict(info).items():
                    if k.endswith('_ref'):
                        refs = True
                        info.update(env['get_by_id'](source, v))
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

def _make_iterfind(env):
    pattern_path = re.compile('([\w/*]*)(\[(\w+[<>=]{1,2}\w+)\])?')
    pattern_cond = re.compile('^\s*(\w+)\s*([<>=]{,2})\s*(\w+)$')
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
                return getattr(operator, func[sign])(d[lhs], type(d[lhs])(rhs))
            return False
        except (AttributeError, KeyError, ValueError):
            raise PyteomicsError('Invalid condition: ' + cond)

    @_keepstate
    def iterfind(source, path, **kwargs):
        """Parse ``source`` and yield info on elements with specified local name
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
        nodes = path.rstrip('/').lower().split('/')
        localname = nodes[0]
        found = False
        for ev, elem in etree.iterparse(source, events=('start', 'end'),
                remove_comments=True):
            name_lc = _local_name(elem).lower()
            if ev == 'start':
                if name_lc == localname or localname == '*':
                    found += True
            else:
                if name_lc == localname or localname == '*':
                    if (absolute and elem.getparent() is None
                            ) or not absolute:
                        for child in get_rel_path(elem, nodes[1:]):
                            info = env['get_info_smart'](source, child, **kwargs)
                            if cond is None or satisfied(info, cond):
                                yield info
                    if not localname == '*':
                        found -= 1
                if not found:
                    elem.clear()
    return iterfind

def _xpath(tree, path, ns=None):
    """Return the results of XPath query with added namespaces.
    Assumes the ns declaration is on the root element or absent."""
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
