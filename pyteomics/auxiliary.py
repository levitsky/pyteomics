"""
auxiliary - common functions and objects 
========================================

Math
----

  :py:func:`linear_regression` - a wrapper for numpy linear regression

Project infrastructure
----------------------

  :py:class:`PyteomicsError` - a pyteomics-specific exception

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

import numpy
from functools import wraps
from lxml import etree
from warnings import warn
from traceback import format_exc
try: # Python 2.7
    from urllib import urlopen
except ImportError: # Python 3.x
    from urllib.request import urlopen

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

    if (a!=None and b==None):
        b = numpy.mean([y[i] - a * x[i] for i in range(len(x))])
    elif (a!=None and b!= None):
        pass
    else:
        a, b = numpy.polyfit(x, y, 1)

    r = numpy.corrcoef(x, y)[0, 1]
    stderr = numpy.std([y[i] - a * x[i] - b for i in range(len(x))])

    return (a, b, r, stderr)

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

def _local_name(element):
    """Strip namespace from the XML element's name"""
    if element.tag.startswith('{'):
        return element.tag.rsplit('}', 1)[1]
    else:
        return element.tag

def _make_version_info(env):
    @_keepstate
    def version_info(source):
        for _, elem in etree.iterparse(source, events=('start',),
                remove_comments=True):
            if _local_name(elem) == env['element']:
                return elem.attrib.get('version'), elem.attrib.get((
                    '{{{}}}'.format(elem.nsmap['xsi'])
                    if 'xsi' in elem.nsmap else '') + 'schemaLocation')
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
    _schema_info_cache = {}
    def _schema_info(source, key):
        if source in _schema_info_cache:
            return _schema_info_cache[source][key]
        
        version, schema = env['version_info'](source)
        if version == env['default_version']:
            ret = env['defaults']
        else:
            ret = {}
            try:
                if not schema:
                    schema_url = ''
                    raise PyteomicsError(
                            'Schema information not found in {}.'.format(source))
                schema_url = schema.split()[-1]
                schema_file = urlopen(schema_url)
                p = etree.XMLParser(remove_comments=True)
                schema_tree = etree.parse(schema_file, parser=p)
                types = {'ints': {'int', 'long', 'nonNegativeInteger'},
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
                                    elname = schema_tree.find(
                                            '//*[@type="{}"]'.format(
                                        anc.attrib['name'])).attrib['name']
                                else:
                                    elname = anc.attrib['name']
                                tuples.add(
                                        (elname, elem.attrib['name']))
                    ret[k] = tuples
                ret['lists'] = set(elem.attrib['name'] for elem in schema_tree.xpath(
                    '//*[local-name()="element"]') if 'name' in elem.attrib and
                    elem.attrib.get( 'maxOccurs', '1') != '1')
            except Exception as e:
                warn("Unknown {} version `{}`. Attempt to use schema\n"
                        "information from <{}> failed.\n{}\n"
                        "Falling back to defaults for {}".format(
                            env['format'], version, schema_url,
                            format_exc(), env['default_version']))
                ret = env['defaults']
        _schema_info_cache[source] = ret
        return ret[key]
    _schema_info.__doc__ = """
        Stores defaults for version {}, tries to retrieve the schema for
        other versions. Keys are: 'floats', 'ints', 'bools', 'lists',
        'intlists', 'floatlists', 'charlists'.""".format(env['default_version'])
    return _schema_info

