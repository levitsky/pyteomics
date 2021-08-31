import re
from collections import defaultdict, Counter
import warnings

try:
    basestring
    PY2 = True
except NameError:
    basestring = (str, bytes)
    PY2 = False


_UNIT_CV_INTERN_TABLE = dict()


def clear_unit_cv_table():
    """Clear the module-level unit name and
    controlled vocabulary accession table.
    """
    _UNIT_CV_INTERN_TABLE.clear()


def _intern_unit_or_cv(unit_or_cv):
    """Intern `unit_or_cv` in :const:`~._UNIT_CV_INTERN_TABLE`, potentially
    keeping a reference to the object stored for the duration of the program.

    Parameters
    ----------
    unit_or_cv : object
        The value to intern

    Returns
    -------
    object:
        The object which `unit_or_cv` hash-equals in :const:`~._UNIT_CV_INTERN_TABLE`.
    """
    if unit_or_cv is None:
        return None
    try:
        return _UNIT_CV_INTERN_TABLE[unit_or_cv]
    except KeyError:
        _UNIT_CV_INTERN_TABLE[unit_or_cv] = unit_or_cv
        return _UNIT_CV_INTERN_TABLE[unit_or_cv]


class PyteomicsError(Exception):
    """Exception raised for errors in Pyteomics library.

    Attributes
    ----------
    message : str
        Error message.
    """

    def __init__(self, msg, *values):
        self.message = msg
        self.values = values

    def __str__(self):
        if not self.values:
            return "Pyteomics error, message: %s" % (repr(self.message),)
        else:
            return "Pyteomics error, message: %s %r" % (repr(self.message), self.values)


class Charge(int):
    """A subclass of :py:class:`int`. Can be constructed from strings in "N+"
    or "N-" format, and the string representation of a :py:class:`Charge` is
    also in that format.
    """
    def __new__(cls, *args, **kwargs):
        try:
            return super(Charge, cls).__new__(cls, *args)
        except ValueError as e:
            if isinstance(args[0], basestring):
                try:
                    num, sign = re.match(r'^(\d+)(\+|-)$', args[0]).groups()
                    return super(Charge, cls).__new__(cls, sign + num, *args[1:], **kwargs)
                except Exception:
                    pass
            raise PyteomicsError(*e.args)

    def __str__(self):
        return str(abs(self)) + '+-'[self < 0]


class Ion(str):
    """Represents an Ion, right now just a subclass of String.
    """
    _pattern = r'([abcxyz]\d+(\-H2O|\-NH3)?)([\+|-]\d+)'  # "y2-H2O+1"

    def __init__(self, *args, **kwargs):
        if args and isinstance(args[0], basestring):
            try:
                self.ion_type, self.neutral_loss, self.charge = re.match(self._pattern, args[0]).groups()
            except Exception:
                raise PyteomicsError("Malformed ion string, must match the regex {!r}".format(self._pattern))


class ChargeList(list):
    """Just a list of :py:class:`Charge`s. When printed, looks like an
    enumeration of the list contents. Can also be constructed from such
    strings (e.g. "2+, 3+ and 4+").
    """

    def __init__(self, *args, **kwargs):
        if args and isinstance(args[0], basestring):
            delim = r'(?:,\s*)|(?:\s*and\s*)'
            self.extend(map(Charge, re.split(delim, args[0])))
        else:
            try:
                super(ChargeList, self).__init__(
                    sorted(set(args[0])), *args[1:], **kwargs)
            except Exception:
                super(ChargeList, self).__init__(*args, **kwargs)
            self[:] = map(Charge, self)

    def __str__(self):
        if len(self) > 1:
            return ', '.join(map(str, self[:-1])) + ' and {}'.format(self[-1])
        elif self:
            return str(self[0])
        return super(ChargeList, self).__str__()


def _parse_charge(s, list_only=False):
    if not list_only:
        try:
            return Charge(s)
        except PyteomicsError:
            pass
    return ChargeList(s)


def _parse_ion(ion_text):
    try:
        return Ion(ion_text)
    except Exception as e:
        warnings.warn('Could not parse ion string: {} ({})'.format(ion_text, e.args[0]))


class BasicComposition(defaultdict, Counter):
    """A generic dictionary for compositions.
    Keys should be strings, values should be integers.
    Allows simple arithmetics."""

    def __init__(self, *args, **kwargs):
        defaultdict.__init__(self, int)
        Counter.__init__(self, *args, **kwargs)
        for k, v in list(self.items()):
            if not v:
                del self[k]

    def __str__(self):
        return '{}({})'.format(type(self).__name__, dict.__repr__(self))

    def __repr__(self):
        return str(self)

    def _repr_pretty_(self, p, cycle):
        if cycle:  # should never happen
            p.text('{} object with a cyclic reference'.format(type(self).__name__))
        p.text(str(self))

    def __add__(self, other):
        result = self.copy()
        for elem, cnt in other.items():
            result[elem] += cnt
        return result

    def __iadd__(self, other):
        for elem, cnt in other.items():
            self[elem] += cnt
        return self

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        result = self.copy()
        for elem, cnt in other.items():
            result[elem] -= cnt
        return result

    def __isub__(self, other):
        for elem, cnt in other.items():
            self[elem] -= cnt
        return self

    def __rsub__(self, other):
        return (self - other) * (-1)

    def __mul__(self, other):
        if not isinstance(other, int):
            raise PyteomicsError('Cannot multiply Composition by non-integer',
                                 other)
        return type(self)({k: v * other for k, v in self.items()})

    def __imul__(self, other):
        if not isinstance(other, int):
            raise PyteomicsError('Cannot multiply Composition by non-integer',
                                 other)
        for elem in self:
            self[elem] *= other
        return self

    def __rmul__(self, other):
        return self * other

    def __eq__(self, other):
        if not isinstance(other, dict):
            return False
        self_items = {i for i in self.items() if i[1]}
        other_items = {i for i in other.items() if i[1]}
        return self_items == other_items

    # override default behavior:
    # we don't want to add 0's to the dictionary
    def __missing__(self, key):
        return 0

    def __setitem__(self, key, value):
        if isinstance(value, float):
            value = int(round(value))
        elif not isinstance(value, int):
            raise PyteomicsError('Only integers allowed as values in '
                                 'Composition, got {}.'.format(type(value).__name__))
        if value:  # reject 0's
            super(BasicComposition, self).__setitem__(key, value)
        elif key in self:
            del self[key]

    def copy(self):
        return type(self)(self)

    def __reduce__(self):
        class_, args, state, list_iterator, dict_iterator = super(
            BasicComposition, self).__reduce__()
        # Override the reduce of defaultdict so we do not provide the
        # `int` type as the first argument
        # which prevents from correctly unpickling the object
        args = ()
        return class_, args, state, list_iterator, dict_iterator


class _MappingOverAttributeProxy(object):
    '''A replacement for __dict__ for unpickling an object which once
    has __slots__ now but did not before.'''

    def __init__(self, obj):
        self.obj = obj

    def __getitem__(self, key):
        return getattr(self.obj, key)

    def __setitem__(self, key, value):
        setattr(self.obj, key, value)

    def __contains__(self, key):
        return hasattr(self.obj, key)

    def __repr__(self):
        return "{self.__class__.__name__}({self.obj})".format(self=self)


class unitint(int):
    '''Represents an integer value with a unit name.

    Behaves identically to a built-in :class:`int` type.

    Attributes
    ----------
    unit_info : :class:`str`
        The name of the unit this value posseses.
    '''
    def __new__(cls, value, unit_info=None):
        inst = int.__new__(cls, value)
        inst.unit_info = unit_info
        return inst

    def __reduce__(self):
        return self.__class__, (int(self), self.unit_info)

    def _repr_pretty_(self, p, cycle):
        base = super(unitint, self).__repr__()
        if self.unit_info:
            string = "%s %s" % (base, self.unit_info)
        else:
            string = base
        p.text(string)


class unitfloat(float):
    '''Represents an float value with a unit name.

    Behaves identically to a built-in :class:`float` type.

    Attributes
    ----------
    unit_info : :class:`str`
        The name of the unit this value posseses.
    '''
    __slots__ = ('unit_info', )

    def __new__(cls, value, unit_info=None):
        inst = float.__new__(cls, value)
        inst.unit_info = unit_info
        return inst

    @property
    def __dict__(self):
        return _MappingOverAttributeProxy(self)

    def __reduce__(self):
        return self.__class__, (float(self), self.unit_info)

    def _repr_pretty_(self, p, cycle):
        base = super(unitfloat, self).__repr__()
        if self.unit_info:
            string = "%s %s" % (base, self.unit_info)
        else:
            string = base
        p.text(string)


class unitstr(str):
    '''Represents an string value with a unit name.

    Behaves identically to a built-in :class:`str` type.

    Attributes
    ----------
    unit_info : :class:`str`
        The name of the unit this value posseses.
    '''
    if not PY2:
        __slots__ = ("unit_info", )

    def __new__(cls, value, unit_info=None):
        if PY2 and isinstance(value, unicode):
            value = value.encode('utf-8')
        inst = str.__new__(cls, value)
        inst.unit_info = unit_info
        return inst

    @property
    def __dict__(self):
        return _MappingOverAttributeProxy(self)

    def __reduce__(self):
        return self.__class__, (str(self), self.unit_info)

    def _repr_pretty_(self, p, cycle):
        base = super(unitstr, self).__repr__()
        if self.unit_info:
            string = "%s %s" % (base, self.unit_info)
        else:
            string = base
        p.text(string)


class cvstr(str):
    '''A helper class to associate a controlled vocabullary accession
    number with an otherwise plain :class:`str` object

    Attributes
    ----------
    accession : str
        The accession number for this parameter, e.g. MS:1000040
    unit_accession : str
        The accession number for the unit of the value, if any
    '''

    if not PY2:
        __slots__ = ('accession', 'unit_accession')

    _cache = {}

    def __new__(cls, value, accession=None, unit_accession=None):
        try:
            inst = cls._cache[value]
            if inst.accession == accession and inst.unit_accession == unit_accession:
                return inst
        except KeyError:
            pass

        if PY2 and isinstance(value, unicode):
            value = value.encode('utf-8')
        inst = str.__new__(cls, value)
        inst.accession = _intern_unit_or_cv(accession)
        inst.unit_accession = _intern_unit_or_cv(unit_accession)
        cls._cache[value] = inst
        return inst

    @property
    def __dict__(self):
        return _MappingOverAttributeProxy(self)

    def __reduce__(self):
        return self.__class__, (str(self), self.accession, self.unit_accession)


class CVQueryEngine(object):
    '''Traverse an arbitrarily nested dictionary looking
    for keys which are :class:`cvstr` instances, or objects
    with an attribute called ``accession``.
    '''

    def _accession(self, key):
        return getattr(key, 'accession', None)

    def _query_dict(self, data, accession):
        for key, value in data.items():
            if self._accession(key) == accession:
                if not isinstance(value, str) or value != '':
                    return value
                else:
                    return key
            elif isinstance(value, dict):
                inner = self._query_dict(value, accession)
                if inner is not None:
                    return inner
            elif isinstance(value, (list, tuple)):
                inner = self._query_sequence(value, accession)
                if inner is not None:
                    return inner
            elif self._accession(value) == accession:
                return value

    def _query_sequence(self, data, accession):
        for value in data:
            if isinstance(value, dict):
                inner = self._query_dict(value, accession)
                if inner is not None:
                    return inner
            elif isinstance(value, (list, tuple)):
                inner = self._query_sequence(value, accession)
                if inner is not None:
                    return inner
            elif self._accession(value) == accession:
                return value

    def query(self, data, accession):
        '''Search ``data`` for a key with the accession
        number ``accession``. Returns :const:`None` if
        not found.
        '''
        if accession is None:
            raise TypeError("`accession` cannot be None")
        return self._query_dict(data, accession)

    def _is_empty(self, value):
        if isinstance(value, basestring):
            return value == ''
        return False

    def _walk_dict(self, data, index):
        for key, value in data.items():
            accession = self._accession(key)
            if accession:
                if not self._is_empty(value):
                    index[accession] = value
                else:
                    index[accession] = key
            elif isinstance(value, dict):
                self._walk_dict(value, index)
            elif isinstance(value, (list, tuple)):
                self._walk_sequence(value, index)
            accession = self._accession(value)
            if accession:
                index[accession] = value
        return index

    def _walk_sequence(self, data, index):
        for value in data:
            if isinstance(value, dict):
                self._walk_dict(value, index)
            elif isinstance(value, (list, tuple)):
                self._walk_sequence(value, index)
            else:
                accession = self._accession(value)
                if accession:
                    index[accession] = value

    def index(self, data):
        '''Construct a flat :class:`dict` whose keys are the
        accession numbers for all qualified keys in ``data``
        and whose values are the mapped values from ``data``.
        '''
        index = self._walk_dict(data, {})
        return index

    def __call__(self, data, accession=None):
        '''If ``accession`` is :const:`None`, calls
        :meth:`index` on ``data``, otherwise calls
        :meth:`query` with ``data`` and ``accession``.
        '''
        if accession is None:
            return self.index(data)
        else:
            return self.query(data, accession)

'''A ready-to-use instance of :class:`~.CVQueryEngine`'''
cvquery = CVQueryEngine()
