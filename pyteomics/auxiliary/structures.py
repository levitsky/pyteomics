import re
from collections import defaultdict, Counter

try:
    basestring
except NameError:
    basestring = (str, bytes)


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
                except Exception:
                    pass
            raise PyteomicsError(*e.args)

    def __str__(self):
        return str(abs(self)) + '+-'[self < 0]


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


class unitint(int):
    def __new__(cls, value, unit_info=None):
        inst = super(unitint, cls).__new__(cls, value)
        inst.unit_info = unit_info
        return inst

    def _repr_pretty_(self, p, cycle):
        base = super(unitint, self).__repr__()
        if self.unit_info:
            string = "%s %s" % (base, self.unit_info)
        else:
            string = base
        p.text(string)


class unitfloat(float):
    def __new__(cls, value, unit_info=None):
        inst = super(unitfloat, cls).__new__(cls, value)
        inst.unit_info = unit_info
        return inst

    def _repr_pretty_(self, p, cycle):
        base = super(unitfloat, self).__repr__()
        if self.unit_info:
            string = "%s %s" % (base, self.unit_info)
        else:
            string = base
        p.text(string)


class unitstr(str):
    def __new__(cls, value, unit_info=None):
        inst = super(unitstr, cls).__new__(cls, value)
        inst.unit_info = unit_info
        return inst

    def _repr_pretty_(self, p, cycle):
        base = super(unitstr, self).__repr__()
        if self.unit_info:
            string = "%s %s" % (base, self.unit_info)
        else:
            string = base
        p.text(string)


class cvstr(str):
    '''A helper class to associate a controlled vocabullary accession
    number with an otherwise plain :class:`str` object'''
    def __new__(cls, value, accession=None, unit_accession=None):
        inst = super(cvstr, cls).__new__(cls, value)
        inst.accession = accession
        inst.unit_accession = unit_accession
        return inst


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

    def _walk_dict(self, data, index):
        for key, value in data.items():
            accession = self._accession(key)
            if accession:
                if value != '':
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


cvquery = CVQueryEngine()
