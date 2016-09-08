"""
auxiliary - common functions and objects
========================================

Math
----

  :py:func:`linear_regression_vertical` - a wrapper for NumPy linear regression,
  minimizes the sum of squares of *y* errors.

  :py:func:`linear_regression` - alias for :py:func:`linear_regression_vertical`.

  :py:func:`linear_regression_perpendicular` - a wrapper for NumPy linear regression,
  minimizes the sum of squares of (perpendicular) distances between the points and the line.


Target-Decoy Approach
---------------------

  :py:func:`qvalues` - estimate q-values for a set of PSMs.

  :py:func:`filter` - filter PSMs to specified FDR level using TDA or given PEPs.

  :py:func:`filter.chain` - a chained version of :py:func:`filter`.

  :py:func:`fdr` - estimate FDR in a set of PSMs using TDA or given PEPs.

Project infrastructure
----------------------

  :py:class:`PyteomicsError` - a pyteomics-specific exception.

Helpers
-------

  :py:class:`Charge` - a subclass of :py:class:`int` for charge states.

  :py:class:`ChargeList` - a subclass of :py:class:`list` for lists of charges.

  :py:func:`print_tree` - display the structure of a complex nested
  :py:class:`dict`.

  :py:func:`memoize` - makes a
  `memoization <http://stackoverflow.com/a/1988826/1258041>`_
  `function decorator <http://stackoverflow.com/a/1594484/1258041>`_.

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
from functools import wraps
from traceback import format_exc
import re
import zlib
import base64
import operator as op
import sys
from contextlib import contextmanager
import types
from bisect import bisect_right
from collections import Counter, defaultdict
try:
    from collections import Container, Sized
except ImportError:
    from collections.abc import Container, Sized
import math
from distutils.version import LooseVersion
try:
    basestring = basestring
except NameError:
    basestring = (str, bytes)
try:
    import pandas as pd
except ImportError:
    pd = None
else:
    if hasattr(pd, '_version'):
        pv = pd._version.get_versions()['version']
    else:
        pv = pd.version.version
    if LooseVersion(pv) < LooseVersion('0.17'):
        pd.DataFrame.sort_values = pd.DataFrame.sort

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

def linear_regression_vertical(x, y=None, a=None, b=None):
    """Calculate coefficients of a linear regression y = a * x + b.
    The fit minimizes *vertical* distances between the points and the line.

    Requires :py:mod:`numpy`.

    Parameters
    ----------
    x, y : array_like of float
        1-D arrays of floats. If `y` is omitted, `x` must be a 2-D array of shape (N, 2).
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

    import numpy as np
    x = np.array(x, copy=False)
    if y is not None:
        y = np.array(y, copy=False)
    else:
        if len(x.shape) != 2 or x.shape[-1] != 2:
            raise PyteomicsError(
                'If `y` is not given, x.shape should be (N, 2), given: {}'.format(x.shape))
        y = x[:, 1]
        x = x[:, 0]
    if (a is not None and b is None):
        b = (y - a * x).mean()
    elif (a is not None and b is not None):
        pass
    else:
        a, b = np.polyfit(x, y, 1)

    r = np.corrcoef(x, y)[0, 1]
    stderr = (y - a * x - b).std()

    return a, b, r, stderr

def linear_regression(x, y=None, a=None, b=None):
    """Alias of :py:func:`linear_regression_vertical`."""
    return linear_regression_vertical(x, y, a, b)

def linear_regression_perpendicular(x, y=None):
    """Calculate coefficients of a linear regression y = a * x + b.
    The fit minimizes *perpendicular* distances between the points and the line.

    Requires :py:mod:`numpy`.

    Parameters
    ----------
    x, y : array_like of float
        1-D arrays of floats. If `y` is omitted, `x` must be a 2-D array of shape (N, 2).

    Returns
    -------
    out : 4-tuple of float
        The structure is (a, b, r, stderr), where
        a -- slope coefficient,
        b -- free term,
        r -- Peason correlation coefficient,
        stderr -- standard deviation.
    """

    import numpy as np
    x = np.array(x, copy=False)
    if y is not None:
        y = np.array(y, copy=False)
        data = np.hstack((x.reshape((-1, 1)), y.reshape((-1, 1))))
    else:
        if len(x.shape) != 2 or x.shape[-1] != 2:
            raise PyteomicsError(
                'If `y` is not given, x.shape should be (N, 2), given: {}'.format(x.shape))
        data = x
    mu = data.mean(axis=0)
    eigenvectors, eigenvalues, V = np.linalg.svd((data-mu).T, full_matrices=False)
    a = eigenvectors[0][1] / eigenvectors[0][0]
    xm, ym = data.mean(axis=0)
    b = ym - a * xm

    r = np.corrcoef(data[:, 0], data[:, 1])[0, 1]
    stderr = ((data[:, 1] - a * data[:, 0] - b) / np.sqrt(a**2 + 1)).std()

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
            try:
                super(ChargeList, self).__init__(
                        sorted(set(args[0])), *args[1:], **kwargs)
            except:
                super(ChargeList, self).__init__(*args, **kwargs)
            self[:] = map(Charge, self)

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
        if cycle: # should never happen
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
        return type(self)({k: v*other for k, v in self.items()})

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
        if isinstance(value, float): value = int(round(value))
        elif not isinstance(value, int):
            raise PyteomicsError('Only integers allowed as values in '
                         'Composition, got {}.'.format(type(value).__name__))
        if value: # reject 0's
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


### Public API ends here ###

### Next section: File reading helpers
def _keepstate(func):
    """Decorator to help keep the position in open files passed as
    positional arguments to functions"""
    @wraps(func)
    def wrapped(*args, **kwargs):
        positions = [getattr(arg, 'seek', None) and
                getattr(arg, 'tell', type(None))() for arg in args]
        for arg, pos in zip(args, positions):
            if pos is not None:
                arg.seek(0)
        res = func(*args, **kwargs)
        for arg, pos in zip(args, positions):
            if pos is not None:
                try:
                    arg.seek(pos)
                except ValueError:
                    pass
        return res
    return wrapped

class _file_obj(object):
    """Check if `f` is a file name and open the file in `mode`.
    A context manager."""
    def __init__(self, f, mode):
        if f is None:
            self.file = {'r': sys.stdin, 'a': sys.stdout, 'w': sys.stdout
                    }[mode[0]]
            self.none = True
        elif isinstance(f, basestring):
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

class IteratorContextManager(object):
    def __init__(self, _func, *args, **kwargs):
        self._func = _func
        self._args = args
        self._kwargs = kwargs
        if type(self) == IteratorContextManager:
            self.reset()

    def reset(self):
        """Resets the iterator to its initial state."""
        try:
            self._reader = self._func(*self._args, **self._kwargs)
        except:
            self.__exit__(*sys.exc_info())
            raise

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        pass

    def __iter__(self):
        return self

    def __next__(self):
        try:
            return next(self._reader)
        except StopIteration:
            self.__exit__(None, None, None)
            raise

    next = __next__

class FileReader(IteratorContextManager):
    """Abstract class implementing context manager protocol
    for file readers.
    """
    def __init__(self, source, mode, func, pass_file, args, kwargs):
        super(FileReader, self).__init__(func, *args, **kwargs)
        self._pass_file = pass_file
        self._source_init = source
        self._mode = mode
        self.reset()

    def reset(self):
        if hasattr(self, '_source'):
            self._source.__exit__(None, None, None)
        self._source = _file_obj(self._source_init, self._mode)
        try:
            if self._pass_file:
                self._reader = self._func(
                        self._source, *self._args, **self._kwargs)
            else:
                self._reader = self._func(*self._args, **self._kwargs)
        except:  # clean up on any error
            self.__exit__(*sys.exc_info())
            raise

    def __exit__(self, *args, **kwargs):
        self._source.__exit__(*args, **kwargs)

    # delegate everything else to file object
    def __getattr__(self, attr):
        if attr == '_source':
            raise AttributeError
        return getattr(self._source, attr)

def _file_reader(_mode='r'):
    # a lot of the code below is borrowed from
    # http://stackoverflow.com/a/14095585/1258041
    def decorator(_func):
        """A decorator implementing the context manager protocol for functions
        that read files.

        Note: 'close' must be in kwargs! Otherwise it won't be respected.
        """
        @wraps(_func)
        def helper(*args, **kwargs):
            return FileReader(args[0], _mode, _func, True, args[1:], kwargs)
        return helper
    return decorator

def _file_writer(_mode='a'):
    def decorator(_func):
        """A decorator that opens output files for writer functions.
        """
        @wraps(_func)
        def helper(*args, **kwargs):
            m = kwargs.pop('file_mode', _mode)
            if len(args) > 1:
                with _file_obj(args[1], m) as out:
                    return _func(args[0], out, *args[2:], **kwargs)
            else:
                with _file_obj(kwargs.pop('output', None), m) as out:
                    return _func(*args, output=out, **kwargs)
        return helper
    return decorator

def _make_chain(reader, readername, full_output=False):

    def concat_results(*args, **kwargs):
        results = [reader(arg, **kwargs) for arg in args]
        if pd is not None and all(isinstance(a, pd.DataFrame) for a in args):
            return pd.concat(results)
        else:
            return np.concatenate(results)

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
        yield chain(*files, **kwargs)

    @contextmanager
    def _from_iterable(files, **kwargs):
        yield from_iterable(files, **kwargs)

    def dispatch(*args, **kwargs):
        return dispatch_from_iterable(args, **kwargs)

    def dispatch_from_iterable(args, **kwargs):
        if kwargs.get('full_output', full_output):
            return concat_results(*args, **kwargs)
        else:
            return _chain(*args, **kwargs)

    dispatch.__doc__ = """Chain :py:func:`{0}` for several files.
        Positional arguments should be file names or file objects.
        Keyword arguments are passed to the :py:func:`{0}` function.
        """.format(readername)
    dispatch_from_iterable.__doc__ = """Chain :py:func:`{0}` for several files.
        Keyword arguments are passed to the :py:func:`{0}` function.

        Parameters
        ----------
        files : iterable
            Iterable of file names or file objects.
        """.format(readername)
    dispatch.from_iterable = dispatch_from_iterable

    return dispatch

### End of file helpers section ###

### Target-decoy approach stuff ###

def _fix_docstring(f, **defaults):
    for argname, v in defaults.items():
        if v is not None:
            f.__doc__ = re.sub('{} : .*'.format(argname), lambda m: m.group() + ', optional', f.__doc__)

def _make_qvalues(read, is_decoy, key):
    """Create a function that reads PSMs from a file and calculates q-values
    for each value of `key`."""

    def qvalues(*args, **kwargs):
        """Read `args` and return a NumPy array with scores and q-values.
        q-values are calculated either using TDA or based on provided values of PEP.

        Requires :py:mod:`numpy` (and optionally :py:mod:`pandas`).

        Parameters
        ----------

        positional args : file or str
            Files to read PSMs from. All positional arguments are treated as
            files. The rest of the arguments must be named.

        key : callable / array-like / iterable / str, keyword only
            If callable, a function used for sorting of PSMs. Should accept
            exactly one argument (PSM) and return a number (the smaller the better).
            If array-like, should contain scores for all given PSMs.
            If string, it is used as a field name (PSMs must be in a record array
            or a :py:class:`DataFrame`).

            .. warning::
                The default function may not work
                with your files, because format flavours are diverse.

        reverse : bool, keyword only, optional
            If :py:const:`True`, then PSMs are sorted in descending order,
            i.e. the value of the key function is higher for better PSMs.
            Default is :py:const:`False`.

        is_decoy : callable / array-like / iterable / str, keyword only
            If callable, a function used to determine if the PSM is decoy or not.
            Should accept exactly one argument (PSM) and return a truthy value if the
            PSM should be considered decoy.
            If array-like, should contain boolean values for all given PSMs.
            If string, it is used as a field name (PSMs must be in a record array
            or a :py:class:`DataFrame`).

            .. warning::
                The default function may not work
                with your files, because format flavours are diverse.

        pep : callable / array-like / iterable / str, keyword only, optional
            If callable, a function used to determine the posterior error probability (PEP).
            Should accept exactly one argument (PSM) and return a float.
            If array-like, should contain float values for all given PSMs.
            If string, it is used as a field name (PSMs must be in a record array
            or a :py:class:`DataFrame`).

            .. note:: If this parameter is given, then PEP values will be used to calculate
               q-values. Otherwise, decoy PSMs will be used instead. This option conflicts with:
               `is_decoy`, `remove_decoy`, `formula`, `ratio`, `correction`.
               `key` can still be provided. Without `key`, PSMs will be sorted by PEP.

        remove_decoy : bool, keyword only, optional
            Defines whether decoy matches should be removed from the output.
            Default is :py:const:`False`.

            .. note:: If set to :py:const:`False`, then by default the decoy
               PSMs will be taken into account when estimating FDR. Refer to the
               documentation of :py:func:`fdr` for math; basically, if
               `remove_decoy` is :py:const:`True`, then formula 1 is used
               to control output FDR, otherwise it's formula 2. This can be
               changed by overriding the `formula` argument.

        formula : int, keyword only, optional
            Can be either 1 or 2, defines which formula should be used for FDR
            estimation. Default is 1 if `remove_decoy` is :py:const:`True`,
            else 2 (see :py:func:`fdr` for definitions).

        ratio : float, keyword only, optional
            The size ratio between the decoy and target databases. Default is
            1. In theory, the "size" of the database is the number of
            theoretical peptides eligible for assignment to spectra that are
            produced by *in silico* cleavage of that database.

        correction : int or float, keyword only, optional
            Possible values are 0, 1 and 2, or floating point numbers between 0 and 1.
            Default is 0 (no correction); 1 accounts for the probability that a false
            positive scores better than the first excluded decoy PSM; 2 also corrects
            that probability for finite size of the sample. If a floating point number
            is given, then instead of the expectation value for the number of false PSMs,
            the confidence value is used. The value of `correction` is then interpreted as
            desired confidence level. E.g., if correction=0.95, then the calculated q-values
            do not exceed the "real" q-values with 95% probability.

        full_output : bool, keyword only, optional
            If :py:const:`True`, then the returned array has PSM objects along
            with scores and q-values. Default is :py:const:`False`.

        **kwargs : passed to the :py:func:`chain` function.

        Returns
        -------
        out : numpy.ndarray
            A sorted array of records with the following fields:

            - 'score': :py:class:`np.float64`
            - 'is decoy': :py:class:`np.bool_`
            - 'q': :py:class:`np.float64`
            - 'psm': :py:class:`np.object_` (if `full_output` is :py:const:`True`)
        """
        import numpy as np
        @_keepstate
        def get_scores(*args, **kwargs):
            scores = []
            with read(*args, **kwargs) as f:
                for i, psm in enumerate(f):
                    row = []
                    for func in (keyf, isdecoy):
                        if callable(func):
                            row.append(func(psm))
                        elif isinstance(func, basestring):
                            row.append(psm[func])
                        else:
                            row.append(func[i])
                    row.append(None)
                    if full:
                        row.append(psm)
                    scores.append(tuple(row))
            return scores

        peps = kwargs.pop('pep', None)
        if peps is not None:
            x = {'is_decoy', 'remove_decoy', 'formula', 'ratio', 'correction'}.intersection(kwargs)
            if x:
                raise PyteomicsError("Can't use these parameters with `pep`: " + ', '.join(x))
        keyf = kwargs.pop('key', key)
        reverse = kwargs.pop('reverse', False)
        if keyf is None:
            keyf = peps
            if reverse:
                raise PyteomicsError('reverse = True when using PEPs for sorting')
        if not callable(keyf) and not isinstance(keyf, (Sized, Container)):
            keyf = np.array(list(keyf))
        ratio = kwargs.pop('ratio', 1)
        remove_decoy = kwargs.pop('remove_decoy', False)
        if peps is None:
            isdecoy = kwargs.pop('is_decoy', is_decoy)
        else:
            isdecoy = peps
        if not callable(isdecoy) and not isinstance(isdecoy, (Sized, Container)):
            isdecoy = np.array(list(isdecoy))
        full = kwargs.pop('full_output', False)
        formula = kwargs.pop('formula', (2, 1)[bool(remove_decoy)])
        correction = kwargs.pop('correction', 0)
        if formula not in {1, 2}:
            raise PyteomicsError('`formula` must be either 1 or 2')

        decoy_or_pep_label = 'is decoy' if peps is None else 'PEP'

        fields = [('score', np.float64),
            (decoy_or_pep_label, np.bool_ if peps is None else np.float64),
            ('q', np.float64)]
        # if all args are NumPy arrays with common dtype, use it in the output
        if full:
            dtypes = {getattr(arg, 'dtype', None) for arg in args}
            if len(dtypes) == 1 and None not in dtypes:
                psm_dtype = dtypes.pop()
            else:
                psm_dtype = np.object_
            dtype = np.dtype(fields + [('psm', psm_dtype)])
        else:
            dtype = np.dtype(fields)

        arr_flag = False
        psms = None
        if callable(keyf) or callable(isdecoy):
            scores = np.array(get_scores(*args, **kwargs), dtype=dtype)
        else:
            if pd is not None and all(isinstance(arg, pd.DataFrame) for arg in args):
                psms = pd.concat(args)
            elif all(isinstance(arg, np.ndarray) for arg in args):
                psms = np.concatenate(args)

            if not isinstance(keyf, basestring):
                keyf = np.array(keyf)
                arr_flag = True
            if not isinstance(isdecoy, basestring):
                isdecoy = np.array(isdecoy)
                arr_flag = True

            if arr_flag:
                scores = np.empty(keyf.size if hasattr(keyf, 'size') else isdecoy.size,
                    dtype=dtype)
                for func, label in zip((keyf, isdecoy), ('score', decoy_or_pep_label)):
                    if not isinstance(func, basestring):
                        scores[label] = func
                    else:
                        scores[label] = psms[func]
            else:
                scores = np.empty(psms.shape[0], dtype=fields)
                scores['score'] = psms[keyf]
                scores[decoy_or_pep_label] = psms[isdecoy]

        if not scores.size:
            if full and psms is not None:
                return psms
            return scores

        if not reverse:
            keys = scores[decoy_or_pep_label], scores['score']
        else:
            keys = scores[decoy_or_pep_label], -scores['score']
        lexsort = np.lexsort(keys)
        scores = scores[lexsort]
        if psms is not None:
            if pd is not None and isinstance(psms, pd.DataFrame):
                if arr_flag:
                    psms = psms.iloc[lexsort]
                else:
                    psms.sort_values([keyf, isdecoy], ascending=[not reverse, True], inplace=True)
            else:
                psms = psms[lexsort]

        cumsum = scores[decoy_or_pep_label].cumsum(dtype=np.float64)
        tfalse = cumsum.copy()
        ind = np.arange(1., scores.shape[0] + 1., dtype=np.float64)

        if peps is not None:
            q = cumsum / ind
        else:
            if isinstance(correction, int):
                if correction == 1:
                    tfalse += 1
                elif correction == 2:
                    p = 1. / (1. + ratio)
                    targ = ind - cumsum
                    for i in range(tfalse.size):
                        tfalse[i] = _expectation(cumsum[i], targ[i], p)
            elif 0 < correction < 1:
                p = 1. / (1. + ratio)
                targ = ind - cumsum
                for i in range(tfalse.size):
                    tfalse[i] = _confidence_value(correction, cumsum[i], targ[i], p)
            elif correction:
                raise PyteomicsError('Invalid value for `correction`.')

            if formula == 1:
                q = tfalse / (ind - cumsum) / ratio
            else:
                q = (cumsum + tfalse / ratio) / ind
        # Make sure that q-values are equal for equal scores (conservatively)
        # and that q-values are monotonic
        for i in range(scores.size-1, 0, -1):
            if (scores['score'][i] == scores['score'][i-1] or
                    q[i-1] > q[i]):
                q[i-1] = q[i]
        scores['q'] = q
        if remove_decoy:
            if psms is not None:
                psms = psms[~scores[decoy_or_pep_label]]
            scores = scores[~scores[decoy_or_pep_label]]

        if full and psms is not None:
            if isinstance(psms, np.ndarray):
                fields = sorted(psms.dtype.fields, key=lambda x: psms.dtype.fields[x][1])
                extra = []
                for func, label in zip((keyf, isdecoy), ('score', decoy_or_pep_label)):
                    if not (isinstance(func, basestring) or label in psms.dtype.fields):
                        extra.append(label)
                    elif label in psms.dtype.fields:
                        psms[label] = scores[label]
                newdt = [(name, psms.dtype.fields[name][0]) for name in fields] + [
                    (name, np.float64) for name in extra] + [('q', np.float64)]
                psms_ = psms
                psms = np.empty_like(psms_, dtype=newdt)
                for f in fields:
                    psms[f] = psms_[f]
                for f in extra:
                    psms[f] = scores[f]
            else:
                for func, label in zip((keyf, isdecoy), ('score', decoy_or_pep_label)):
                    if not isinstance(label, basestring):
                        psms[label] = scores[label]
            psms['q'] = scores['q']
            return psms
        return scores
    
    _fix_docstring(qvalues, is_decoy=is_decoy, key=key)
    if read is _iter:
        qvalues.__doc__ = qvalues.__doc__.replace("""positional args : file or str
            Files to read PSMs from. All positional arguments are treated as
            files.""", """positional args : iterables
            Iterables to read PSMs from. All positional arguments are chained."""
            ).replace(""".. warning::
                The default function may not work
                with your files, because format flavours are diverse.\n""", "")

    return qvalues


def _make_filter(read, is_decoy, key, qvalues):
    """Create a function that reads PSMs from a file and filters them to
    the desired FDR level (estimated by TDA), returning the top PSMs
    sorted by `key`.
    """
    def filter(*args, **kwargs):
        import numpy as np
        try:
            fdr = kwargs.pop('fdr')
        except KeyError:
            raise PyteomicsError('Keyword argument required: fdr')

        args = [list(arg) if not isinstance(arg, (Container, Sized))
                else arg for arg in args]
        peps = kwargs.get('pep')
        if peps is None:
            remove_decoy = kwargs.pop('remove_decoy', True)
            scores = qvalues(*args, remove_decoy=remove_decoy, **kwargs)
        else:
            scores = qvalues(*args, **kwargs)
        keyf = kwargs.pop('key', key)
        if keyf is None: keyf = peps
        reverse = kwargs.pop('reverse', False)
        better = [op.lt, op.gt][bool(reverse)]
        isdecoy = kwargs.pop('is_decoy', is_decoy)
        kwargs.pop('formula', None)
        try:
            i = scores['q'].searchsorted(fdr, side='right')
            if isinstance(i, Sized):
                i = i[0]
        except AttributeError:
            i = bisect_right(scores['q'], fdr)
        if kwargs.pop('full_output', False):
            if pd is not None and isinstance(scores, pd.DataFrame):
                return scores.iloc[:i]
            elif callable(keyf) or callable(isdecoy):
                return scores['psm'][:i]
            else:
                return scores[:i]
        elif not scores.size:
            return (_ for _ in ())
        cutoff = scores['score'][i] if i < scores.size else (
                scores['score'][-1] + (1, -1)[bool(reverse)])
        def out():
            with read(*args, **kwargs) as f:
                for p, s in zip(f, scores):
                    if peps is not None or not remove_decoy or not s['is decoy']:
                        if better(s['score'], cutoff):
                            yield p
        return out()

    def _filter(*args, **kwargs):
        """Read `args` and yield only the PSMs that form a set with
        estimated false discovery rate (FDR) not exceeding `fdr`.

        Requires :py:mod:`numpy` and, optionally, :py:mod:`pandas`.

        Parameters
        ----------
        positional args : file or str
            Files to read PSMs from. All positional arguments are treated as
            files. The rest of the arguments must be named.

        fdr : float, keyword only, 0 <= fdr <= 1
            Desired FDR level.

        key : callable, keyword only
            A function used for sorting of PSMs. Should accept exactly one
            argument (PSM) and return a number (the smaller the better). The
            default is a function that tries to extract e-value from the PSM.

            .. warning::
                The default function may not work
                with your files, because format flavours are diverse.

        reverse : bool, keyword only, optional
            If :py:const:`True`, then PSMs are sorted in descending order,
            i.e. the value of the key function is higher for better PSMs.
            Default is :py:const:`False`.

        is_decoy : callable, keyword only
            A function used to determine if the PSM is decoy or not. Should
            accept exactly one argument (PSM) and return a truthy value if the
            PSM should be considered decoy.

            .. warning::
                The default function may not work
                with your files, because format flavours are diverse.

        remove_decoy : bool, keyword only, optional
            Defines whether decoy matches should be removed from the output.
            Default is :py:const:`True`.

            .. note:: If set to :py:const:`False`, then by default the decoy
               PSMs will be taken into account when estimating FDR. Refer to the
               documentation of :py:func:`fdr` for math; basically, if
               `remove_decoy` is :py:const:`True`, then formula 1 is used
               to control output FDR, otherwise it's formula 2. This can be
               changed by overriding the `formula` argument.

        formula : int, keyword only, optional
            Can be either 1 or 2, defines which formula should be used for FDR
            estimation. Default is 1 if `remove_decoy` is :py:const:`True`,
            else 2 (see :py:func:`fdr` for definitions).

        ratio : float, keyword only, optional
            The size ratio between the decoy and target databases. Default is
            1. In theory, the "size" of the database is the number of
            theoretical peptides eligible for assignment to spectra that are
            produced by *in silico* cleavage of that database.

        correction : int or float, keyword only, optional
            Possible values are 0, 1 and 2, or floating point numbers between 0 and 1.
            Default is 0 (no correction); 1 accounts for the probability that a false
            positive scores better than the first excluded decoy PSM; 2 also corrects
            that probability for finite size of the sample. If a floating point number
            is given, then instead of the expectation value for the number of false PSMs,
            the confidence value is used. The value of `correction` is then interpreted as
            desired confidence level. E.g., if correction=0.95, then the calculated q-values
            do not exceed the "real" q-values with 95% probability.

        pep : callable / array-like / iterable / str, keyword only, optional
            If callable, a function used to determine the posterior error probability (PEP).
            Should accept exactly one argument (PSM) and return a float.
            If array-like, should contain float values for all given PSMs.
            If string, it is used as a field name (PSMs must be in a record array
            or a :py:class:`DataFrame`).

            .. note:: If this parameter is given, then PEP values will be used to calculate
               q-values. Otherwise, decoy PSMs will be used instead. This option conflicts with:
               `is_decoy`, `remove_decoy`, `formula`, `ratio`, `correction`.
               `key` can still be provided. Without `key`, PSMs will be sorted by PEP.

        full_output : bool, keyword only, optional
            If :py:const:`True`, then an array of PSM objects is returned.
            Otherwise, an iterator / context manager object is returned, and the
            files are parsed twice. This saves some RAM, but is ~2x slower.
            Default is :py:const:`True`.

            .. note:: The name for the parameter comes from the fact that it is
                      internally passed to :py:func:`qvalues`.

        **kwargs : passed to the :py:func:`chain` function.

        Returns
        -------
        out : iterator or :py:class:`numpy.ndarray` or :py:class:`pandas.DataFrame`
        """
        if kwargs.pop('full_output', True):
            return filter(*args, full_output=True, **kwargs)
        return IteratorContextManager(filter, *args, **kwargs)

    _fix_docstring(_filter, is_decoy=is_decoy, key=key)
    if read is _iter:
        _filter.__doc__ = _filter.__doc__.replace("""positional args : file or str
            Files to read PSMs from. All positional arguments are treated as
            files.""", """positional args : iterables
            Iterables to read PSMs from. All positional arguments are chained."""
            ).replace(""".. warning::
                The default function may not work
                with your files, because format flavours are diverse.\n""", "")
    return _filter

@contextmanager
def _itercontext(x, **kw):
    try:
        yield (row for i, row in x.iterrows())
    except AttributeError:
        yield x
_iter = _make_chain(_itercontext, 'iter')
qvalues = _make_qvalues(_iter, None, None)

filter = _make_filter(_iter, None, None, qvalues)
filter.chain = _make_chain(filter, 'filter', True)

try:
    import numpy as np
    _precalc_fact = np.log([math.factorial(n) for n in range(20)])
    def log_factorial(x):
        x = np.array(x)
        pf = _precalc_fact
        m = (x >= pf.size)
        out = np.empty(x.shape)
        out[~m] = pf[x[~m].astype(int)]
        x = x[m]
        out[m] = x * np.log(x) - x + 0.5 * np.log(2 * np.pi * x)
        return out

    def _expectation(d, T, p=0.5):
        if T is None: return d+1
        T = np.array(T, dtype=int)
        m = np.arange(T.max()+1, dtype=int)
        pi = np.exp(_log_pi(d, m, p))
        return ((m * pi).cumsum() / pi.cumsum())[T]

    def _confidence_value(conf, d, T, p=0.5):
        if T is not None:
            T = np.array(T, dtype=int)
            m = np.arange(T.max()+1, dtype=int)
        else:
            m = np.arange(max(50*d, 10000))
        log_pi = _log_pi(d, m, p)
        pics = np.exp(log_pi).cumsum()
        return np.searchsorted(pics, conf * (pics[T] if T is not None else 1))

except ImportError:
    def log_factorial(n):
        if n > 10:
            return n * math.log(n) - n + 0.5 * math.log(2*math.pi*n)
        else:
            return math.log(math.factorial(n))

def _log_pi_r(d, k, p=0.5):
    return k * math.log(p) + log_factorial(k + d) - log_factorial(k) - log_factorial(d)

def _log_pi(d, k, p=0.5):
    return _log_pi_r(d, k, p) + (d + 1) * math.log(1 - p)

def _make_fdr(is_decoy):
    def fdr(psms=None, formula=1, is_decoy=is_decoy, ratio=1, correction=0, pep=None):
        """Estimate FDR of a data set using TDA or given PEP values.
        Two formulas can be used. The first one (default) is:

        .. math::

                FDR = \\frac{N_{decoy}}{N_{target} * ratio}

        The second formula is:

        .. math::

                FDR = \\frac{N_{decoy} * (1 + \\frac{1}{ratio})}{N_{total}}

        .. note::
            This function is less versatile than :py:func:`qvalues`. To obtain FDR,
            you can call :py:func:`qvalues` and take the last q-value. This function
            can be used (with `correction = 0` or `1`) when :py:mod:`numpy` is not available.

        Parameters
        ----------
        psms : iterable, optional
            An iterable of PSMs, e.g. as returned by :py:func:`read`.
            Not needed if `is_decoy` is an iterable.

        formula : int, optional
            Can be either 1 or 2, defines which formula should be used for FDR
            estimation. Default is 1.

        is_decoy : callable, iterable, or str
            If callable, should accept exactly one argument (PSM) and return a truthy value
            if the PSM is considered decoy. Default is :py:func:`is_decoy`.
            If array-like, should contain float values for all given PSMs.
            If string, it is used as a field name (PSMs must be in a record array
            or a :py:class:`pandas.DataFrame`).

            .. warning::
                The default function may not work
                with your files, because format flavours are diverse.

        pep : callable, iterable, or str, optional
            If callable, a function used to determine the posterior error probability (PEP).
            Should accept exactly one argument (PSM) and return a float.
            If array-like, should contain float values for all given PSMs.
            If string, it is used as a field name (PSMs must be in a record array
            or a :py:class:`pandas.DataFrame`).

            .. note:: If this parameter is given, then PEP values will be used to calculate FDR.
               Otherwise, decoy PSMs will be used instead. This option conflicts with:
               `is_decoy`, `formula`, `ratio`, `correction`.

        ratio : float, optional
            The size ratio between the decoy and target databases. Default is 1.
            In theory, the "size" of the database is the number of
            theoretical peptides eligible for assignment to spectra that are
            produced by *in silico* cleavage of that database.

        correction : int or float, optional
            Possible values are 0, 1 and 2, or floating point numbers between 0 and 1.
            Default is 0 (no correction); 1 accounts for the probability that a false
            positive scores better than the first excluded decoy PSM; 2 also corrects
            that probability for finite size of the sample. If a floating point number
            is given, then instead of the expectation value for the number of false PSMs,
            the confidence value is used. The value of `correction` is then interpreted as
            desired confidence level. E.g., if correction=0.95, then the calculated q-values
            do not exceed the "real" q-values with 95% probability.

            .. note::
                Requires :py:mod:`numpy`, if `correction` is a float or 2.

            .. note::
                Correction is only needed if the PSM set at hand was obtained using TDA-based
                filtering based on decoy counting (as done by using :py:func:`!filter` without
                `correction`).

        Returns
        -------
        out : float
            The estimation of FDR, (roughly) between 0 and 1.
        """
        if formula not in {1, 2}:
            raise PyteomicsError('`formula` must be either 1 or 2.')
        total, decoy = 0, 0
        if pep is not None:
            is_decoy = pep
        if isinstance(is_decoy, basestring):
            decoy = psms[is_decoy].sum()
            total = psms.shape[0]
        elif callable(is_decoy):
            for psm in psms:
                total += 1
                d = is_decoy(psm)
                decoy += d if pep is not None else bool(d)
        else:
            if not isinstance(is_decoy, (Sized, Container)):
                is_decoy = list(is_decoy)
            if pep is not None:
                decoy = sum(is_decoy)
            else:
                decoy = sum(map(bool, is_decoy))
            total = len(is_decoy)
        if pep is not None:
            return float(decoy) / total
        tfalse = decoy
        if correction == 1 or (correction == 2 and total/decoy > 10):
            tfalse += 1
        elif correction == 2:
            p = 1. / (1. + ratio)
            tfalse = _expectation(decoy, total-decoy, p)
        elif 0 < correction < 1:
            p = 1. / (1. + ratio)
            tfalse = _confidence_value(correction, decoy, total-decoy, p)
        if formula == 1:
            return float(tfalse) / (total - decoy) / ratio
        return (decoy + tfalse / ratio) / total
    
    _fix_docstring(fdr, is_decoy=is_decoy)
    if is_decoy is None:
        fdr.__doc__ = fdr.__doc__.replace(""".. warning::
                The default function may not work
                with your files, because format flavours are diverse.\n""", "")
    return fdr

fdr = _make_fdr(None)

### Miscellaneous ###

def _parse_charge(s, list_only=False):
    if not list_only:
        try:
            return Charge(s)
        except PyteomicsError:
            pass
    return ChargeList(s)

def _decode_base64_data_array(source, dtype, is_compressed):
    """Read a base64-encoded binary array.

    Parameters
    ----------
    source : str
        A binary array encoded with base64.
    dtype : dtype
        The type of the array in numpy dtype notation.
    is_compressed : bool
        If True then the array will be decompressed with zlib.

    Returns
    -------
    out : numpy.array
    """

    decoded_source = base64.b64decode(source.encode('ascii'))
    if is_compressed:
        decoded_source = zlib.decompress(decoded_source)
    output = np.frombuffer(decoded_source, dtype=dtype)
    return output

### Bulky constants for other modules are defined below.

_nist_mass = {'Ac': {0: (227, 1.0),
  206: (206.0145, 0.0),
  207: (207.01195, 0.0),
  208: (208.01155, 0.0),
  209: (209.00949, 0.0),
  210: (210.00944, 0.0),
  211: (211.00773, 0.0),
  212: (212.00781, 0.0),
  213: (213.00661, 0.0),
  214: (214.006902, 0.0),
  215: (215.006454, 0.0),
  216: (216.00872, 0.0),
  217: (217.009347, 0.0),
  218: (218.01164, 0.0),
  219: (219.01242, 0.0),
  220: (220.014763, 0.0),
  221: (221.01559, 0.0),
  222: (222.017844, 0.0),
  223: (223.019137, 0.0),
  224: (224.021723, 0.0),
  225: (225.02323, 0.0),
  226: (226.026098, 0.0),
  227: (227.0277521, 0.0),
  228: (228.0310211, 0.0),
  229: (229.03302, 0.0),
  230: (230.03629, 0.0),
  231: (231.03856, 0.0),
  232: (232.04203, 0.0),
  233: (233.04455, 0.0),
  234: (234.04842, 0.0),
  235: (235.05123, 0.0),
  236: (236.0553, 0.0)},
 'Ag': {0: (106.905097, 1.0),
  93: (92.94978, 0.0),
  94: (93.94278, 0.0),
  95: (94.93548, 0.0),
  96: (95.93068, 0.0),
  97: (96.92397, 0.0),
  98: (97.92157, 0.0),
  99: (98.9176, 0.0),
  100: (99.9161, 0.0),
  101: (100.9128, 0.0),
  102: (101.91169, 0.0),
  103: (102.908973, 0.0),
  104: (103.908629, 0.0),
  105: (104.906529, 0.0),
  106: (105.906669, 0.0),
  107: (106.905097, 0.51839),
  108: (107.905956, 0.0),
  109: (108.904752, 0.48161),
  110: (109.906107, 0.0),
  111: (110.905291, 0.0),
  112: (111.907005, 0.0),
  113: (112.906567, 0.0),
  114: (113.908804, 0.0),
  115: (114.90876, 0.0),
  116: (115.91136, 0.0),
  117: (116.91168, 0.0),
  118: (117.91458, 0.0),
  119: (118.91567, 0.0),
  120: (119.91879, 0.0),
  121: (120.91985, 0.0),
  122: (121.92353, 0.0),
  123: (122.9249, 0.0),
  124: (123.92864, 0.0),
  125: (124.93043, 0.0),
  126: (125.9345, 0.0),
  127: (126.93677, 0.0),
  128: (127.94117, 0.0),
  129: (128.94369, 0.0),
  130: (129.95045, 0.0)},
 'Al': {0: (26.98153863, 1.0),
  21: (21.02804, 0.0),
  22: (22.01952, 0.0),
  23: (23.007267, 0.0),
  24: (23.9999389, 0.0),
  25: (24.9904281, 0.0),
  26: (25.98689169, 0.0),
  27: (26.98153863, 1.0),
  28: (27.98191031, 0.0),
  29: (28.980445, 0.0),
  30: (29.98296, 0.0),
  31: (30.983947, 0.0),
  32: (31.98812, 0.0),
  33: (32.99084, 0.0),
  34: (33.99685, 0.0),
  35: (34.99986, 0.0),
  36: (36.00621, 0.0),
  37: (37.01068, 0.0),
  38: (38.01723, 0.0),
  39: (39.02297, 0.0),
  40: (40.03145, 0.0),
  41: (41.03833, 0.0),
  42: (42.04689, 0.0)},
 'Am': {0: (243, 1.0),
  231: (231.04556, 0.0),
  232: (232.04659, 0.0),
  233: (233.04635, 0.0),
  234: (234.04781, 0.0),
  235: (235.04795, 0.0),
  236: (236.04958, 0.0),
  237: (237.05, 0.0),
  238: (238.05198, 0.0),
  239: (239.0530245, 0.0),
  240: (240.0553, 0.0),
  241: (241.0568291, 0.0),
  242: (242.0595492, 0.0),
  243: (243.0613811, 0.0),
  244: (244.0642848, 0.0),
  245: (245.066452, 0.0),
  246: (246.069775, 0.0),
  247: (247.07209, 0.0),
  248: (248.07575, 0.0),
  249: (249.07848, 0.0)},
 'Ar': {0: (39.9623831225, 1.0),
  30: (30.02156, 0.0),
  31: (31.01212, 0.0),
  32: (31.997638, 0.0),
  33: (32.9899257, 0.0),
  34: (33.9802712, 0.0),
  35: (34.9752576, 0.0),
  36: (35.967545106, 0.003365),
  37: (36.96677632, 0.0),
  38: (37.9627324, 0.000632),
  39: (38.964313, 0.0),
  40: (39.9623831225, 0.996003),
  41: (40.9645006, 0.0),
  42: (41.963046, 0.0),
  43: (42.965636, 0.0),
  44: (43.964924, 0.0),
  45: (44.96804, 0.0),
  46: (45.96809, 0.0),
  47: (46.97219, 0.0),
  48: (47.97454, 0.0),
  49: (48.98052, 0.0),
  50: (49.98443, 0.0),
  51: (50.99163, 0.0),
  52: (51.99678, 0.0),
  53: (53.00494, 0.0)},
 'As': {0: (74.9215965, 1.0),
  60: (59.99313, 0.0),
  61: (60.98062, 0.0),
  62: (61.9732, 0.0),
  63: (62.96369, 0.0),
  64: (63.95757, 0.0),
  65: (64.94956, 0.0),
  66: (65.94471, 0.0),
  67: (66.93919, 0.0),
  68: (67.93677, 0.0),
  69: (68.93227, 0.0),
  70: (69.93092, 0.0),
  71: (70.927112, 0.0),
  72: (71.926752, 0.0),
  73: (72.923825, 0.0),
  74: (73.9239287, 0.0),
  75: (74.9215965, 1.0),
  76: (75.922394, 0.0),
  77: (76.9206473, 0.0),
  78: (77.921827, 0.0),
  79: (78.920948, 0.0),
  80: (79.922534, 0.0),
  81: (80.922132, 0.0),
  82: (81.9245, 0.0),
  83: (82.92498, 0.0),
  84: (83.92906, 0.0),
  85: (84.93202, 0.0),
  86: (85.9365, 0.0),
  87: (86.9399, 0.0),
  88: (87.94494, 0.0),
  89: (88.94939, 0.0),
  90: (89.9555, 0.0),
  91: (90.96043, 0.0),
  92: (91.9668, 0.0)},
 'At': {0: (210, 1.0),
  193: (192.99984, 0.0),
  194: (193.99873, 0.0),
  195: (194.996268, 0.0),
  196: (195.99579, 0.0),
  197: (196.99319, 0.0),
  198: (197.99284, 0.0),
  199: (198.99053, 0.0),
  200: (199.990351, 0.0),
  201: (200.988417, 0.0),
  202: (201.98863, 0.0),
  203: (202.986942, 0.0),
  204: (203.987251, 0.0),
  205: (204.986074, 0.0),
  206: (205.986667, 0.0),
  207: (206.985784, 0.0),
  208: (207.98659, 0.0),
  209: (208.986173, 0.0),
  210: (209.987148, 0.0),
  211: (210.9874963, 0.0),
  212: (211.990745, 0.0),
  213: (212.992937, 0.0),
  214: (213.996372, 0.0),
  215: (214.998653, 0.0),
  216: (216.002423, 0.0),
  217: (217.004719, 0.0),
  218: (218.008694, 0.0),
  219: (219.011162, 0.0),
  220: (220.01541, 0.0),
  221: (221.01805, 0.0),
  222: (222.02233, 0.0),
  223: (223.02519, 0.0)},
 'Au': {0: (196.9665687, 1.0),
  169: (168.99808, 0.0),
  170: (169.99612, 0.0),
  171: (170.991879, 0.0),
  172: (171.99004, 0.0),
  173: (172.986237, 0.0),
  174: (173.98476, 0.0),
  175: (174.98127, 0.0),
  176: (175.9801, 0.0),
  177: (176.976865, 0.0),
  178: (177.97603, 0.0),
  179: (178.973213, 0.0),
  180: (179.972521, 0.0),
  181: (180.970079, 0.0),
  182: (181.969618, 0.0),
  183: (182.967593, 0.0),
  184: (183.967452, 0.0),
  185: (184.965789, 0.0),
  186: (185.965953, 0.0),
  187: (186.964568, 0.0),
  188: (187.965324, 0.0),
  189: (188.963948, 0.0),
  190: (189.9647, 0.0),
  191: (190.9637, 0.0),
  192: (191.964813, 0.0),
  193: (192.96415, 0.0),
  194: (193.965365, 0.0),
  195: (194.9650346, 0.0),
  196: (195.96657, 0.0),
  197: (196.9665687, 1.0),
  198: (197.9682423, 0.0),
  199: (198.9687652, 0.0),
  200: (199.97073, 0.0),
  201: (200.971657, 0.0),
  202: (201.97381, 0.0),
  203: (202.975155, 0.0),
  204: (203.97772, 0.0),
  205: (204.97987, 0.0)},
 'B': {0: (11.0093054, 1.0),
  6: (6.04681, 0.0),
  7: (7.02992, 0.0),
  8: (8.0246072, 0.0),
  9: (9.0133288, 0.0),
  10: (10.012937, 0.199),
  11: (11.0093054, 0.801),
  12: (12.0143521, 0.0),
  13: (13.0177802, 0.0),
  14: (14.025404, 0.0),
  15: (15.031103, 0.0),
  16: (16.03981, 0.0),
  17: (17.04699, 0.0),
  18: (18.05617, 0.0),
  19: (19.06373, 0.0)},
 'Ba': {0: (137.9052472, 1.0),
  114: (113.95068, 0.0),
  115: (114.94737, 0.0),
  116: (115.94138, 0.0),
  117: (116.9385, 0.0),
  118: (117.93304, 0.0),
  119: (118.93066, 0.0),
  120: (119.92604, 0.0),
  121: (120.92405, 0.0),
  122: (121.9199, 0.0),
  123: (122.918781, 0.0),
  124: (123.915094, 0.0),
  125: (124.914473, 0.0),
  126: (125.91125, 0.0),
  127: (126.911094, 0.0),
  128: (127.908318, 0.0),
  129: (128.908679, 0.0),
  130: (129.9063208, 0.00106),
  131: (130.906941, 0.0),
  132: (131.9050613, 0.00101),
  133: (132.9060075, 0.0),
  134: (133.9045084, 0.02417),
  135: (134.9056886, 0.06592),
  136: (135.9045759, 0.07854),
  137: (136.9058274, 0.11232),
  138: (137.9052472, 0.71698),
  139: (138.9088413, 0.0),
  140: (139.910605, 0.0),
  141: (140.914411, 0.0),
  142: (141.916453, 0.0),
  143: (142.920627, 0.0),
  144: (143.922953, 0.0),
  145: (144.92763, 0.0),
  146: (145.93022, 0.0),
  147: (146.93495, 0.0),
  148: (147.93772, 0.0),
  149: (148.94258, 0.0),
  150: (149.94568, 0.0),
  151: (150.95081, 0.0),
  152: (151.95427, 0.0),
  153: (152.95961, 0.0)},
 'Be': {0: (9.0121822, 1.0),
  5: (5.04079, 0.0),
  6: (6.019726, 0.0),
  7: (7.01692983, 0.0),
  8: (8.0053051, 0.0),
  9: (9.0121822, 1.0),
  10: (10.0135338, 0.0),
  11: (11.021658, 0.0),
  12: (12.026921, 0.0),
  13: (13.03569, 0.0),
  14: (14.04289, 0.0),
  15: (15.05346, 0.0),
  16: (16.06192, 0.0)},
 'Bh': {0: (272, 1.0),
  260: (260.12197, 0.0),
  261: (261.12166, 0.0),
  262: (262.12289, 0.0),
  263: (263.12304, 0.0),
  264: (264.1246, 0.0),
  265: (265.12515, 0.0),
  266: (266.12694, 0.0),
  267: (267.12765, 0.0),
  268: (268.12976, 0.0),
  269: (269.13069, 0.0),
  270: (270.13362, 0.0),
  271: (271.13518, 0.0),
  272: (272.13803, 0.0),
  273: (273.13962, 0.0),
  274: (274.14244, 0.0),
  275: (275.14425, 0.0)},
 'Bi': {0: (208.9803987, 1.0),
  184: (184.00112, 0.0),
  185: (184.99763, 0.0),
  186: (185.9966, 0.0),
  187: (186.993158, 0.0),
  188: (187.99227, 0.0),
  189: (188.9892, 0.0),
  190: (189.9883, 0.0),
  191: (190.985786, 0.0),
  192: (191.98546, 0.0),
  193: (192.98296, 0.0),
  194: (193.98283, 0.0),
  195: (194.980651, 0.0),
  196: (195.980667, 0.0),
  197: (196.978864, 0.0),
  198: (197.97921, 0.0),
  199: (198.977672, 0.0),
  200: (199.978132, 0.0),
  201: (200.977009, 0.0),
  202: (201.977742, 0.0),
  203: (202.976876, 0.0),
  204: (203.977813, 0.0),
  205: (204.977389, 0.0),
  206: (205.978499, 0.0),
  207: (206.9784707, 0.0),
  208: (207.9797422, 0.0),
  209: (208.9803987, 1.0),
  210: (209.9841204, 0.0),
  211: (210.987269, 0.0),
  212: (211.9912857, 0.0),
  213: (212.994385, 0.0),
  214: (213.998712, 0.0),
  215: (215.00177, 0.0),
  216: (216.006306, 0.0),
  217: (217.00947, 0.0),
  218: (218.01432, 0.0)},
 'Bk': {0: (247, 1.0),
  235: (235.05658, 0.0),
  236: (236.05733, 0.0),
  237: (237.057, 0.0),
  238: (238.05828, 0.0),
  239: (239.05828, 0.0),
  240: (240.05976, 0.0),
  241: (241.06023, 0.0),
  242: (242.06198, 0.0),
  243: (243.063008, 0.0),
  244: (244.065181, 0.0),
  245: (245.0663616, 0.0),
  246: (246.06867, 0.0),
  247: (247.070307, 0.0),
  248: (248.07309, 0.0),
  249: (249.0749867, 0.0),
  250: (250.078317, 0.0),
  251: (251.08076, 0.0),
  252: (252.08431, 0.0),
  253: (253.08688, 0.0),
  254: (254.0906, 0.0)},
 'Br': {0: (78.9183371, 1.0),
  67: (66.96479, 0.0),
  68: (67.95852, 0.0),
  69: (68.95011, 0.0),
  70: (69.94479, 0.0),
  71: (70.93874, 0.0),
  72: (71.93664, 0.0),
  73: (72.93169, 0.0),
  74: (73.929891, 0.0),
  75: (74.925776, 0.0),
  76: (75.924541, 0.0),
  77: (76.921379, 0.0),
  78: (77.921146, 0.0),
  79: (78.9183371, 0.5069),
  80: (79.9185293, 0.0),
  81: (80.9162906, 0.4931),
  82: (81.9168041, 0.0),
  83: (82.91518, 0.0),
  84: (83.916479, 0.0),
  85: (84.915608, 0.0),
  86: (85.918798, 0.0),
  87: (86.920711, 0.0),
  88: (87.92407, 0.0),
  89: (88.92639, 0.0),
  90: (89.93063, 0.0),
  91: (90.93397, 0.0),
  92: (91.93926, 0.0),
  93: (92.94305, 0.0),
  94: (93.94868, 0.0),
  95: (94.95287, 0.0),
  96: (95.95853, 0.0),
  97: (96.9628, 0.0)},
 'C': {0: (12.0, 1.0),
  8: (8.037675, 0.0),
  9: (9.0310367, 0.0),
  10: (10.0168532, 0.0),
  11: (11.0114336, 0.0),
  12: (12.0, 0.9893),
  13: (13.0033548378, 0.0107),
  14: (14.003241989, 0.0),
  15: (15.0105993, 0.0),
  16: (16.014701, 0.0),
  17: (17.022586, 0.0),
  18: (18.02676, 0.0),
  19: (19.03481, 0.0),
  20: (20.04032, 0.0),
  21: (21.04934, 0.0),
  22: (22.0572, 0.0)},
 'Ca': {0: (39.96259098, 1.0),
  34: (34.01412, 0.0),
  35: (35.00494, 0.0),
  36: (35.99309, 0.0),
  37: (36.98587, 0.0),
  38: (37.976318, 0.0),
  39: (38.9707197, 0.0),
  40: (39.96259098, 0.96941),
  41: (40.96227806, 0.0),
  42: (41.95861801, 0.00647),
  43: (42.9587666, 0.00135),
  44: (43.9554818, 0.02086),
  45: (44.9561866, 0.0),
  46: (45.9536926, 4e-05),
  47: (46.954546, 0.0),
  48: (47.952534, 0.00187),
  49: (48.955674, 0.0),
  50: (49.957519, 0.0),
  51: (50.9615, 0.0),
  52: (51.9651, 0.0),
  53: (52.97005, 0.0),
  54: (53.97435, 0.0),
  55: (54.98055, 0.0),
  56: (55.98557, 0.0),
  57: (56.99236, 0.0)},
 'Cd': {0: (113.9033585, 1.0),
  95: (94.94987, 0.0),
  96: (95.93977, 0.0),
  97: (96.93494, 0.0),
  98: (97.9274, 0.0),
  99: (98.92501, 0.0),
  100: (99.92029, 0.0),
  101: (100.91868, 0.0),
  102: (101.91446, 0.0),
  103: (102.913419, 0.0),
  104: (103.909849, 0.0),
  105: (104.909468, 0.0),
  106: (105.906459, 0.0125),
  107: (106.906618, 0.0),
  108: (107.904184, 0.0089),
  109: (108.904982, 0.0),
  110: (109.9030021, 0.1249),
  111: (110.9041781, 0.128),
  112: (111.9027578, 0.2413),
  113: (112.9044017, 0.1222),
  114: (113.9033585, 0.2873),
  115: (114.905431, 0.0),
  116: (115.904756, 0.0749),
  117: (116.907219, 0.0),
  118: (117.906915, 0.0),
  119: (118.90992, 0.0),
  120: (119.90985, 0.0),
  121: (120.91298, 0.0),
  122: (121.91333, 0.0),
  123: (122.917, 0.0),
  124: (123.91765, 0.0),
  125: (124.92125, 0.0),
  126: (125.92235, 0.0),
  127: (126.92644, 0.0),
  128: (127.92776, 0.0),
  129: (128.93215, 0.0),
  130: (129.9339, 0.0),
  131: (130.94067, 0.0),
  132: (131.94555, 0.0)},
 'Ce': {0: (139.9054387, 1.0),
  119: (118.95276, 0.0),
  120: (119.94664, 0.0),
  121: (120.94342, 0.0),
  122: (121.93791, 0.0),
  123: (122.9354, 0.0),
  124: (123.93041, 0.0),
  125: (124.92844, 0.0),
  126: (125.92397, 0.0),
  127: (126.92273, 0.0),
  128: (127.91891, 0.0),
  129: (128.9181, 0.0),
  130: (129.91474, 0.0),
  131: (130.91442, 0.0),
  132: (131.91146, 0.0),
  133: (132.911515, 0.0),
  134: (133.908925, 0.0),
  135: (134.909151, 0.0),
  136: (135.907172, 0.00185),
  137: (136.907806, 0.0),
  138: (137.905991, 0.00251),
  139: (138.906653, 0.0),
  140: (139.9054387, 0.8845),
  141: (140.9082763, 0.0),
  142: (141.909244, 0.11114),
  143: (142.912386, 0.0),
  144: (143.913647, 0.0),
  145: (144.91723, 0.0),
  146: (145.91876, 0.0),
  147: (146.92267, 0.0),
  148: (147.92443, 0.0),
  149: (148.9284, 0.0),
  150: (149.93041, 0.0),
  151: (150.93398, 0.0),
  152: (151.93654, 0.0),
  153: (152.94058, 0.0),
  154: (153.94342, 0.0),
  155: (154.94804, 0.0),
  156: (155.95126, 0.0),
  157: (156.95634, 0.0)},
 'Cf': {0: (251, 1.0),
  237: (237.06207, 0.0),
  238: (238.06141, 0.0),
  239: (239.06242, 0.0),
  240: (240.0623, 0.0),
  241: (241.06373, 0.0),
  242: (242.0637, 0.0),
  243: (243.06543, 0.0),
  244: (244.066001, 0.0),
  245: (245.068049, 0.0),
  246: (246.0688053, 0.0),
  247: (247.071001, 0.0),
  248: (248.072185, 0.0),
  249: (249.0748535, 0.0),
  250: (250.0764061, 0.0),
  251: (251.079587, 0.0),
  252: (252.081626, 0.0),
  253: (253.085133, 0.0),
  254: (254.087323, 0.0),
  255: (255.09105, 0.0),
  256: (256.09344, 0.0)},
 'Cl': {0: (34.96885268, 1.0),
  28: (28.02851, 0.0),
  29: (29.01411, 0.0),
  30: (30.00477, 0.0),
  31: (30.99241, 0.0),
  32: (31.98569, 0.0),
  33: (32.9774519, 0.0),
  34: (33.97376282, 0.0),
  35: (34.96885268, 0.7576),
  36: (35.96830698, 0.0),
  37: (36.96590259, 0.2424),
  38: (37.96801043, 0.0),
  39: (38.9680082, 0.0),
  40: (39.97042, 0.0),
  41: (40.97068, 0.0),
  42: (41.97325, 0.0),
  43: (42.97405, 0.0),
  44: (43.97828, 0.0),
  45: (44.98029, 0.0),
  46: (45.98421, 0.0),
  47: (46.98871, 0.0),
  48: (47.99495, 0.0),
  49: (49.00032, 0.0),
  50: (50.00784, 0.0),
  51: (51.01449, 0.0)},
 'Cm': {0: (247, 1.0),
  233: (233.05077, 0.0),
  234: (234.05016, 0.0),
  235: (235.05143, 0.0),
  236: (236.05141, 0.0),
  237: (237.0529, 0.0),
  238: (238.05303, 0.0),
  239: (239.05496, 0.0),
  240: (240.0555295, 0.0),
  241: (241.057653, 0.0),
  242: (242.0588358, 0.0),
  243: (243.0613891, 0.0),
  244: (244.0627526, 0.0),
  245: (245.0654912, 0.0),
  246: (246.0672237, 0.0),
  247: (247.070354, 0.0),
  248: (248.072349, 0.0),
  249: (249.075953, 0.0),
  250: (250.078357, 0.0),
  251: (251.082285, 0.0),
  252: (252.08487, 0.0)},
 'Cn': {0: (285, 1.0),
  277: (277.16394, 0.0),
  278: (278.16431, 0.0),
  279: (279.16655, 0.0),
  280: (280.16704, 0.0),
  281: (281.16929, 0.0),
  282: (282.16977, 0.0),
  283: (283.17179, 0.0),
  284: (284.17238, 0.0),
  285: (285.17411, 0.0)},
 'Co': {0: (58.933195, 1.0),
  47: (47.01149, 0.0),
  48: (48.00176, 0.0),
  49: (48.98972, 0.0),
  50: (49.98154, 0.0),
  51: (50.97072, 0.0),
  52: (51.96359, 0.0),
  53: (52.954219, 0.0),
  54: (53.9484596, 0.0),
  55: (54.941999, 0.0),
  56: (55.9398393, 0.0),
  57: (56.9362914, 0.0),
  58: (57.9357528, 0.0),
  59: (58.933195, 1.0),
  60: (59.9338171, 0.0),
  61: (60.9324758, 0.0),
  62: (61.934051, 0.0),
  63: (62.933612, 0.0),
  64: (63.93581, 0.0),
  65: (64.936478, 0.0),
  66: (65.93976, 0.0),
  67: (66.94089, 0.0),
  68: (67.94487, 0.0),
  69: (68.94632, 0.0),
  70: (69.951, 0.0),
  71: (70.9529, 0.0),
  72: (71.95781, 0.0),
  73: (72.96024, 0.0),
  74: (73.96538, 0.0),
  75: (74.96833, 0.0)},
 'Cr': {0: (51.9405075, 1.0),
  42: (42.00643, 0.0),
  43: (42.99771, 0.0),
  44: (43.98555, 0.0),
  45: (44.97964, 0.0),
  46: (45.968359, 0.0),
  47: (46.9629, 0.0),
  48: (47.954032, 0.0),
  49: (48.9513357, 0.0),
  50: (49.9460442, 0.04345),
  51: (50.9447674, 0.0),
  52: (51.9405075, 0.83789),
  53: (52.9406494, 0.09501),
  54: (53.9388804, 0.02365),
  55: (54.9408397, 0.0),
  56: (55.9406531, 0.0),
  57: (56.943613, 0.0),
  58: (57.94435, 0.0),
  59: (58.94859, 0.0),
  60: (59.95008, 0.0),
  61: (60.95472, 0.0),
  62: (61.95661, 0.0),
  63: (62.96186, 0.0),
  64: (63.96441, 0.0),
  65: (64.97016, 0.0),
  66: (65.97338, 0.0),
  67: (66.97955, 0.0)},
 'Cs': {0: (132.905451933, 1.0),
  112: (111.9503, 0.0),
  113: (112.94449, 0.0),
  114: (113.94145, 0.0),
  115: (114.93591, 0.0),
  116: (115.93337, 0.0),
  117: (116.92867, 0.0),
  118: (117.926559, 0.0),
  119: (118.922377, 0.0),
  120: (119.920677, 0.0),
  121: (120.917229, 0.0),
  122: (121.91611, 0.0),
  123: (122.912996, 0.0),
  124: (123.912258, 0.0),
  125: (124.909728, 0.0),
  126: (125.909452, 0.0),
  127: (126.907418, 0.0),
  128: (127.907749, 0.0),
  129: (128.906064, 0.0),
  130: (129.906709, 0.0),
  131: (130.905464, 0.0),
  132: (131.9064343, 0.0),
  133: (132.905451933, 1.0),
  134: (133.906718475, 0.0),
  135: (134.905977, 0.0),
  136: (135.9073116, 0.0),
  137: (136.9070895, 0.0),
  138: (137.911017, 0.0),
  139: (138.913364, 0.0),
  140: (139.917282, 0.0),
  141: (140.920046, 0.0),
  142: (141.924299, 0.0),
  143: (142.927352, 0.0),
  144: (143.932077, 0.0),
  145: (144.935526, 0.0),
  146: (145.94029, 0.0),
  147: (146.94416, 0.0),
  148: (147.94922, 0.0),
  149: (148.95293, 0.0),
  150: (149.95817, 0.0),
  151: (150.96219, 0.0)},
 'Cu': {0: (62.9295975, 1.0),
  52: (51.99718, 0.0),
  53: (52.98555, 0.0),
  54: (53.97671, 0.0),
  55: (54.96605, 0.0),
  56: (55.95856, 0.0),
  57: (56.949211, 0.0),
  58: (57.9445385, 0.0),
  59: (58.939498, 0.0),
  60: (59.937365, 0.0),
  61: (60.9334578, 0.0),
  62: (61.932584, 0.0),
  63: (62.9295975, 0.6915),
  64: (63.9297642, 0.0),
  65: (64.9277895, 0.3085),
  66: (65.9288688, 0.0),
  67: (66.9277303, 0.0),
  68: (67.9296109, 0.0),
  69: (68.9294293, 0.0),
  70: (69.9323923, 0.0),
  71: (70.9326768, 0.0),
  72: (71.9358203, 0.0),
  73: (72.936675, 0.0),
  74: (73.939875, 0.0),
  75: (74.9419, 0.0),
  76: (75.945275, 0.0),
  77: (76.94785, 0.0),
  78: (77.95196, 0.0),
  79: (78.95456, 0.0),
  80: (79.96087, 0.0)},
 'Db': {0: (268, 1.0),
  255: (255.1074, 0.0),
  256: (256.10813, 0.0),
  257: (257.10772, 0.0),
  258: (258.10923, 0.0),
  259: (259.10961, 0.0),
  260: (260.1113, 0.0),
  261: (261.11206, 0.0),
  262: (262.11408, 0.0),
  263: (263.11499, 0.0),
  264: (264.1174, 0.0),
  265: (265.1186, 0.0),
  266: (266.12103, 0.0),
  267: (267.12238, 0.0),
  268: (268.12545, 0.0),
  269: (269.12746, 0.0),
  270: (270.13071, 0.0)},
 'Ds': {0: (281, 1.0),
  267: (267.14434, 0.0),
  268: (268.1438, 0.0),
  269: (269.14512, 0.0),
  270: (270.14472, 0.0),
  271: (271.14606, 0.0),
  272: (272.14632, 0.0),
  273: (273.14886, 0.0),
  274: (274.14949, 0.0),
  275: (275.15218, 0.0),
  276: (276.15303, 0.0),
  277: (277.15565, 0.0),
  278: (278.15647, 0.0),
  279: (279.15886, 0.0),
  280: (280.1598, 0.0),
  281: (281.16206, 0.0)},
 'Dy': {0: (163.9291748, 1.0),
  138: (137.96249, 0.0),
  139: (138.95954, 0.0),
  140: (139.95401, 0.0),
  141: (140.95135, 0.0),
  142: (141.94637, 0.0),
  143: (142.94383, 0.0),
  144: (143.93925, 0.0),
  145: (144.93743, 0.0),
  146: (145.932845, 0.0),
  147: (146.931092, 0.0),
  148: (147.92715, 0.0),
  149: (148.927305, 0.0),
  150: (149.925585, 0.0),
  151: (150.926185, 0.0),
  152: (151.924718, 0.0),
  153: (152.925765, 0.0),
  154: (153.924424, 0.0),
  155: (154.925754, 0.0),
  156: (155.924283, 0.00056),
  157: (156.925466, 0.0),
  158: (157.924409, 0.00095),
  159: (158.9257392, 0.0),
  160: (159.9251975, 0.02329),
  161: (160.9269334, 0.18889),
  162: (161.9267984, 0.25475),
  163: (162.9287312, 0.24896),
  164: (163.9291748, 0.2826),
  165: (164.9317033, 0.0),
  166: (165.9328067, 0.0),
  167: (166.93566, 0.0),
  168: (167.93713, 0.0),
  169: (168.94031, 0.0),
  170: (169.94239, 0.0),
  171: (170.9462, 0.0),
  172: (171.94876, 0.0),
  173: (172.953, 0.0)},
 'Er': {0: (165.9302931, 1.0),
  143: (142.96634, 0.0),
  144: (143.96038, 0.0),
  145: (144.95739, 0.0),
  146: (145.952, 0.0),
  147: (146.94949, 0.0),
  148: (147.94455, 0.0),
  149: (148.94231, 0.0),
  150: (149.937914, 0.0),
  151: (150.937449, 0.0),
  152: (151.93505, 0.0),
  153: (152.935063, 0.0),
  154: (153.932783, 0.0),
  155: (154.933209, 0.0),
  156: (155.931065, 0.0),
  157: (156.93192, 0.0),
  158: (157.929893, 0.0),
  159: (158.930684, 0.0),
  160: (159.929083, 0.0),
  161: (160.929995, 0.0),
  162: (161.928778, 0.00139),
  163: (162.930033, 0.0),
  164: (163.9292, 0.01601),
  165: (164.930726, 0.0),
  166: (165.9302931, 0.33503),
  167: (166.9320482, 0.22869),
  168: (167.9323702, 0.26978),
  169: (168.9345904, 0.0),
  170: (169.9354643, 0.1491),
  171: (170.9380298, 0.0),
  172: (171.939356, 0.0),
  173: (172.9424, 0.0),
  174: (173.94423, 0.0),
  175: (174.94777, 0.0),
  176: (175.95008, 0.0),
  177: (176.95405, 0.0)},
 'Es': {0: (252, 1.0),
  240: (240.06892, 0.0),
  241: (241.06854, 0.0),
  242: (242.06975, 0.0),
  243: (243.06955, 0.0),
  244: (244.07088, 0.0),
  245: (245.07132, 0.0),
  246: (246.0729, 0.0),
  247: (247.07366, 0.0),
  248: (248.07547, 0.0),
  249: (249.07641, 0.0),
  250: (250.07861, 0.0),
  251: (251.079992, 0.0),
  252: (252.08298, 0.0),
  253: (253.0848247, 0.0),
  254: (254.088022, 0.0),
  255: (255.090273, 0.0),
  256: (256.0936, 0.0),
  257: (257.09598, 0.0),
  258: (258.09952, 0.0)},
 'Eu': {0: (152.9212303, 1.0),
  130: (129.96357, 0.0),
  131: (130.95775, 0.0),
  132: (131.95437, 0.0),
  133: (132.94924, 0.0),
  134: (133.94651, 0.0),
  135: (134.94182, 0.0),
  136: (135.9396, 0.0),
  137: (136.93557, 0.0),
  138: (137.93371, 0.0),
  139: (138.929792, 0.0),
  140: (139.92809, 0.0),
  141: (140.924931, 0.0),
  142: (141.92343, 0.0),
  143: (142.920298, 0.0),
  144: (143.918817, 0.0),
  145: (144.916265, 0.0),
  146: (145.917206, 0.0),
  147: (146.916746, 0.0),
  148: (147.918086, 0.0),
  149: (148.917931, 0.0),
  150: (149.919702, 0.0),
  151: (150.9198502, 0.4781),
  152: (151.9217445, 0.0),
  153: (152.9212303, 0.5219),
  154: (153.9229792, 0.0),
  155: (154.9228933, 0.0),
  156: (155.924752, 0.0),
  157: (156.925424, 0.0),
  158: (157.92785, 0.0),
  159: (158.929089, 0.0),
  160: (159.93197, 0.0),
  161: (160.93368, 0.0),
  162: (161.93704, 0.0),
  163: (162.93921, 0.0),
  164: (163.94299, 0.0),
  165: (164.94572, 0.0),
  166: (165.94997, 0.0),
  167: (166.95321, 0.0)},
 'F': {0: (18.99840322, 1.0),
  14: (14.03506, 0.0),
  15: (15.01801, 0.0),
  16: (16.011466, 0.0),
  17: (17.00209524, 0.0),
  18: (18.000938, 0.0),
  19: (18.99840322, 1.0),
  20: (19.99998132, 0.0),
  21: (20.999949, 0.0),
  22: (22.002999, 0.0),
  23: (23.00357, 0.0),
  24: (24.00812, 0.0),
  25: (25.0121, 0.0),
  26: (26.01962, 0.0),
  27: (27.02676, 0.0),
  28: (28.03567, 0.0),
  29: (29.04326, 0.0),
  30: (30.0525, 0.0),
  31: (31.06043, 0.0)},
 'Fe': {0: (55.9349375, 1.0),
  45: (45.01458, 0.0),
  46: (46.00081, 0.0),
  47: (46.99289, 0.0),
  48: (47.9805, 0.0),
  49: (48.97361, 0.0),
  50: (49.96299, 0.0),
  51: (50.95682, 0.0),
  52: (51.948114, 0.0),
  53: (52.9453079, 0.0),
  54: (53.9396105, 0.05845),
  55: (54.9382934, 0.0),
  56: (55.9349375, 0.91754),
  57: (56.935394, 0.02119),
  58: (57.9332756, 0.00282),
  59: (58.9348755, 0.0),
  60: (59.934072, 0.0),
  61: (60.936745, 0.0),
  62: (61.936767, 0.0),
  63: (62.94037, 0.0),
  64: (63.9412, 0.0),
  65: (64.94538, 0.0),
  66: (65.94678, 0.0),
  67: (66.95095, 0.0),
  68: (67.9537, 0.0),
  69: (68.95878, 0.0),
  70: (69.96146, 0.0),
  71: (70.96672, 0.0),
  72: (71.96962, 0.0)},
 'Fm': {0: (257, 1.0),
  242: (242.07343, 0.0),
  243: (243.07435, 0.0),
  244: (244.07408, 0.0),
  245: (245.07539, 0.0),
  246: (246.0753, 0.0),
  247: (247.07685, 0.0),
  248: (248.077195, 0.0),
  249: (249.07903, 0.0),
  250: (250.079521, 0.0),
  251: (251.081575, 0.0),
  252: (252.082467, 0.0),
  253: (253.085185, 0.0),
  254: (254.0868542, 0.0),
  255: (255.089962, 0.0),
  256: (256.091773, 0.0),
  257: (257.095105, 0.0),
  258: (258.09708, 0.0),
  259: (259.1006, 0.0),
  260: (260.10268, 0.0)},
 'Fr': {0: (223, 1.0),
  199: (199.00726, 0.0),
  200: (200.00657, 0.0),
  201: (201.00386, 0.0),
  202: (202.00337, 0.0),
  203: (203.000925, 0.0),
  204: (204.000653, 0.0),
  205: (204.998594, 0.0),
  206: (205.99867, 0.0),
  207: (206.99695, 0.0),
  208: (207.99714, 0.0),
  209: (208.995954, 0.0),
  210: (209.996408, 0.0),
  211: (210.995537, 0.0),
  212: (211.996202, 0.0),
  213: (212.996189, 0.0),
  214: (213.998971, 0.0),
  215: (215.000341, 0.0),
  216: (216.003198, 0.0),
  217: (217.004632, 0.0),
  218: (218.007578, 0.0),
  219: (219.009252, 0.0),
  220: (220.012327, 0.0),
  221: (221.014255, 0.0),
  222: (222.017552, 0.0),
  223: (223.0197359, 0.0),
  224: (224.02325, 0.0),
  225: (225.02557, 0.0),
  226: (226.02939, 0.0),
  227: (227.03184, 0.0),
  228: (228.03573, 0.0),
  229: (229.03845, 0.0),
  230: (230.04251, 0.0),
  231: (231.04544, 0.0),
  232: (232.04977, 0.0)},
 'Ga': {0: (68.9255736, 1.0),
  56: (55.99491, 0.0),
  57: (56.98293, 0.0),
  58: (57.97425, 0.0),
  59: (58.96337, 0.0),
  60: (59.95706, 0.0),
  61: (60.94945, 0.0),
  62: (61.944175, 0.0),
  63: (62.9392942, 0.0),
  64: (63.9368387, 0.0),
  65: (64.9327348, 0.0),
  66: (65.931589, 0.0),
  67: (66.9282017, 0.0),
  68: (67.9279801, 0.0),
  69: (68.9255736, 0.60108),
  70: (69.926022, 0.0),
  71: (70.9247013, 0.39892),
  72: (71.9263663, 0.0),
  73: (72.9251747, 0.0),
  74: (73.926946, 0.0),
  75: (74.9265002, 0.0),
  76: (75.9288276, 0.0),
  77: (76.9291543, 0.0),
  78: (77.9316082, 0.0),
  79: (78.93289, 0.0),
  80: (79.93652, 0.0),
  81: (80.93775, 0.0),
  82: (81.94299, 0.0),
  83: (82.94698, 0.0),
  84: (83.95265, 0.0),
  85: (84.957, 0.0),
  86: (85.96312, 0.0)},
 'Gd': {0: (157.9241039, 1.0),
  134: (133.95537, 0.0),
  135: (134.95257, 0.0),
  136: (135.94734, 0.0),
  137: (136.94502, 0.0),
  138: (137.94012, 0.0),
  139: (138.93824, 0.0),
  140: (139.93367, 0.0),
  141: (140.932126, 0.0),
  142: (141.92812, 0.0),
  143: (142.92675, 0.0),
  144: (143.92296, 0.0),
  145: (144.921709, 0.0),
  146: (145.918311, 0.0),
  147: (146.919094, 0.0),
  148: (147.918115, 0.0),
  149: (148.919341, 0.0),
  150: (149.918659, 0.0),
  151: (150.920348, 0.0),
  152: (151.919791, 0.002),
  153: (152.9217495, 0.0),
  154: (153.9208656, 0.0218),
  155: (154.922622, 0.148),
  156: (155.9221227, 0.2047),
  157: (156.9239601, 0.1565),
  158: (157.9241039, 0.2484),
  159: (158.9263887, 0.0),
  160: (159.9270541, 0.2186),
  161: (160.9296692, 0.0),
  162: (161.930985, 0.0),
  163: (162.93399, 0.0),
  164: (163.93586, 0.0),
  165: (164.93938, 0.0),
  166: (165.9416, 0.0),
  167: (166.94557, 0.0),
  168: (167.94836, 0.0),
  169: (168.95287, 0.0)},
 'Ge': {0: (73.9211778, 1.0),
  58: (57.99101, 0.0),
  59: (58.98175, 0.0),
  60: (59.97019, 0.0),
  61: (60.96379, 0.0),
  62: (61.95465, 0.0),
  63: (62.94964, 0.0),
  64: (63.94165, 0.0),
  65: (64.93944, 0.0),
  66: (65.93384, 0.0),
  67: (66.932734, 0.0),
  68: (67.928094, 0.0),
  69: (68.9279645, 0.0),
  70: (69.9242474, 0.2038),
  71: (70.924951, 0.0),
  72: (71.9220758, 0.2731),
  73: (72.9234589, 0.0776),
  74: (73.9211778, 0.3672),
  75: (74.9228589, 0.0),
  76: (75.9214026, 0.0783),
  77: (76.9235486, 0.0),
  78: (77.922853, 0.0),
  79: (78.9254, 0.0),
  80: (79.92537, 0.0),
  81: (80.92882, 0.0),
  82: (81.92955, 0.0),
  83: (82.93462, 0.0),
  84: (83.93747, 0.0),
  85: (84.94303, 0.0),
  86: (85.94649, 0.0),
  87: (86.95251, 0.0),
  88: (87.95691, 0.0),
  89: (88.96383, 0.0)},
 'H': {0: (1.00782503207, 1.0),
  1: (1.00782503207, 0.999885),
  2: (2.0141017778, 0.000115),
  3: (3.0160492777, 0.0),
  4: (4.02781, 0.0),
  5: (5.03531, 0.0),
  6: (6.04494, 0.0),
  7: (7.05275, 0.0)},
 'H+': {0: (1.00727646677, 1.0), 1: (1.00727646677, 1.0)},
 'He': {0: (4.00260325415, 1.0),
  3: (3.0160293191, 1.34e-06),
  4: (4.00260325415, 0.99999866),
  5: (5.01222, 0.0),
  6: (6.0188891, 0.0),
  7: (7.028021, 0.0),
  8: (8.033922, 0.0),
  9: (9.04395, 0.0),
  10: (10.0524, 0.0)},
 'Hf': {0: (179.94655, 1.0),
  153: (152.97069, 0.0),
  154: (153.96486, 0.0),
  155: (154.96339, 0.0),
  156: (155.95936, 0.0),
  157: (156.9584, 0.0),
  158: (157.954799, 0.0),
  159: (158.953995, 0.0),
  160: (159.950684, 0.0),
  161: (160.950275, 0.0),
  162: (161.94721, 0.0),
  163: (162.94709, 0.0),
  164: (163.944367, 0.0),
  165: (164.94457, 0.0),
  166: (165.94218, 0.0),
  167: (166.9426, 0.0),
  168: (167.94057, 0.0),
  169: (168.94126, 0.0),
  170: (169.93961, 0.0),
  171: (170.94049, 0.0),
  172: (171.939448, 0.0),
  173: (172.94051, 0.0),
  174: (173.940046, 0.0016),
  175: (174.941509, 0.0),
  176: (175.9414086, 0.0526),
  177: (176.9432207, 0.186),
  178: (177.9436988, 0.2728),
  179: (178.9458161, 0.1362),
  180: (179.94655, 0.3508),
  181: (180.9491012, 0.0),
  182: (181.950554, 0.0),
  183: (182.95353, 0.0),
  184: (183.95545, 0.0),
  185: (184.95882, 0.0),
  186: (185.96089, 0.0),
  187: (186.96459, 0.0),
  188: (187.96685, 0.0)},
 'Hg': {0: (201.970643, 1.0),
  171: (171.00376, 0.0),
  172: (171.99883, 0.0),
  173: (172.99724, 0.0),
  174: (173.992864, 0.0),
  175: (174.99142, 0.0),
  176: (175.987355, 0.0),
  177: (176.98628, 0.0),
  178: (177.982483, 0.0),
  179: (178.981834, 0.0),
  180: (179.978266, 0.0),
  181: (180.977819, 0.0),
  182: (181.97469, 0.0),
  183: (182.97445, 0.0),
  184: (183.971713, 0.0),
  185: (184.971899, 0.0),
  186: (185.969362, 0.0),
  187: (186.969814, 0.0),
  188: (187.967577, 0.0),
  189: (188.96819, 0.0),
  190: (189.966322, 0.0),
  191: (190.967157, 0.0),
  192: (191.965634, 0.0),
  193: (192.966665, 0.0),
  194: (193.965439, 0.0),
  195: (194.96672, 0.0),
  196: (195.965833, 0.0015),
  197: (196.967213, 0.0),
  198: (197.966769, 0.0997),
  199: (198.9682799, 0.1687),
  200: (199.968326, 0.231),
  201: (200.9703023, 0.1318),
  202: (201.970643, 0.2986),
  203: (202.9728725, 0.0),
  204: (203.9734939, 0.0687),
  205: (204.976073, 0.0),
  206: (205.977514, 0.0),
  207: (206.98259, 0.0),
  208: (207.98594, 0.0),
  209: (208.99104, 0.0),
  210: (209.99451, 0.0)},
 'Ho': {0: (164.9303221, 1.0),
  140: (139.96854, 0.0),
  141: (140.9631, 0.0),
  142: (141.95977, 0.0),
  143: (142.95461, 0.0),
  144: (143.95148, 0.0),
  145: (144.9472, 0.0),
  146: (145.94464, 0.0),
  147: (146.94006, 0.0),
  148: (147.93772, 0.0),
  149: (148.933775, 0.0),
  150: (149.933496, 0.0),
  151: (150.931688, 0.0),
  152: (151.931714, 0.0),
  153: (152.930199, 0.0),
  154: (153.930602, 0.0),
  155: (154.929103, 0.0),
  156: (155.92984, 0.0),
  157: (156.928256, 0.0),
  158: (157.928941, 0.0),
  159: (158.927712, 0.0),
  160: (159.928729, 0.0),
  161: (160.927855, 0.0),
  162: (161.929096, 0.0),
  163: (162.9287339, 0.0),
  164: (163.9302335, 0.0),
  165: (164.9303221, 1.0),
  166: (165.9322842, 0.0),
  167: (166.933133, 0.0),
  168: (167.93552, 0.0),
  169: (168.936872, 0.0),
  170: (169.93962, 0.0),
  171: (170.94147, 0.0),
  172: (171.94482, 0.0),
  173: (172.94729, 0.0),
  174: (173.95115, 0.0),
  175: (174.95405, 0.0)},
 'Hs': {0: (270, 1.0),
  263: (263.12856, 0.0),
  264: (264.12839, 0.0),
  265: (265.13009, 0.0),
  266: (266.1301, 0.0),
  267: (267.13179, 0.0),
  268: (268.13216, 0.0),
  269: (269.13406, 0.0),
  270: (270.13465, 0.0),
  271: (271.13766, 0.0),
  272: (272.13905, 0.0),
  273: (273.14199, 0.0),
  274: (274.14313, 0.0),
  275: (275.14595, 0.0),
  276: (276.14721, 0.0),
  277: (277.14984, 0.0)},
 'I': {0: (126.904473, 1.0),
  108: (107.94348, 0.0),
  109: (108.93815, 0.0),
  110: (109.93524, 0.0),
  111: (110.93028, 0.0),
  112: (111.92797, 0.0),
  113: (112.92364, 0.0),
  114: (113.92185, 0.0),
  115: (114.91805, 0.0),
  116: (115.91681, 0.0),
  117: (116.91365, 0.0),
  118: (117.913074, 0.0),
  119: (118.91007, 0.0),
  120: (119.910048, 0.0),
  121: (120.907367, 0.0),
  122: (121.907589, 0.0),
  123: (122.905589, 0.0),
  124: (123.9062099, 0.0),
  125: (124.9046302, 0.0),
  126: (125.905624, 0.0),
  127: (126.904473, 1.0),
  128: (127.905809, 0.0),
  129: (128.904988, 0.0),
  130: (129.906674, 0.0),
  131: (130.9061246, 0.0),
  132: (131.907997, 0.0),
  133: (132.907797, 0.0),
  134: (133.909744, 0.0),
  135: (134.910048, 0.0),
  136: (135.91465, 0.0),
  137: (136.917871, 0.0),
  138: (137.92235, 0.0),
  139: (138.9261, 0.0),
  140: (139.931, 0.0),
  141: (140.93503, 0.0),
  142: (141.94018, 0.0),
  143: (142.94456, 0.0),
  144: (143.94999, 0.0)},
 'In': {0: (114.903878, 1.0),
  97: (96.94954, 0.0),
  98: (97.94214, 0.0),
  99: (98.93422, 0.0),
  100: (99.93111, 0.0),
  101: (100.92634, 0.0),
  102: (101.92409, 0.0),
  103: (102.919914, 0.0),
  104: (103.9183, 0.0),
  105: (104.914674, 0.0),
  106: (105.913465, 0.0),
  107: (106.910295, 0.0),
  108: (107.909698, 0.0),
  109: (108.907151, 0.0),
  110: (109.907165, 0.0),
  111: (110.905103, 0.0),
  112: (111.905532, 0.0),
  113: (112.904058, 0.0429),
  114: (113.904914, 0.0),
  115: (114.903878, 0.9571),
  116: (115.90526, 0.0),
  117: (116.904514, 0.0),
  118: (117.906354, 0.0),
  119: (118.905845, 0.0),
  120: (119.90796, 0.0),
  121: (120.907846, 0.0),
  122: (121.91028, 0.0),
  123: (122.910438, 0.0),
  124: (123.91318, 0.0),
  125: (124.9136, 0.0),
  126: (125.91646, 0.0),
  127: (126.91735, 0.0),
  128: (127.92017, 0.0),
  129: (128.9217, 0.0),
  130: (129.92497, 0.0),
  131: (130.92685, 0.0),
  132: (131.93299, 0.0),
  133: (132.93781, 0.0),
  134: (133.94415, 0.0),
  135: (134.94933, 0.0)},
 'Ir': {0: (192.9629264, 1.0),
  164: (163.9922, 0.0),
  165: (164.98752, 0.0),
  166: (165.98582, 0.0),
  167: (166.981665, 0.0),
  168: (167.97988, 0.0),
  169: (168.976295, 0.0),
  170: (169.97497, 0.0),
  171: (170.97163, 0.0),
  172: (171.97046, 0.0),
  173: (172.967502, 0.0),
  174: (173.966861, 0.0),
  175: (174.964113, 0.0),
  176: (175.963649, 0.0),
  177: (176.961302, 0.0),
  178: (177.961082, 0.0),
  179: (178.959122, 0.0),
  180: (179.959229, 0.0),
  181: (180.957625, 0.0),
  182: (181.958076, 0.0),
  183: (182.956846, 0.0),
  184: (183.95748, 0.0),
  185: (184.9567, 0.0),
  186: (185.957946, 0.0),
  187: (186.957363, 0.0),
  188: (187.958853, 0.0),
  189: (188.958719, 0.0),
  190: (189.960546, 0.0),
  191: (190.960594, 0.373),
  192: (191.962605, 0.0),
  193: (192.9629264, 0.627),
  194: (193.9650784, 0.0),
  195: (194.9659796, 0.0),
  196: (195.9684, 0.0),
  197: (196.969653, 0.0),
  198: (197.97228, 0.0),
  199: (198.9738, 0.0)},
 'K': {0: (38.96370668, 1.0),
  32: (32.02192, 0.0),
  33: (33.00726, 0.0),
  34: (33.99841, 0.0),
  35: (34.98801, 0.0),
  36: (35.981292, 0.0),
  37: (36.97337589, 0.0),
  38: (37.9690812, 0.0),
  39: (38.96370668, 0.932581),
  40: (39.96399848, 0.000117),
  41: (40.96182576, 0.067302),
  42: (41.96240281, 0.0),
  43: (42.960716, 0.0),
  44: (43.96156, 0.0),
  45: (44.960699, 0.0),
  46: (45.961977, 0.0),
  47: (46.961678, 0.0),
  48: (47.965514, 0.0),
  49: (48.96745, 0.0),
  50: (49.97278, 0.0),
  51: (50.97638, 0.0),
  52: (51.98261, 0.0),
  53: (52.98712, 0.0),
  54: (53.9942, 0.0),
  55: (54.99971, 0.0)},
 'Kr': {0: (83.911507, 1.0),
  69: (68.96518, 0.0),
  70: (69.95526, 0.0),
  71: (70.94963, 0.0),
  72: (71.942092, 0.0),
  73: (72.939289, 0.0),
  74: (73.9330844, 0.0),
  75: (74.930946, 0.0),
  76: (75.92591, 0.0),
  77: (76.92467, 0.0),
  78: (77.9203648, 0.00355),
  79: (78.920082, 0.0),
  80: (79.916379, 0.02286),
  81: (80.916592, 0.0),
  82: (81.9134836, 0.11593),
  83: (82.914136, 0.115),
  84: (83.911507, 0.56987),
  85: (84.9125273, 0.0),
  86: (85.91061073, 0.17279),
  87: (86.91335486, 0.0),
  88: (87.914447, 0.0),
  89: (88.91763, 0.0),
  90: (89.919517, 0.0),
  91: (90.92345, 0.0),
  92: (91.926156, 0.0),
  93: (92.93127, 0.0),
  94: (93.93436, 0.0),
  95: (94.93984, 0.0),
  96: (95.94307, 0.0),
  97: (96.94856, 0.0),
  98: (97.95191, 0.0),
  99: (98.9576, 0.0),
  100: (99.96114, 0.0)},
 'La': {0: (138.9063533, 1.0),
  117: (116.95007, 0.0),
  118: (117.94673, 0.0),
  119: (118.94099, 0.0),
  120: (119.93807, 0.0),
  121: (120.93301, 0.0),
  122: (121.93071, 0.0),
  123: (122.92624, 0.0),
  124: (123.92457, 0.0),
  125: (124.920816, 0.0),
  126: (125.91951, 0.0),
  127: (126.916375, 0.0),
  128: (127.91559, 0.0),
  129: (128.912693, 0.0),
  130: (129.912369, 0.0),
  131: (130.91007, 0.0),
  132: (131.9101, 0.0),
  133: (132.90822, 0.0),
  134: (133.908514, 0.0),
  135: (134.906977, 0.0),
  136: (135.90764, 0.0),
  137: (136.906494, 0.0),
  138: (137.907112, 0.0009),
  139: (138.9063533, 0.9991),
  140: (139.9094776, 0.0),
  141: (140.910962, 0.0),
  142: (141.914079, 0.0),
  143: (142.916063, 0.0),
  144: (143.9196, 0.0),
  145: (144.92165, 0.0),
  146: (145.92579, 0.0),
  147: (146.92824, 0.0),
  148: (147.93223, 0.0),
  149: (148.93473, 0.0),
  150: (149.93877, 0.0),
  151: (150.94172, 0.0),
  152: (151.94625, 0.0),
  153: (152.94962, 0.0),
  154: (153.9545, 0.0),
  155: (154.95835, 0.0)},
 'Li': {0: (7.01600455, 1.0),
  3: (3.03078, 0.0),
  4: (4.02719, 0.0),
  5: (5.01254, 0.0),
  6: (6.015122795, 0.0759),
  7: (7.01600455, 0.9241),
  8: (8.02248736, 0.0),
  9: (9.0267895, 0.0),
  10: (10.035481, 0.0),
  11: (11.043798, 0.0),
  12: (12.05378, 0.0)},
 'Lr': {0: (262, 1.0),
  251: (251.09436, 0.0),
  252: (252.09537, 0.0),
  253: (253.09521, 0.0),
  254: (254.09645, 0.0),
  255: (255.09668, 0.0),
  256: (256.09863, 0.0),
  257: (257.09956, 0.0),
  258: (258.10181, 0.0),
  259: (259.1029, 0.0),
  260: (260.1055, 0.0),
  261: (261.10688, 0.0),
  262: (262.10963, 0.0),
  263: (263.11129, 0.0),
  264: (264.11404, 0.0),
  265: (265.11584, 0.0),
  266: (266.11931, 0.0)},
 'Lu': {0: (174.9407718, 1.0),
  150: (149.97323, 0.0),
  151: (150.96758, 0.0),
  152: (151.96412, 0.0),
  153: (152.95877, 0.0),
  154: (153.95752, 0.0),
  155: (154.954316, 0.0),
  156: (155.95303, 0.0),
  157: (156.950098, 0.0),
  158: (157.949313, 0.0),
  159: (158.94663, 0.0),
  160: (159.94603, 0.0),
  161: (160.94357, 0.0),
  162: (161.94328, 0.0),
  163: (162.94118, 0.0),
  164: (163.94134, 0.0),
  165: (164.939407, 0.0),
  166: (165.93986, 0.0),
  167: (166.93827, 0.0),
  168: (167.93874, 0.0),
  169: (168.937651, 0.0),
  170: (169.938475, 0.0),
  171: (170.9379131, 0.0),
  172: (171.939086, 0.0),
  173: (172.9389306, 0.0),
  174: (173.9403375, 0.0),
  175: (174.9407718, 0.9741),
  176: (175.9426863, 0.0259),
  177: (176.9437581, 0.0),
  178: (177.945955, 0.0),
  179: (178.947327, 0.0),
  180: (179.94988, 0.0),
  181: (180.95197, 0.0),
  182: (181.95504, 0.0),
  183: (182.95757, 0.0),
  184: (183.96091, 0.0)},
 'Md': {0: (258, 1.0),
  245: (245.08083, 0.0),
  246: (246.08189, 0.0),
  247: (247.08164, 0.0),
  248: (248.08282, 0.0),
  249: (249.08301, 0.0),
  250: (250.08442, 0.0),
  251: (251.08484, 0.0),
  252: (252.08656, 0.0),
  253: (253.08728, 0.0),
  254: (254.08966, 0.0),
  255: (255.091083, 0.0),
  256: (256.09406, 0.0),
  257: (257.095541, 0.0),
  258: (258.098431, 0.0),
  259: (259.10051, 0.0),
  260: (260.10365, 0.0),
  261: (261.10572, 0.0),
  262: (262.10887, 0.0)},
 'Mg': {0: (23.9850417, 1.0),
  19: (19.03547, 0.0),
  20: (20.018863, 0.0),
  21: (21.011713, 0.0),
  22: (21.9995738, 0.0),
  23: (22.9941237, 0.0),
  24: (23.9850417, 0.7899),
  25: (24.98583692, 0.1),
  26: (25.982592929, 0.1101),
  27: (26.98434059, 0.0),
  28: (27.9838768, 0.0),
  29: (28.9886, 0.0),
  30: (29.990434, 0.0),
  31: (30.996546, 0.0),
  32: (31.998975, 0.0),
  33: (33.005254, 0.0),
  34: (34.00946, 0.0),
  35: (35.01734, 0.0),
  36: (36.023, 0.0),
  37: (37.0314, 0.0),
  38: (38.03757, 0.0),
  39: (39.04677, 0.0),
  40: (40.05393, 0.0)},
 'Mn': {0: (54.9380451, 1.0),
  44: (44.00687, 0.0),
  45: (44.99451, 0.0),
  46: (45.98672, 0.0),
  47: (46.9761, 0.0),
  48: (47.96852, 0.0),
  49: (48.959618, 0.0),
  50: (49.9542382, 0.0),
  51: (50.9482108, 0.0),
  52: (51.9455655, 0.0),
  53: (52.9412901, 0.0),
  54: (53.9403589, 0.0),
  55: (54.9380451, 1.0),
  56: (55.9389049, 0.0),
  57: (56.9382854, 0.0),
  58: (57.93998, 0.0),
  59: (58.94044, 0.0),
  60: (59.94291, 0.0),
  61: (60.94465, 0.0),
  62: (61.94843, 0.0),
  63: (62.95024, 0.0),
  64: (63.95425, 0.0),
  65: (64.95634, 0.0),
  66: (65.96108, 0.0),
  67: (66.96414, 0.0),
  68: (67.9693, 0.0),
  69: (68.97284, 0.0)},
 'Mo': {0: (97.9054082, 1.0),
  83: (82.94874, 0.0),
  84: (83.94009, 0.0),
  85: (84.93655, 0.0),
  86: (85.9307, 0.0),
  87: (86.92733, 0.0),
  88: (87.921953, 0.0),
  89: (88.91948, 0.0),
  90: (89.913937, 0.0),
  91: (90.91175, 0.0),
  92: (91.906811, 0.1477),
  93: (92.906813, 0.0),
  94: (93.9050883, 0.0923),
  95: (94.9058421, 0.159),
  96: (95.9046795, 0.1668),
  97: (96.9060215, 0.0956),
  98: (97.9054082, 0.2419),
  99: (98.9077119, 0.0),
  100: (99.907477, 0.0967),
  101: (100.910347, 0.0),
  102: (101.910297, 0.0),
  103: (102.91321, 0.0),
  104: (103.91376, 0.0),
  105: (104.91697, 0.0),
  106: (105.918137, 0.0),
  107: (106.92169, 0.0),
  108: (107.92345, 0.0),
  109: (108.92781, 0.0),
  110: (109.92973, 0.0),
  111: (110.93441, 0.0),
  112: (111.93684, 0.0),
  113: (112.94188, 0.0),
  114: (113.94492, 0.0),
  115: (114.95029, 0.0)},
 'Mt': {0: (276, 1.0),
  265: (265.13615, 0.0),
  266: (266.1373, 0.0),
  267: (267.13731, 0.0),
  268: (268.13873, 0.0),
  269: (269.13906, 0.0),
  270: (270.14066, 0.0),
  271: (271.14114, 0.0),
  272: (272.14374, 0.0),
  273: (273.14491, 0.0),
  274: (274.14749, 0.0),
  275: (275.14865, 0.0),
  276: (276.15116, 0.0),
  277: (277.15242, 0.0),
  278: (278.15481, 0.0),
  279: (279.15619, 0.0)},
 'N': {0: (14.0030740048, 1.0),
  10: (10.04165, 0.0),
  11: (11.02609, 0.0),
  12: (12.0186132, 0.0),
  13: (13.00573861, 0.0),
  14: (14.0030740048, 0.99636),
  15: (15.0001088982, 0.00364),
  16: (16.0061017, 0.0),
  17: (17.00845, 0.0),
  18: (18.014079, 0.0),
  19: (19.017029, 0.0),
  20: (20.02337, 0.0),
  21: (21.02711, 0.0),
  22: (22.03439, 0.0),
  23: (23.04122, 0.0),
  24: (24.05104, 0.0),
  25: (25.06066, 0.0)},
 'Na': {0: (22.9897692809, 1.0),
  18: (18.02597, 0.0),
  19: (19.013877, 0.0),
  20: (20.007351, 0.0),
  21: (20.9976552, 0.0),
  22: (21.9944364, 0.0),
  23: (22.9897692809, 1.0),
  24: (23.99096278, 0.0),
  25: (24.989954, 0.0),
  26: (25.992633, 0.0),
  27: (26.994077, 0.0),
  28: (27.998938, 0.0),
  29: (29.002861, 0.0),
  30: (30.008976, 0.0),
  31: (31.01359, 0.0),
  32: (32.02047, 0.0),
  33: (33.02672, 0.0),
  34: (34.03517, 0.0),
  35: (35.04249, 0.0),
  36: (36.05148, 0.0),
  37: (37.05934, 0.0)},
 'Nb': {0: (92.9063781, 1.0),
  81: (80.94903, 0.0),
  82: (81.94313, 0.0),
  83: (82.93671, 0.0),
  84: (83.93357, 0.0),
  85: (84.92791, 0.0),
  86: (85.92504, 0.0),
  87: (86.92036, 0.0),
  88: (87.91833, 0.0),
  89: (88.913418, 0.0),
  90: (89.911265, 0.0),
  91: (90.906996, 0.0),
  92: (91.907194, 0.0),
  93: (92.9063781, 1.0),
  94: (93.9072839, 0.0),
  95: (94.9068358, 0.0),
  96: (95.908101, 0.0),
  97: (96.9080986, 0.0),
  98: (97.910328, 0.0),
  99: (98.911618, 0.0),
  100: (99.914182, 0.0),
  101: (100.915252, 0.0),
  102: (101.91804, 0.0),
  103: (102.91914, 0.0),
  104: (103.92246, 0.0),
  105: (104.92394, 0.0),
  106: (105.92797, 0.0),
  107: (106.93031, 0.0),
  108: (107.93484, 0.0),
  109: (108.93763, 0.0),
  110: (109.94244, 0.0),
  111: (110.94565, 0.0),
  112: (111.95083, 0.0),
  113: (112.9547, 0.0)},
 'Nd': {0: (141.9077233, 1.0),
  124: (123.95223, 0.0),
  125: (124.94888, 0.0),
  126: (125.94322, 0.0),
  127: (126.9405, 0.0),
  128: (127.93539, 0.0),
  129: (128.93319, 0.0),
  130: (129.92851, 0.0),
  131: (130.92725, 0.0),
  132: (131.923321, 0.0),
  133: (132.92235, 0.0),
  134: (133.91879, 0.0),
  135: (134.918181, 0.0),
  136: (135.914976, 0.0),
  137: (136.914567, 0.0),
  138: (137.91195, 0.0),
  139: (138.911978, 0.0),
  140: (139.90955, 0.0),
  141: (140.90961, 0.0),
  142: (141.9077233, 0.272),
  143: (142.9098143, 0.122),
  144: (143.9100873, 0.238),
  145: (144.9125736, 0.083),
  146: (145.9131169, 0.172),
  147: (146.9161004, 0.0),
  148: (147.916893, 0.057),
  149: (148.920149, 0.0),
  150: (149.920891, 0.056),
  151: (150.923829, 0.0),
  152: (151.924682, 0.0),
  153: (152.927698, 0.0),
  154: (153.92948, 0.0),
  155: (154.93293, 0.0),
  156: (155.93502, 0.0),
  157: (156.93903, 0.0),
  158: (157.9416, 0.0),
  159: (158.94609, 0.0),
  160: (159.94909, 0.0),
  161: (160.95388, 0.0)},
 'Ne': {0: (19.9924401754, 1.0),
  16: (16.025761, 0.0),
  17: (17.017672, 0.0),
  18: (18.0057082, 0.0),
  19: (19.0018802, 0.0),
  20: (19.9924401754, 0.9048),
  21: (20.99384668, 0.0027),
  22: (21.991385114, 0.0925),
  23: (22.9944669, 0.0),
  24: (23.9936108, 0.0),
  25: (24.997737, 0.0),
  26: (26.000461, 0.0),
  27: (27.00759, 0.0),
  28: (28.01207, 0.0),
  29: (29.01939, 0.0),
  30: (30.0248, 0.0),
  31: (31.03311, 0.0),
  32: (32.04002, 0.0),
  33: (33.04938, 0.0),
  34: (34.05703, 0.0)},
 'Ni': {0: (57.9353429, 1.0),
  48: (48.01975, 0.0),
  49: (49.00966, 0.0),
  50: (49.99593, 0.0),
  51: (50.98772, 0.0),
  52: (51.97568, 0.0),
  53: (52.96847, 0.0),
  54: (53.95791, 0.0),
  55: (54.95133, 0.0),
  56: (55.942132, 0.0),
  57: (56.9397935, 0.0),
  58: (57.9353429, 0.680769),
  59: (58.9343467, 0.0),
  60: (59.9307864, 0.262231),
  61: (60.931056, 0.011399),
  62: (61.9283451, 0.036345),
  63: (62.9296694, 0.0),
  64: (63.927966, 0.009256),
  65: (64.9300843, 0.0),
  66: (65.9291393, 0.0),
  67: (66.931569, 0.0),
  68: (67.931869, 0.0),
  69: (68.93561, 0.0),
  70: (69.9365, 0.0),
  71: (70.94074, 0.0),
  72: (71.94209, 0.0),
  73: (72.94647, 0.0),
  74: (73.94807, 0.0),
  75: (74.95287, 0.0),
  76: (75.95533, 0.0),
  77: (76.96055, 0.0),
  78: (77.96318, 0.0)},
 'No': {0: (259, 1.0),
  248: (248.0866, 0.0),
  249: (249.08783, 0.0),
  250: (250.08751, 0.0),
  251: (251.08901, 0.0),
  252: (252.088977, 0.0),
  253: (253.09068, 0.0),
  254: (254.090955, 0.0),
  255: (255.093241, 0.0),
  256: (256.094283, 0.0),
  257: (257.096877, 0.0),
  258: (258.09821, 0.0),
  259: (259.10103, 0.0),
  260: (260.10264, 0.0),
  261: (261.10575, 0.0),
  262: (262.1073, 0.0),
  263: (263.11055, 0.0),
  264: (264.11235, 0.0)},
 'Np': {0: (237, 1.0),
  225: (225.03391, 0.0),
  226: (226.03515, 0.0),
  227: (227.03496, 0.0),
  228: (228.03618, 0.0),
  229: (229.03626, 0.0),
  230: (230.03783, 0.0),
  231: (231.03825, 0.0),
  232: (232.04011, 0.0),
  233: (233.04074, 0.0),
  234: (234.042895, 0.0),
  235: (235.0440633, 0.0),
  236: (236.04657, 0.0),
  237: (237.0481734, 0.0),
  238: (238.0509464, 0.0),
  239: (239.052939, 0.0),
  240: (240.056162, 0.0),
  241: (241.05825, 0.0),
  242: (242.06164, 0.0),
  243: (243.06428, 0.0),
  244: (244.06785, 0.0)},
 'O': {0: (15.99491461956, 1.0),
  12: (12.034405, 0.0),
  13: (13.024812, 0.0),
  14: (14.00859625, 0.0),
  15: (15.0030656, 0.0),
  16: (15.99491461956, 0.99757),
  17: (16.9991317, 0.00038),
  18: (17.999161, 0.00205),
  19: (19.00358, 0.0),
  20: (20.0040767, 0.0),
  21: (21.008656, 0.0),
  22: (22.00997, 0.0),
  23: (23.01569, 0.0),
  24: (24.02047, 0.0),
  25: (25.02946, 0.0),
  26: (26.03834, 0.0),
  27: (27.04826, 0.0),
  28: (28.05781, 0.0)},
 'Os': {0: (191.9614807, 1.0),
  162: (161.98443, 0.0),
  163: (162.98269, 0.0),
  164: (163.97804, 0.0),
  165: (164.97676, 0.0),
  166: (165.972691, 0.0),
  167: (166.97155, 0.0),
  168: (167.967804, 0.0),
  169: (168.967019, 0.0),
  170: (169.963577, 0.0),
  171: (170.963185, 0.0),
  172: (171.960023, 0.0),
  173: (172.959808, 0.0),
  174: (173.957062, 0.0),
  175: (174.956946, 0.0),
  176: (175.95481, 0.0),
  177: (176.954965, 0.0),
  178: (177.953251, 0.0),
  179: (178.953816, 0.0),
  180: (179.952379, 0.0),
  181: (180.95324, 0.0),
  182: (181.95211, 0.0),
  183: (182.95313, 0.0),
  184: (183.9524891, 0.0002),
  185: (184.9540423, 0.0),
  186: (185.9538382, 0.0159),
  187: (186.9557505, 0.0196),
  188: (187.9558382, 0.1324),
  189: (188.9581475, 0.1615),
  190: (189.958447, 0.2626),
  191: (190.9609297, 0.0),
  192: (191.9614807, 0.4078),
  193: (192.9641516, 0.0),
  194: (193.9651821, 0.0),
  195: (194.96813, 0.0),
  196: (195.96964, 0.0)},
 'P': {0: (30.97376163, 1.0),
  24: (24.03435, 0.0),
  25: (25.02026, 0.0),
  26: (26.01178, 0.0),
  27: (26.99923, 0.0),
  28: (27.992315, 0.0),
  29: (28.9818006, 0.0),
  30: (29.9783138, 0.0),
  31: (30.97376163, 1.0),
  32: (31.97390727, 0.0),
  33: (32.9717255, 0.0),
  34: (33.973636, 0.0),
  35: (34.9733141, 0.0),
  36: (35.97826, 0.0),
  37: (36.97961, 0.0),
  38: (37.98416, 0.0),
  39: (38.98618, 0.0),
  40: (39.9913, 0.0),
  41: (40.99434, 0.0),
  42: (42.00101, 0.0),
  43: (43.00619, 0.0),
  44: (44.01299, 0.0),
  45: (45.01922, 0.0),
  46: (46.02738, 0.0)},
 'Pa': {0: (231.035884, 1.0),
  212: (212.0232, 0.0),
  213: (213.02111, 0.0),
  214: (214.02092, 0.0),
  215: (215.01919, 0.0),
  216: (216.01911, 0.0),
  217: (217.01832, 0.0),
  218: (218.020042, 0.0),
  219: (219.01988, 0.0),
  220: (220.02188, 0.0),
  221: (221.02188, 0.0),
  222: (222.02374, 0.0),
  223: (223.02396, 0.0),
  224: (224.025626, 0.0),
  225: (225.02613, 0.0),
  226: (226.027948, 0.0),
  227: (227.028805, 0.0),
  228: (228.031051, 0.0),
  229: (229.0320968, 0.0),
  230: (230.034541, 0.0),
  231: (231.035884, 1.0),
  232: (232.038592, 0.0),
  233: (233.0402473, 0.0),
  234: (234.043308, 0.0),
  235: (235.04544, 0.0),
  236: (236.04868, 0.0),
  237: (237.05115, 0.0),
  238: (238.0545, 0.0),
  239: (239.05726, 0.0),
  240: (240.06098, 0.0)},
 'Pb': {0: (207.9766521, 1.0),
  178: (178.00383, 0.0),
  179: (179.00215, 0.0),
  180: (179.997918, 0.0),
  181: (180.99662, 0.0),
  182: (181.992672, 0.0),
  183: (182.99187, 0.0),
  184: (183.988142, 0.0),
  185: (184.98761, 0.0),
  186: (185.984239, 0.0),
  187: (186.983918, 0.0),
  188: (187.980874, 0.0),
  189: (188.98081, 0.0),
  190: (189.978082, 0.0),
  191: (190.97827, 0.0),
  192: (191.975785, 0.0),
  193: (192.97617, 0.0),
  194: (193.974012, 0.0),
  195: (194.974542, 0.0),
  196: (195.972774, 0.0),
  197: (196.973431, 0.0),
  198: (197.972034, 0.0),
  199: (198.972917, 0.0),
  200: (199.971827, 0.0),
  201: (200.972885, 0.0),
  202: (201.972159, 0.0),
  203: (202.973391, 0.0),
  204: (203.9730436, 0.014),
  205: (204.9744818, 0.0),
  206: (205.9744653, 0.241),
  207: (206.9758969, 0.221),
  208: (207.9766521, 0.524),
  209: (208.9810901, 0.0),
  210: (209.9841885, 0.0),
  211: (210.988737, 0.0),
  212: (211.9918975, 0.0),
  213: (212.996581, 0.0),
  214: (213.9998054, 0.0),
  215: (215.00481, 0.0)},
 'Pd': {0: (105.903486, 1.0),
  91: (90.94911, 0.0),
  92: (91.94042, 0.0),
  93: (92.93591, 0.0),
  94: (93.92877, 0.0),
  95: (94.92469, 0.0),
  96: (95.91816, 0.0),
  97: (96.91648, 0.0),
  98: (97.912721, 0.0),
  99: (98.911768, 0.0),
  100: (99.908506, 0.0),
  101: (100.908289, 0.0),
  102: (101.905609, 0.0102),
  103: (102.906087, 0.0),
  104: (103.904036, 0.1114),
  105: (104.905085, 0.2233),
  106: (105.903486, 0.2733),
  107: (106.905133, 0.0),
  108: (107.903892, 0.2646),
  109: (108.90595, 0.0),
  110: (109.905153, 0.1172),
  111: (110.907671, 0.0),
  112: (111.907314, 0.0),
  113: (112.91015, 0.0),
  114: (113.910363, 0.0),
  115: (114.91368, 0.0),
  116: (115.91416, 0.0),
  117: (116.91784, 0.0),
  118: (117.91898, 0.0),
  119: (118.92311, 0.0),
  120: (119.92469, 0.0),
  121: (120.92887, 0.0),
  122: (121.93055, 0.0),
  123: (122.93493, 0.0),
  124: (123.93688, 0.0)},
 'Pm': {0: (145, 1.0),
  126: (125.95752, 0.0),
  127: (126.95163, 0.0),
  128: (127.94842, 0.0),
  129: (128.94316, 0.0),
  130: (129.94045, 0.0),
  131: (130.93587, 0.0),
  132: (131.93375, 0.0),
  133: (132.92978, 0.0),
  134: (133.92835, 0.0),
  135: (134.92488, 0.0),
  136: (135.92357, 0.0),
  137: (136.920479, 0.0),
  138: (137.919548, 0.0),
  139: (138.916804, 0.0),
  140: (139.91604, 0.0),
  141: (140.913555, 0.0),
  142: (141.912874, 0.0),
  143: (142.910933, 0.0),
  144: (143.912591, 0.0),
  145: (144.912749, 0.0),
  146: (145.914696, 0.0),
  147: (146.9151385, 0.0),
  148: (147.917475, 0.0),
  149: (148.918334, 0.0),
  150: (149.920984, 0.0),
  151: (150.921207, 0.0),
  152: (151.923497, 0.0),
  153: (152.924117, 0.0),
  154: (153.92646, 0.0),
  155: (154.9281, 0.0),
  156: (155.93106, 0.0),
  157: (156.93304, 0.0),
  158: (157.93656, 0.0),
  159: (158.93897, 0.0),
  160: (159.94299, 0.0),
  161: (160.94586, 0.0),
  162: (161.95029, 0.0),
  163: (162.95368, 0.0)},
 'Po': {0: (209, 1.0),
  188: (187.999422, 0.0),
  189: (188.998481, 0.0),
  190: (189.995101, 0.0),
  191: (190.994574, 0.0),
  192: (191.991335, 0.0),
  193: (192.99103, 0.0),
  194: (193.988186, 0.0),
  195: (194.98811, 0.0),
  196: (195.985535, 0.0),
  197: (196.98566, 0.0),
  198: (197.983389, 0.0),
  199: (198.983666, 0.0),
  200: (199.981799, 0.0),
  201: (200.98226, 0.0),
  202: (201.980758, 0.0),
  203: (202.98142, 0.0),
  204: (203.980318, 0.0),
  205: (204.981203, 0.0),
  206: (205.980481, 0.0),
  207: (206.981593, 0.0),
  208: (207.9812457, 0.0),
  209: (208.9824304, 0.0),
  210: (209.9828737, 0.0),
  211: (210.9866532, 0.0),
  212: (211.988868, 0.0),
  213: (212.992857, 0.0),
  214: (213.9952014, 0.0),
  215: (214.99942, 0.0),
  216: (216.001915, 0.0),
  217: (217.006335, 0.0),
  218: (218.008973, 0.0),
  219: (219.01374, 0.0),
  220: (220.0166, 0.0)},
 'Pr': {0: (140.9076528, 1.0),
  121: (120.95536, 0.0),
  122: (121.95181, 0.0),
  123: (122.94596, 0.0),
  124: (123.94296, 0.0),
  125: (124.93783, 0.0),
  126: (125.93531, 0.0),
  127: (126.93083, 0.0),
  128: (127.92879, 0.0),
  129: (128.9251, 0.0),
  130: (129.92359, 0.0),
  131: (130.92026, 0.0),
  132: (131.91926, 0.0),
  133: (132.916331, 0.0),
  134: (133.91571, 0.0),
  135: (134.913112, 0.0),
  136: (135.912692, 0.0),
  137: (136.910705, 0.0),
  138: (137.910755, 0.0),
  139: (138.908938, 0.0),
  140: (139.909076, 0.0),
  141: (140.9076528, 1.0),
  142: (141.9100448, 0.0),
  143: (142.9108169, 0.0),
  144: (143.913305, 0.0),
  145: (144.914512, 0.0),
  146: (145.91764, 0.0),
  147: (146.918996, 0.0),
  148: (147.922135, 0.0),
  149: (148.92372, 0.0),
  150: (149.926673, 0.0),
  151: (150.928319, 0.0),
  152: (151.9315, 0.0),
  153: (152.93384, 0.0),
  154: (153.93752, 0.0),
  155: (154.94012, 0.0),
  156: (155.94427, 0.0),
  157: (156.94743, 0.0),
  158: (157.95198, 0.0),
  159: (158.9555, 0.0)},
 'Pt': {0: (194.9647911, 1.0),
  166: (165.99486, 0.0),
  167: (166.99298, 0.0),
  168: (167.98815, 0.0),
  169: (168.98672, 0.0),
  170: (169.982495, 0.0),
  171: (170.98124, 0.0),
  172: (171.977347, 0.0),
  173: (172.97644, 0.0),
  174: (173.972819, 0.0),
  175: (174.972421, 0.0),
  176: (175.968945, 0.0),
  177: (176.968469, 0.0),
  178: (177.965649, 0.0),
  179: (178.965363, 0.0),
  180: (179.963031, 0.0),
  181: (180.963097, 0.0),
  182: (181.961171, 0.0),
  183: (182.961597, 0.0),
  184: (183.959922, 0.0),
  185: (184.96062, 0.0),
  186: (185.959351, 0.0),
  187: (186.96059, 0.0),
  188: (187.959395, 0.0),
  189: (188.960834, 0.0),
  190: (189.959932, 0.00014),
  191: (190.961677, 0.0),
  192: (191.961038, 0.00782),
  193: (192.9629874, 0.0),
  194: (193.9626803, 0.32967),
  195: (194.9647911, 0.33832),
  196: (195.9649515, 0.25242),
  197: (196.9673402, 0.0),
  198: (197.967893, 0.07163),
  199: (198.970593, 0.0),
  200: (199.971441, 0.0),
  201: (200.97451, 0.0),
  202: (201.97574, 0.0)},
 'Pu': {0: (244, 1.0),
  228: (228.03874, 0.0),
  229: (229.04015, 0.0),
  230: (230.03965, 0.0),
  231: (231.041101, 0.0),
  232: (232.041187, 0.0),
  233: (233.043, 0.0),
  234: (234.043317, 0.0),
  235: (235.045286, 0.0),
  236: (236.046058, 0.0),
  237: (237.0484097, 0.0),
  238: (238.0495599, 0.0),
  239: (239.0521634, 0.0),
  240: (240.0538135, 0.0),
  241: (241.0568515, 0.0),
  242: (242.0587426, 0.0),
  243: (243.062003, 0.0),
  244: (244.064204, 0.0),
  245: (245.067747, 0.0),
  246: (246.070205, 0.0),
  247: (247.07407, 0.0)},
 'Ra': {0: (226, 1.0),
  202: (202.00989, 0.0),
  203: (203.00927, 0.0),
  204: (204.0065, 0.0),
  205: (205.00627, 0.0),
  206: (206.003827, 0.0),
  207: (207.0038, 0.0),
  208: (208.00184, 0.0),
  209: (209.00199, 0.0),
  210: (210.000495, 0.0),
  211: (211.000898, 0.0),
  212: (211.999794, 0.0),
  213: (213.000384, 0.0),
  214: (214.000108, 0.0),
  215: (215.00272, 0.0),
  216: (216.003533, 0.0),
  217: (217.00632, 0.0),
  218: (218.00714, 0.0),
  219: (219.010085, 0.0),
  220: (220.011028, 0.0),
  221: (221.013917, 0.0),
  222: (222.015375, 0.0),
  223: (223.0185022, 0.0),
  224: (224.0202118, 0.0),
  225: (225.023612, 0.0),
  226: (226.0254098, 0.0),
  227: (227.0291778, 0.0),
  228: (228.0310703, 0.0),
  229: (229.034958, 0.0),
  230: (230.037056, 0.0),
  231: (231.04122, 0.0),
  232: (232.04364, 0.0),
  233: (233.04806, 0.0),
  234: (234.0507, 0.0)},
 'Rb': {0: (84.911789738, 1.0),
  71: (70.96532, 0.0),
  72: (71.95908, 0.0),
  73: (72.95056, 0.0),
  74: (73.944265, 0.0),
  75: (74.93857, 0.0),
  76: (75.9350722, 0.0),
  77: (76.930408, 0.0),
  78: (77.928141, 0.0),
  79: (78.923989, 0.0),
  80: (79.922519, 0.0),
  81: (80.918996, 0.0),
  82: (81.9182086, 0.0),
  83: (82.91511, 0.0),
  84: (83.914385, 0.0),
  85: (84.911789738, 0.7217),
  86: (85.91116742, 0.0),
  87: (86.909180527, 0.2783),
  88: (87.91131559, 0.0),
  89: (88.912278, 0.0),
  90: (89.914802, 0.0),
  91: (90.916537, 0.0),
  92: (91.919729, 0.0),
  93: (92.922042, 0.0),
  94: (93.926405, 0.0),
  95: (94.929303, 0.0),
  96: (95.93427, 0.0),
  97: (96.93735, 0.0),
  98: (97.94179, 0.0),
  99: (98.94538, 0.0),
  100: (99.94987, 0.0),
  101: (100.9532, 0.0),
  102: (101.95887, 0.0)},
 'Re': {0: (186.9557531, 1.0),
  160: (159.98212, 0.0),
  161: (160.97759, 0.0),
  162: (161.976, 0.0),
  163: (162.972081, 0.0),
  164: (163.97032, 0.0),
  165: (164.967089, 0.0),
  166: (165.96581, 0.0),
  167: (166.9626, 0.0),
  168: (167.96157, 0.0),
  169: (168.95879, 0.0),
  170: (169.95822, 0.0),
  171: (170.95572, 0.0),
  172: (171.95542, 0.0),
  173: (172.95324, 0.0),
  174: (173.95312, 0.0),
  175: (174.95138, 0.0),
  176: (175.95162, 0.0),
  177: (176.95033, 0.0),
  178: (177.95099, 0.0),
  179: (178.949988, 0.0),
  180: (179.950789, 0.0),
  181: (180.950068, 0.0),
  182: (181.95121, 0.0),
  183: (182.95082, 0.0),
  184: (183.952521, 0.0),
  185: (184.952955, 0.374),
  186: (185.9549861, 0.0),
  187: (186.9557531, 0.626),
  188: (187.9581144, 0.0),
  189: (188.959229, 0.0),
  190: (189.96182, 0.0),
  191: (190.963125, 0.0),
  192: (191.96596, 0.0),
  193: (192.96747, 0.0),
  194: (193.97042, 0.0)},
 'Rf': {0: (265, 1.0),
  253: (253.10069, 0.0),
  254: (254.10018, 0.0),
  255: (255.10134, 0.0),
  256: (256.101166, 0.0),
  257: (257.10299, 0.0),
  258: (258.10349, 0.0),
  259: (259.10564, 0.0),
  260: (260.10644, 0.0),
  261: (261.10877, 0.0),
  262: (262.10993, 0.0),
  263: (263.11255, 0.0),
  264: (264.11399, 0.0),
  265: (265.1167, 0.0),
  266: (266.11796, 0.0),
  267: (267.12153, 0.0),
  268: (268.12364, 0.0)},
 'Rg': {0: (280, 1.0),
  272: (272.15362, 0.0),
  273: (273.15368, 0.0),
  274: (274.15571, 0.0),
  275: (275.15614, 0.0),
  276: (276.15849, 0.0),
  277: (277.15952, 0.0),
  278: (278.1616, 0.0),
  279: (279.16247, 0.0),
  280: (280.16447, 0.0),
  281: (281.16537, 0.0),
  282: (282.16749, 0.0),
  283: (283.16842, 0.0)},
 'Rh': {0: (102.905504, 1.0),
  89: (88.94884, 0.0),
  90: (89.94287, 0.0),
  91: (90.93655, 0.0),
  92: (91.93198, 0.0),
  93: (92.92574, 0.0),
  94: (93.9217, 0.0),
  95: (94.9159, 0.0),
  96: (95.914461, 0.0),
  97: (96.91134, 0.0),
  98: (97.910708, 0.0),
  99: (98.908132, 0.0),
  100: (99.908122, 0.0),
  101: (100.906164, 0.0),
  102: (101.906843, 0.0),
  103: (102.905504, 1.0),
  104: (103.906656, 0.0),
  105: (104.905694, 0.0),
  106: (105.907287, 0.0),
  107: (106.906748, 0.0),
  108: (107.90873, 0.0),
  109: (108.908737, 0.0),
  110: (109.91114, 0.0),
  111: (110.91159, 0.0),
  112: (111.91439, 0.0),
  113: (112.91553, 0.0),
  114: (113.91881, 0.0),
  115: (114.92033, 0.0),
  116: (115.92406, 0.0),
  117: (116.92598, 0.0),
  118: (117.93007, 0.0),
  119: (118.93211, 0.0),
  120: (119.93641, 0.0),
  121: (120.93872, 0.0),
  122: (121.94321, 0.0)},
 'Rn': {0: (222, 1.0),
  195: (195.00544, 0.0),
  196: (196.002115, 0.0),
  197: (197.00158, 0.0),
  198: (197.998679, 0.0),
  199: (198.99837, 0.0),
  200: (199.995699, 0.0),
  201: (200.99563, 0.0),
  202: (201.993263, 0.0),
  203: (202.993387, 0.0),
  204: (203.991429, 0.0),
  205: (204.99172, 0.0),
  206: (205.990214, 0.0),
  207: (206.990734, 0.0),
  208: (207.989642, 0.0),
  209: (208.990415, 0.0),
  210: (209.989696, 0.0),
  211: (210.990601, 0.0),
  212: (211.990704, 0.0),
  213: (212.993883, 0.0),
  214: (213.995363, 0.0),
  215: (214.998745, 0.0),
  216: (216.000274, 0.0),
  217: (217.003928, 0.0),
  218: (218.0056013, 0.0),
  219: (219.0094802, 0.0),
  220: (220.011394, 0.0),
  221: (221.015537, 0.0),
  222: (222.0175777, 0.0),
  223: (223.02179, 0.0),
  224: (224.02409, 0.0),
  225: (225.02844, 0.0),
  226: (226.03089, 0.0),
  227: (227.03541, 0.0),
  228: (228.03799, 0.0)},
 'Ru': {0: (101.9043493, 1.0),
  87: (86.94918, 0.0),
  88: (87.94026, 0.0),
  89: (88.93611, 0.0),
  90: (89.92989, 0.0),
  91: (90.92629, 0.0),
  92: (91.92012, 0.0),
  93: (92.91705, 0.0),
  94: (93.91136, 0.0),
  95: (94.910413, 0.0),
  96: (95.907598, 0.0554),
  97: (96.907555, 0.0),
  98: (97.905287, 0.0187),
  99: (98.9059393, 0.1276),
  100: (99.9042195, 0.126),
  101: (100.9055821, 0.1706),
  102: (101.9043493, 0.3155),
  103: (102.9063238, 0.0),
  104: (103.905433, 0.1862),
  105: (104.907753, 0.0),
  106: (105.907329, 0.0),
  107: (106.90991, 0.0),
  108: (107.91017, 0.0),
  109: (108.9132, 0.0),
  110: (109.91414, 0.0),
  111: (110.9177, 0.0),
  112: (111.91897, 0.0),
  113: (112.92249, 0.0),
  114: (113.92428, 0.0),
  115: (114.92869, 0.0),
  116: (115.93081, 0.0),
  117: (116.93558, 0.0),
  118: (117.93782, 0.0),
  119: (118.94284, 0.0),
  120: (119.94531, 0.0)},
 'S': {0: (31.972071, 1.0),
  26: (26.02788, 0.0),
  27: (27.01883, 0.0),
  28: (28.00437, 0.0),
  29: (28.99661, 0.0),
  30: (29.984903, 0.0),
  31: (30.9795547, 0.0),
  32: (31.972071, 0.9499),
  33: (32.97145876, 0.0075),
  34: (33.9678669, 0.0425),
  35: (34.96903216, 0.0),
  36: (35.96708076, 0.0001),
  37: (36.97112557, 0.0),
  38: (37.971163, 0.0),
  39: (38.97513, 0.0),
  40: (39.97545, 0.0),
  41: (40.97958, 0.0),
  42: (41.98102, 0.0),
  43: (42.98715, 0.0),
  44: (43.99021, 0.0),
  45: (44.99651, 0.0),
  46: (46.00075, 0.0),
  47: (47.00859, 0.0),
  48: (48.01417, 0.0),
  49: (49.02362, 0.0)},
 'Sb': {0: (120.9038157, 1.0),
  103: (102.93969, 0.0),
  104: (103.93647, 0.0),
  105: (104.93149, 0.0),
  106: (105.92879, 0.0),
  107: (106.92415, 0.0),
  108: (107.92216, 0.0),
  109: (108.918132, 0.0),
  110: (109.91675, 0.0),
  111: (110.91316, 0.0),
  112: (111.912398, 0.0),
  113: (112.909372, 0.0),
  114: (113.90927, 0.0),
  115: (114.906598, 0.0),
  116: (115.906794, 0.0),
  117: (116.904836, 0.0),
  118: (117.905529, 0.0),
  119: (118.903942, 0.0),
  120: (119.905072, 0.0),
  121: (120.9038157, 0.5721),
  122: (121.9051737, 0.0),
  123: (122.904214, 0.4279),
  124: (123.9059357, 0.0),
  125: (124.9052538, 0.0),
  126: (125.90725, 0.0),
  127: (126.906924, 0.0),
  128: (127.909169, 0.0),
  129: (128.909148, 0.0),
  130: (129.911656, 0.0),
  131: (130.911982, 0.0),
  132: (131.914467, 0.0),
  133: (132.915252, 0.0),
  134: (133.92038, 0.0),
  135: (134.92517, 0.0),
  136: (135.93035, 0.0),
  137: (136.93531, 0.0),
  138: (137.94079, 0.0),
  139: (138.94598, 0.0)},
 'Sc': {0: (44.9559119, 1.0),
  36: (36.01492, 0.0),
  37: (37.00305, 0.0),
  38: (37.9947, 0.0),
  39: (38.98479, 0.0),
  40: (39.977967, 0.0),
  41: (40.96925113, 0.0),
  42: (41.96551643, 0.0),
  43: (42.9611507, 0.0),
  44: (43.9594028, 0.0),
  45: (44.9559119, 1.0),
  46: (45.9551719, 0.0),
  47: (46.9524075, 0.0),
  48: (47.952231, 0.0),
  49: (48.950024, 0.0),
  50: (49.952188, 0.0),
  51: (50.953603, 0.0),
  52: (51.95668, 0.0),
  53: (52.95961, 0.0),
  54: (53.96326, 0.0),
  55: (54.96824, 0.0),
  56: (55.97287, 0.0),
  57: (56.97779, 0.0),
  58: (57.98371, 0.0),
  59: (58.98922, 0.0),
  60: (59.99571, 0.0)},
 'Se': {0: (79.9165213, 1.0),
  65: (64.96466, 0.0),
  66: (65.95521, 0.0),
  67: (66.95009, 0.0),
  68: (67.9418, 0.0),
  69: (68.93956, 0.0),
  70: (69.93339, 0.0),
  71: (70.93224, 0.0),
  72: (71.927112, 0.0),
  73: (72.926765, 0.0),
  74: (73.9224764, 0.0089),
  75: (74.9225234, 0.0),
  76: (75.9192136, 0.0937),
  77: (76.919914, 0.0763),
  78: (77.9173091, 0.2377),
  79: (78.9184991, 0.0),
  80: (79.9165213, 0.4961),
  81: (80.9179925, 0.0),
  82: (81.9166994, 0.0873),
  83: (82.919118, 0.0),
  84: (83.918462, 0.0),
  85: (84.92225, 0.0),
  86: (85.924272, 0.0),
  87: (86.92852, 0.0),
  88: (87.93142, 0.0),
  89: (88.93645, 0.0),
  90: (89.93996, 0.0),
  91: (90.94596, 0.0),
  92: (91.94992, 0.0),
  93: (92.95629, 0.0),
  94: (93.96049, 0.0)},
 'Sg': {0: (271, 1.0),
  258: (258.11317, 0.0),
  259: (259.1145, 0.0),
  260: (260.11442, 0.0),
  261: (261.11612, 0.0),
  262: (262.1164, 0.0),
  263: (263.11832, 0.0),
  264: (264.11893, 0.0),
  265: (265.12111, 0.0),
  266: (266.12207, 0.0),
  267: (267.12443, 0.0),
  268: (268.12561, 0.0),
  269: (269.12876, 0.0),
  270: (270.13033, 0.0),
  271: (271.13347, 0.0),
  272: (272.13516, 0.0),
  273: (273.13822, 0.0)},
 'Si': {0: (27.9769265325, 1.0),
  22: (22.03453, 0.0),
  23: (23.02552, 0.0),
  24: (24.011546, 0.0),
  25: (25.004106, 0.0),
  26: (25.99233, 0.0),
  27: (26.98670491, 0.0),
  28: (27.9769265325, 0.92223),
  29: (28.9764947, 0.04685),
  30: (29.97377017, 0.03092),
  31: (30.97536323, 0.0),
  32: (31.97414808, 0.0),
  33: (32.978, 0.0),
  34: (33.978576, 0.0),
  35: (34.98458, 0.0),
  36: (35.9866, 0.0),
  37: (36.99294, 0.0),
  38: (37.99563, 0.0),
  39: (39.00207, 0.0),
  40: (40.00587, 0.0),
  41: (41.01456, 0.0),
  42: (42.01979, 0.0),
  43: (43.02866, 0.0),
  44: (44.03526, 0.0)},
 'Sm': {0: (151.9197324, 1.0),
  128: (127.95808, 0.0),
  129: (128.95464, 0.0),
  130: (129.94892, 0.0),
  131: (130.94611, 0.0),
  132: (131.94069, 0.0),
  133: (132.93867, 0.0),
  134: (133.93397, 0.0),
  135: (134.93252, 0.0),
  136: (135.928276, 0.0),
  137: (136.92697, 0.0),
  138: (137.923244, 0.0),
  139: (138.922297, 0.0),
  140: (139.918995, 0.0),
  141: (140.918476, 0.0),
  142: (141.915198, 0.0),
  143: (142.914628, 0.0),
  144: (143.911999, 0.0307),
  145: (144.91341, 0.0),
  146: (145.913041, 0.0),
  147: (146.9148979, 0.1499),
  148: (147.9148227, 0.1124),
  149: (148.9171847, 0.1382),
  150: (149.9172755, 0.0738),
  151: (150.9199324, 0.0),
  152: (151.9197324, 0.2675),
  153: (152.9220974, 0.0),
  154: (153.9222093, 0.2275),
  155: (154.9246402, 0.0),
  156: (155.925528, 0.0),
  157: (156.92836, 0.0),
  158: (157.92999, 0.0),
  159: (158.93321, 0.0),
  160: (159.93514, 0.0),
  161: (160.93883, 0.0),
  162: (161.94122, 0.0),
  163: (162.94536, 0.0),
  164: (163.94828, 0.0),
  165: (164.95298, 0.0)},
 'Sn': {0: (119.9021947, 1.0),
  99: (98.94933, 0.0),
  100: (99.93904, 0.0),
  101: (100.93606, 0.0),
  102: (101.9303, 0.0),
  103: (102.9281, 0.0),
  104: (103.92314, 0.0),
  105: (104.92135, 0.0),
  106: (105.91688, 0.0),
  107: (106.91564, 0.0),
  108: (107.911925, 0.0),
  109: (108.911283, 0.0),
  110: (109.907843, 0.0),
  111: (110.907734, 0.0),
  112: (111.904818, 0.0097),
  113: (112.905171, 0.0),
  114: (113.902779, 0.0066),
  115: (114.903342, 0.0034),
  116: (115.901741, 0.1454),
  117: (116.902952, 0.0768),
  118: (117.901603, 0.2422),
  119: (118.903308, 0.0859),
  120: (119.9021947, 0.3258),
  121: (120.9042355, 0.0),
  122: (121.903439, 0.0463),
  123: (122.9057208, 0.0),
  124: (123.9052739, 0.0579),
  125: (124.9077841, 0.0),
  126: (125.907653, 0.0),
  127: (126.91036, 0.0),
  128: (127.910537, 0.0),
  129: (128.91348, 0.0),
  130: (129.913967, 0.0),
  131: (130.917, 0.0),
  132: (131.917816, 0.0),
  133: (132.92383, 0.0),
  134: (133.92829, 0.0),
  135: (134.93473, 0.0),
  136: (135.93934, 0.0),
  137: (136.94599, 0.0)},
 'Sr': {0: (87.9056121, 1.0),
  73: (72.96597, 0.0),
  74: (73.95631, 0.0),
  75: (74.94995, 0.0),
  76: (75.94177, 0.0),
  77: (76.937945, 0.0),
  78: (77.93218, 0.0),
  79: (78.929708, 0.0),
  80: (79.924521, 0.0),
  81: (80.923212, 0.0),
  82: (81.918402, 0.0),
  83: (82.917557, 0.0),
  84: (83.913425, 0.0056),
  85: (84.912933, 0.0),
  86: (85.9092602, 0.0986),
  87: (86.9088771, 0.07),
  88: (87.9056121, 0.8258),
  89: (88.9074507, 0.0),
  90: (89.907738, 0.0),
  91: (90.910203, 0.0),
  92: (91.911038, 0.0),
  93: (92.914026, 0.0),
  94: (93.915361, 0.0),
  95: (94.919359, 0.0),
  96: (95.921697, 0.0),
  97: (96.926153, 0.0),
  98: (97.928453, 0.0),
  99: (98.93324, 0.0),
  100: (99.93535, 0.0),
  101: (100.94052, 0.0),
  102: (101.94302, 0.0),
  103: (102.94895, 0.0),
  104: (103.95233, 0.0),
  105: (104.95858, 0.0)},
 'Ta': {0: (180.9479958, 1.0),
  155: (154.97459, 0.0),
  156: (155.9723, 0.0),
  157: (156.96819, 0.0),
  158: (157.9667, 0.0),
  159: (158.963018, 0.0),
  160: (159.96149, 0.0),
  161: (160.95842, 0.0),
  162: (161.95729, 0.0),
  163: (162.95433, 0.0),
  164: (163.95353, 0.0),
  165: (164.950773, 0.0),
  166: (165.95051, 0.0),
  167: (166.94809, 0.0),
  168: (167.94805, 0.0),
  169: (168.94601, 0.0),
  170: (169.94618, 0.0),
  171: (170.94448, 0.0),
  172: (171.9449, 0.0),
  173: (172.94375, 0.0),
  174: (173.94445, 0.0),
  175: (174.94374, 0.0),
  176: (175.94486, 0.0),
  177: (176.944472, 0.0),
  178: (177.945778, 0.0),
  179: (178.9459295, 0.0),
  180: (179.9474648, 0.00012),
  181: (180.9479958, 0.99988),
  182: (181.9501518, 0.0),
  183: (182.9513726, 0.0),
  184: (183.954008, 0.0),
  185: (184.955559, 0.0),
  186: (185.95855, 0.0),
  187: (186.96053, 0.0),
  188: (187.9637, 0.0),
  189: (188.96583, 0.0),
  190: (189.96923, 0.0)},
 'Tb': {0: (158.9253468, 1.0),
  136: (135.96138, 0.0),
  137: (136.95598, 0.0),
  138: (137.95316, 0.0),
  139: (138.94829, 0.0),
  140: (139.94581, 0.0),
  141: (140.94145, 0.0),
  142: (141.93874, 0.0),
  143: (142.93512, 0.0),
  144: (143.93305, 0.0),
  145: (144.92927, 0.0),
  146: (145.92725, 0.0),
  147: (146.924045, 0.0),
  148: (147.924272, 0.0),
  149: (148.923246, 0.0),
  150: (149.92366, 0.0),
  151: (150.923103, 0.0),
  152: (151.92407, 0.0),
  153: (152.923435, 0.0),
  154: (153.92468, 0.0),
  155: (154.923505, 0.0),
  156: (155.924747, 0.0),
  157: (156.9240246, 0.0),
  158: (157.9254131, 0.0),
  159: (158.9253468, 1.0),
  160: (159.9271676, 0.0),
  161: (160.9275699, 0.0),
  162: (161.92949, 0.0),
  163: (162.930648, 0.0),
  164: (163.93335, 0.0),
  165: (164.93488, 0.0),
  166: (165.93799, 0.0),
  167: (166.94005, 0.0),
  168: (167.94364, 0.0),
  169: (168.94622, 0.0),
  170: (169.95025, 0.0),
  171: (170.9533, 0.0)},
 'Tc': {0: (98, 1.0),
  85: (84.94883, 0.0),
  86: (85.94288, 0.0),
  87: (86.93653, 0.0),
  88: (87.93268, 0.0),
  89: (88.92717, 0.0),
  90: (89.92356, 0.0),
  91: (90.91843, 0.0),
  92: (91.91526, 0.0),
  93: (92.910249, 0.0),
  94: (93.909657, 0.0),
  95: (94.907657, 0.0),
  96: (95.907871, 0.0),
  97: (96.906365, 0.0),
  98: (97.907216, 0.0),
  99: (98.9062547, 0.0),
  100: (99.9076578, 0.0),
  101: (100.907315, 0.0),
  102: (101.909215, 0.0),
  103: (102.909181, 0.0),
  104: (103.91145, 0.0),
  105: (104.91166, 0.0),
  106: (105.914358, 0.0),
  107: (106.91508, 0.0),
  108: (107.91846, 0.0),
  109: (108.91998, 0.0),
  110: (109.92382, 0.0),
  111: (110.92569, 0.0),
  112: (111.92915, 0.0),
  113: (112.93159, 0.0),
  114: (113.93588, 0.0),
  115: (114.93869, 0.0),
  116: (115.94337, 0.0),
  117: (116.94648, 0.0),
  118: (117.95148, 0.0)},
 'Te': {0: (129.9062244, 1.0),
  105: (104.94364, 0.0),
  106: (105.9375, 0.0),
  107: (106.93501, 0.0),
  108: (107.92944, 0.0),
  109: (108.92742, 0.0),
  110: (109.92241, 0.0),
  111: (110.92111, 0.0),
  112: (111.91701, 0.0),
  113: (112.91589, 0.0),
  114: (113.91209, 0.0),
  115: (114.9119, 0.0),
  116: (115.90846, 0.0),
  117: (116.908645, 0.0),
  118: (117.905828, 0.0),
  119: (118.906404, 0.0),
  120: (119.90402, 0.0009),
  121: (120.904936, 0.0),
  122: (121.9030439, 0.0255),
  123: (122.90427, 0.0089),
  124: (123.9028179, 0.0474),
  125: (124.9044307, 0.0707),
  126: (125.9033117, 0.1884),
  127: (126.9052263, 0.0),
  128: (127.9044631, 0.3174),
  129: (128.9065982, 0.0),
  130: (129.9062244, 0.3408),
  131: (130.9085239, 0.0),
  132: (131.908553, 0.0),
  133: (132.910955, 0.0),
  134: (133.911369, 0.0),
  135: (134.91645, 0.0),
  136: (135.9201, 0.0),
  137: (136.92532, 0.0),
  138: (137.92922, 0.0),
  139: (138.93473, 0.0),
  140: (139.93885, 0.0),
  141: (140.94465, 0.0),
  142: (141.94908, 0.0)},
 'Th': {0: (232.0380553, 1.0),
  209: (209.01772, 0.0),
  210: (210.015075, 0.0),
  211: (211.01493, 0.0),
  212: (212.01298, 0.0),
  213: (213.01301, 0.0),
  214: (214.0115, 0.0),
  215: (215.01173, 0.0),
  216: (216.011062, 0.0),
  217: (217.013114, 0.0),
  218: (218.013284, 0.0),
  219: (219.01554, 0.0),
  220: (220.015748, 0.0),
  221: (221.018184, 0.0),
  222: (222.018468, 0.0),
  223: (223.020811, 0.0),
  224: (224.021467, 0.0),
  225: (225.023951, 0.0),
  226: (226.024903, 0.0),
  227: (227.0277041, 0.0),
  228: (228.0287411, 0.0),
  229: (229.031762, 0.0),
  230: (230.0331338, 0.0),
  231: (231.0363043, 0.0),
  232: (232.0380553, 1.0),
  233: (233.0415818, 0.0),
  234: (234.043601, 0.0),
  235: (235.04751, 0.0),
  236: (236.04987, 0.0),
  237: (237.05389, 0.0),
  238: (238.0565, 0.0)},
 'Ti': {0: (47.9479463, 1.0),
  38: (38.00977, 0.0),
  39: (39.00161, 0.0),
  40: (39.9905, 0.0),
  41: (40.98315, 0.0),
  42: (41.973031, 0.0),
  43: (42.968522, 0.0),
  44: (43.9596901, 0.0),
  45: (44.9581256, 0.0),
  46: (45.9526316, 0.0825),
  47: (46.9517631, 0.0744),
  48: (47.9479463, 0.7372),
  49: (48.94787, 0.0541),
  50: (49.9447912, 0.0518),
  51: (50.946615, 0.0),
  52: (51.946897, 0.0),
  53: (52.94973, 0.0),
  54: (53.95105, 0.0),
  55: (54.95527, 0.0),
  56: (55.9582, 0.0),
  57: (56.96399, 0.0),
  58: (57.96697, 0.0),
  59: (58.97293, 0.0),
  60: (59.97676, 0.0),
  61: (60.9832, 0.0),
  62: (61.98749, 0.0),
  63: (62.99442, 0.0)},
 'Tl': {0: (204.9744275, 1.0),
  176: (176.00059, 0.0),
  177: (176.996427, 0.0),
  178: (177.9949, 0.0),
  179: (178.99109, 0.0),
  180: (179.98991, 0.0),
  181: (180.986257, 0.0),
  182: (181.98567, 0.0),
  183: (182.982193, 0.0),
  184: (183.98187, 0.0),
  185: (184.97879, 0.0),
  186: (185.97833, 0.0),
  187: (186.975906, 0.0),
  188: (187.97601, 0.0),
  189: (188.973588, 0.0),
  190: (189.97388, 0.0),
  191: (190.971786, 0.0),
  192: (191.97223, 0.0),
  193: (192.97067, 0.0),
  194: (193.9712, 0.0),
  195: (194.969774, 0.0),
  196: (195.970481, 0.0),
  197: (196.969575, 0.0),
  198: (197.97048, 0.0),
  199: (198.96988, 0.0),
  200: (199.970963, 0.0),
  201: (200.970819, 0.0),
  202: (201.972106, 0.0),
  203: (202.9723442, 0.2952),
  204: (203.9738635, 0.0),
  205: (204.9744275, 0.7048),
  206: (205.9761103, 0.0),
  207: (206.977419, 0.0),
  208: (207.9820187, 0.0),
  209: (208.985359, 0.0),
  210: (209.990074, 0.0),
  211: (210.99348, 0.0),
  212: (211.99823, 0.0)},
 'Tm': {0: (168.9342133, 1.0),
  145: (144.97007, 0.0),
  146: (145.96643, 0.0),
  147: (146.96096, 0.0),
  148: (147.95784, 0.0),
  149: (148.95272, 0.0),
  150: (149.94996, 0.0),
  151: (150.945483, 0.0),
  152: (151.94442, 0.0),
  153: (152.942012, 0.0),
  154: (153.941568, 0.0),
  155: (154.939199, 0.0),
  156: (155.93898, 0.0),
  157: (156.93697, 0.0),
  158: (157.93698, 0.0),
  159: (158.93498, 0.0),
  160: (159.93526, 0.0),
  161: (160.93355, 0.0),
  162: (161.933995, 0.0),
  163: (162.932651, 0.0),
  164: (163.93356, 0.0),
  165: (164.932435, 0.0),
  166: (165.933554, 0.0),
  167: (166.9328516, 0.0),
  168: (167.934173, 0.0),
  169: (168.9342133, 1.0),
  170: (169.9358014, 0.0),
  171: (170.9364294, 0.0),
  172: (171.9384, 0.0),
  173: (172.939604, 0.0),
  174: (173.94217, 0.0),
  175: (174.94384, 0.0),
  176: (175.94699, 0.0),
  177: (176.94904, 0.0),
  178: (177.95264, 0.0),
  179: (178.95534, 0.0)},
 'U': {0: (238.0507882, 1.0),
  217: (217.02437, 0.0),
  218: (218.02354, 0.0),
  219: (219.02492, 0.0),
  220: (220.02472, 0.0),
  221: (221.0264, 0.0),
  222: (222.02609, 0.0),
  223: (223.02774, 0.0),
  224: (224.027605, 0.0),
  225: (225.029391, 0.0),
  226: (226.029339, 0.0),
  227: (227.031156, 0.0),
  228: (228.031374, 0.0),
  229: (229.033506, 0.0),
  230: (230.03394, 0.0),
  231: (231.036294, 0.0),
  232: (232.0371562, 0.0),
  233: (233.0396352, 0.0),
  234: (234.0409521, 5.4e-05),
  235: (235.0439299, 0.007204),
  236: (236.045568, 0.0),
  237: (237.0487302, 0.0),
  238: (238.0507882, 0.992742),
  239: (239.0542933, 0.0),
  240: (240.056592, 0.0),
  241: (241.06033, 0.0),
  242: (242.06293, 0.0)},
 'Uuh': {0: (293, 1.0),
  289: (289.19886, 0.0),
  290: (290.19859, 0.0),
  291: (291.20001, 0.0),
  292: (292.19979, 0.0)},
 'Uuo': {0: (294, 1.0), 293: (293.21467, 0.0)},
 'Uup': {0: (288, 1.0),
  287: (287.19119, 0.0),
  288: (288.19249, 0.0),
  289: (289.19272, 0.0),
  290: (290.19414, 0.0),
  291: (291.19438, 0.0)},
 'Uuq': {0: (289, 1.0),
  285: (285.1837, 0.0),
  286: (286.18386, 0.0),
  287: (287.1856, 0.0),
  288: (288.18569, 0.0),
  289: (289.18728, 0.0)},
 'Uus': {0: (292, 1.0), 291: (291.20656, 0.0), 292: (292.20755, 0.0)},
 'Uut': {0: (284, 1.0),
  283: (283.17645, 0.0),
  284: (284.17808, 0.0),
  285: (285.17873, 0.0),
  286: (286.18048, 0.0),
  287: (287.18105, 0.0)},
 'V': {0: (50.9439595, 1.0),
  40: (40.01109, 0.0),
  41: (40.99978, 0.0),
  42: (41.99123, 0.0),
  43: (42.98065, 0.0),
  44: (43.97411, 0.0),
  45: (44.965776, 0.0),
  46: (45.9602005, 0.0),
  47: (46.9549089, 0.0),
  48: (47.9522537, 0.0),
  49: (48.9485161, 0.0),
  50: (49.9471585, 0.0025),
  51: (50.9439595, 0.9975),
  52: (51.9447755, 0.0),
  53: (52.944338, 0.0),
  54: (53.94644, 0.0),
  55: (54.94723, 0.0),
  56: (55.95053, 0.0),
  57: (56.95256, 0.0),
  58: (57.95683, 0.0),
  59: (58.96021, 0.0),
  60: (59.96503, 0.0),
  61: (60.96848, 0.0),
  62: (61.97378, 0.0),
  63: (62.97755, 0.0),
  64: (63.98347, 0.0),
  65: (64.98792, 0.0)},
 'W': {0: (183.9509312, 1.0),
  158: (157.97456, 0.0),
  159: (158.97292, 0.0),
  160: (159.96848, 0.0),
  161: (160.96736, 0.0),
  162: (161.963497, 0.0),
  163: (162.96252, 0.0),
  164: (163.958954, 0.0),
  165: (164.95828, 0.0),
  166: (165.955027, 0.0),
  167: (166.954816, 0.0),
  168: (167.951808, 0.0),
  169: (168.951779, 0.0),
  170: (169.949228, 0.0),
  171: (170.94945, 0.0),
  172: (171.94729, 0.0),
  173: (172.94769, 0.0),
  174: (173.94608, 0.0),
  175: (174.94672, 0.0),
  176: (175.94563, 0.0),
  177: (176.94664, 0.0),
  178: (177.945876, 0.0),
  179: (178.94707, 0.0),
  180: (179.946704, 0.0012),
  181: (180.948197, 0.0),
  182: (181.9482042, 0.265),
  183: (182.950223, 0.1431),
  184: (183.9509312, 0.3064),
  185: (184.9534193, 0.0),
  186: (185.9543641, 0.2843),
  187: (186.9571605, 0.0),
  188: (187.958489, 0.0),
  189: (188.96191, 0.0),
  190: (189.96318, 0.0),
  191: (190.9666, 0.0),
  192: (191.96817, 0.0)},
 'Xe': {0: (131.9041535, 1.0),
  110: (109.94428, 0.0),
  111: (110.9416, 0.0),
  112: (111.93562, 0.0),
  113: (112.93334, 0.0),
  114: (113.92798, 0.0),
  115: (114.926294, 0.0),
  116: (115.921581, 0.0),
  117: (116.920359, 0.0),
  118: (117.916179, 0.0),
  119: (118.915411, 0.0),
  120: (119.911784, 0.0),
  121: (120.911462, 0.0),
  122: (121.908368, 0.0),
  123: (122.908482, 0.0),
  124: (123.905893, 0.000952),
  125: (124.9063955, 0.0),
  126: (125.904274, 0.00089),
  127: (126.905184, 0.0),
  128: (127.9035313, 0.019102),
  129: (128.9047794, 0.264006),
  130: (129.903508, 0.04071),
  131: (130.9050824, 0.212324),
  132: (131.9041535, 0.269086),
  133: (132.9059107, 0.0),
  134: (133.9053945, 0.104357),
  135: (134.907227, 0.0),
  136: (135.907219, 0.088573),
  137: (136.911562, 0.0),
  138: (137.91395, 0.0),
  139: (138.918793, 0.0),
  140: (139.92164, 0.0),
  141: (140.92665, 0.0),
  142: (141.92971, 0.0),
  143: (142.93511, 0.0),
  144: (143.93851, 0.0),
  145: (144.94407, 0.0),
  146: (145.94775, 0.0),
  147: (146.95356, 0.0)},
 'Y': {0: (88.9058483, 1.0),
  76: (75.95845, 0.0),
  77: (76.94965, 0.0),
  78: (77.94361, 0.0),
  79: (78.93735, 0.0),
  80: (79.93428, 0.0),
  81: (80.92913, 0.0),
  82: (81.92679, 0.0),
  83: (82.92235, 0.0),
  84: (83.92039, 0.0),
  85: (84.916433, 0.0),
  86: (85.914886, 0.0),
  87: (86.9108757, 0.0),
  88: (87.9095011, 0.0),
  89: (88.9058483, 1.0),
  90: (89.9071519, 0.0),
  91: (90.907305, 0.0),
  92: (91.908949, 0.0),
  93: (92.909583, 0.0),
  94: (93.911595, 0.0),
  95: (94.912821, 0.0),
  96: (95.915891, 0.0),
  97: (96.918134, 0.0),
  98: (97.922203, 0.0),
  99: (98.924636, 0.0),
  100: (99.92776, 0.0),
  101: (100.93031, 0.0),
  102: (101.93356, 0.0),
  103: (102.93673, 0.0),
  104: (103.94105, 0.0),
  105: (104.94487, 0.0),
  106: (105.94979, 0.0),
  107: (106.95414, 0.0),
  108: (107.95948, 0.0)},
 'Yb': {0: (173.9388621, 1.0),
  148: (147.96742, 0.0),
  149: (148.96404, 0.0),
  150: (149.95842, 0.0),
  151: (150.9554, 0.0),
  152: (151.95029, 0.0),
  153: (152.94948, 0.0),
  154: (153.946394, 0.0),
  155: (154.945782, 0.0),
  156: (155.942818, 0.0),
  157: (156.942628, 0.0),
  158: (157.939866, 0.0),
  159: (158.94005, 0.0),
  160: (159.937552, 0.0),
  161: (160.937902, 0.0),
  162: (161.935768, 0.0),
  163: (162.936334, 0.0),
  164: (163.934489, 0.0),
  165: (164.93528, 0.0),
  166: (165.933882, 0.0),
  167: (166.93495, 0.0),
  168: (167.933897, 0.0013),
  169: (168.93519, 0.0),
  170: (169.9347618, 0.0304),
  171: (170.9363258, 0.1428),
  172: (171.9363815, 0.2183),
  173: (172.9382108, 0.1613),
  174: (173.9388621, 0.3183),
  175: (174.9412765, 0.0),
  176: (175.9425717, 0.1276),
  177: (176.9452608, 0.0),
  178: (177.946647, 0.0),
  179: (178.95017, 0.0),
  180: (179.95233, 0.0),
  181: (180.95615, 0.0)},
 'Zn': {0: (63.9291422, 1.0),
  54: (53.99295, 0.0),
  55: (54.98398, 0.0),
  56: (55.97238, 0.0),
  57: (56.96479, 0.0),
  58: (57.95459, 0.0),
  59: (58.94926, 0.0),
  60: (59.941827, 0.0),
  61: (60.939511, 0.0),
  62: (61.93433, 0.0),
  63: (62.9332116, 0.0),
  64: (63.9291422, 0.48268),
  65: (64.929241, 0.0),
  66: (65.9260334, 0.27975),
  67: (66.9271273, 0.04102),
  68: (67.9248442, 0.19024),
  69: (68.9265503, 0.0),
  70: (69.9253193, 0.00631),
  71: (70.927722, 0.0),
  72: (71.926858, 0.0),
  73: (72.92978, 0.0),
  74: (73.92946, 0.0),
  75: (74.93294, 0.0),
  76: (75.93329, 0.0),
  77: (76.93696, 0.0),
  78: (77.93844, 0.0),
  79: (78.94265, 0.0),
  80: (79.94434, 0.0),
  81: (80.95048, 0.0),
  82: (81.95442, 0.0),
  83: (82.96103, 0.0)},
 'Zr': {0: (89.9047044, 1.0),
  78: (77.95523, 0.0),
  79: (78.94916, 0.0),
  80: (79.9404, 0.0),
  81: (80.93721, 0.0),
  82: (81.93109, 0.0),
  83: (82.92865, 0.0),
  84: (83.92325, 0.0),
  85: (84.92147, 0.0),
  86: (85.91647, 0.0),
  87: (86.914816, 0.0),
  88: (87.910227, 0.0),
  89: (88.90889, 0.0),
  90: (89.9047044, 0.5145),
  91: (90.9056458, 0.1122),
  92: (91.9050408, 0.1715),
  93: (92.906476, 0.0),
  94: (93.9063152, 0.1738),
  95: (94.9080426, 0.0),
  96: (95.9082734, 0.028),
  97: (96.9109531, 0.0),
  98: (97.912735, 0.0),
  99: (98.916512, 0.0),
  100: (99.91776, 0.0),
  101: (100.92114, 0.0),
  102: (101.92298, 0.0),
  103: (102.9266, 0.0),
  104: (103.92878, 0.0),
  105: (104.93305, 0.0),
  106: (105.93591, 0.0),
  107: (106.94075, 0.0),
  108: (107.94396, 0.0),
  109: (108.94924, 0.0),
  110: (109.95287, 0.0)},
 'e*': {0: (0.00054857990943, 1.0)}}
