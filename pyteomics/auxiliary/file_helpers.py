import sys
import codecs

from functools import wraps
from contextlib import contextmanager


try:
    basestring
except NameError:
    basestring = (str, bytes)

try:
    import pandas as pd
except ImportError:
    pd = None

try:
    import numpy as np
except ImportError:
    np = None


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


def _keepstate_method(func):
    """Decorator for :py:class:`FileReader` methods to help keep the position
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


class _file_obj(object):
    """Check if `f` is a file name and open the file in `mode`.
    A context manager."""

    def __init__(self, f, mode, encoding=None):
        if f is None:
            self.file = {'r': sys.stdin, 'a': sys.stdout, 'w': sys.stdout
                         }[mode[0]]
            self.none = True
        elif isinstance(f, basestring):
            self.file = codecs.open(f, mode, encoding)
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
        except Exception:
            self.__exit__(*sys.exc_info())
            raise

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        pass

    def __iter__(self):
        return self

    def __next__(self):
        # try:
        return next(self._reader)
        # except StopIteration:
        # self.__exit__(None, None, None)
        # raise

    next = __next__


class FileReader(IteratorContextManager):
    """Abstract class implementing context manager protocol
    for file readers.
    """

    def __init__(self, source, mode, func, pass_file, args, kwargs, encoding=None):
        super(FileReader, self).__init__(func, *args, **kwargs)
        self._pass_file = pass_file
        self._source_init = source
        self._mode = mode
        self._encoding = encoding
        self.reset()

    def reset(self):
        if hasattr(self, '_source'):
            self._source.__exit__(None, None, None)
        self._source = _file_obj(self._source_init, self._mode, self._encoding)
        try:
            if self._pass_file:
                self._reader = self._func(
                    self._source, *self._args, **self._kwargs)
            else:
                self._reader = self._func(*self._args, **self._kwargs)
        except Exception:  # clean up on any error
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
            if args:
                return FileReader(args[0], _mode, _func, True, args[1:], kwargs, kwargs.pop('encoding', None))
            source = kwargs.pop('source', None)
            return FileReader(source, _mode, _func, True, (), kwargs, kwargs.pop('encoding', None))
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
