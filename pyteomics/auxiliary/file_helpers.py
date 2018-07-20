import sys
import codecs
import re
from functools import wraps
from contextlib import contextmanager
from collections import OrderedDict
import warnings
warnings.formatwarning = lambda msg, *args, **kw: str(msg) + '\n'

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

def remove_bom(bstr):
    return bstr.replace(codecs.BOM_LE, b'').lstrip(b"\x00")

class IndexedTextReader(FileReader):
    """Abstract class for text file readers that keep an index of records for random access.
    This requires reading the file in binary mode."""

    delimiter = None
    label = None
    block_size = 1000000
    label_group = 1

    def __init__(self, source, func, pass_file, args, kwargs, encoding='utf-8', block_size=None,
        delimiter=None, label=None, label_group=None):
        # the underlying _file_obj gets None as encoding
        # to avoid transparent decoding of StreamReader on read() calls
        super(IndexedTextReader, self).__init__(source, 'rb', func, pass_file, args, kwargs, encoding=None)
        self.encoding = encoding
        if delimiter is not None:
            self.delimiter = delimiter
        if label is not None:
            self.label = label
        if block_size is not None:
            self.block_size = block_size
        if label_group is not None:
            self.label_group = label_group
        self._offset_index = self.build_byte_index()

    def _chunk_iterator(self):
        fh = self._source.file
        delim = remove_bom(self.delimiter.encode(self.encoding))
        buff = fh.read(self.block_size)
        parts = buff.split(delim)
        started_with_delim = buff.startswith(delim)
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
            buff = fh.read(self.block_size)
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
        yield delim + tail

    def _generate_offsets(self):
        i = 0
        pattern = re.compile(remove_bom(self.label.encode(self.encoding)))
        for chunk in self._chunk_iterator():
            match = pattern.search(chunk)
            if match:
                label = match.group(self.label_group)
                yield i, label.decode(self.encoding), match
            i += len(chunk)
        yield i, None, None

    def build_byte_index(self):
        index = OrderedDict()
        g = self._generate_offsets()
        last_offset = 0
        last_label = None
        for offset, label, keyline in g:
            if last_label is not None:
                index[last_label] = (last_offset, offset)
            last_label = label
            last_offset = offset
        assert last_label is None
        return index

    def _read_lines_from_offsets(self, start, end):
        self._source.seek(start)
        lines = self._source.read(end-start).decode(self.encoding).split('\n')
        return lines

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
            enc = kwargs.pop('encoding', None)
            if len(args) > 1:
                with _file_obj(args[1], m, encoding=enc) as out:
                    return _func(args[0], out, *args[2:], **kwargs)
            else:
                with _file_obj(kwargs.pop('output', None), m, encoding=enc) as out:
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

def _check_use_index(source, use_index, default):
    if use_index is not None:
        use_index = bool(use_index)
    if 'b' not in getattr(source, 'mode', 'b'):
        if use_index is True:
            warnings.warn('use_index is True, but the file mode is not binary. '
                'Setting use_index to False')
        use_index = False
    elif 'b' in getattr(source, 'mode', ''):
        if use_index is False:
            warnings.warn('use_index is False, but the file mode is binary. '
                'Setting use_index to True')
        use_index = True
    if use_index is None:
        use_index = default
    return use_index