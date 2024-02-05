import sys
import codecs
import re
from functools import wraps
from contextlib import contextmanager
from collections import OrderedDict, defaultdict
import json
import multiprocessing as mp
import threading
import warnings
import os
from abc import ABCMeta

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

try:
    import dill
except ImportError:
    dill = None
    try:
        import cPickle as pickle
    except ImportError:
        import pickle
    serializer = pickle
else:
    serializer = dill

try:
    from queue import Empty
except ImportError:
    from Queue import Empty

try:
    from collections.abc import Sequence
except ImportError:
    from collections import Sequence

from .structures import PyteomicsError
from .utils import add_metaclass


def _keepstate(func):
    """Decorator to help keep the position in open files passed as
    positional arguments to functions"""
    @wraps(func)
    def wrapped(*args, **kwargs):
        positions = [getattr(arg, 'seek', None) and getattr(arg, 'tell', type(None))() for arg in args]
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
        self._file_spec = None
        self.mode = mode
        if f is None:
            self.file = {'r': sys.stdin, 'a': sys.stdout, 'w': sys.stdout
                         }[mode[0]]
            self._file_spec = None
        elif isinstance(f, basestring):
            self.file = codecs.open(f, mode, encoding)
            self._file_spec = f
        else:
            self._file_spec = f
            self.file = f
        self.encoding = getattr(self.file, 'encoding', encoding)
        self.close_file = (self.file is not f)

    def __enter__(self):
        return self

    def __reduce_ex__(self, protocol):
        return self.__class__, (self._file_spec, self.mode, self.encoding)

    def __exit__(self, *args, **kwargs):
        if (not self.close_file) or self._file_spec is None:
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


class NoOpBaseReader(object):
    def __init__(self, *args, **kwargs):
        pass


class IteratorContextManager(NoOpBaseReader):
    def __init__(self, *args, **kwargs):
        self._func = kwargs.pop('parser_func')
        self._args = args
        self._kwargs = kwargs
        if type(self) == IteratorContextManager:
            self.reset()
        super(IteratorContextManager, self).__init__(*args, **kwargs)

    def __getstate__(self):
        state = {}
        state['_iterator_args'] = self._args
        state['_iterator_kwargs'] = self._kwargs
        return state

    def __setstate__(self, state):
        self._args = state['_iterator_args']
        self._kwargs = state['_iterator_kwargs']

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


@add_metaclass(ABCMeta)
class FileReader(IteratorContextManager):
    """Abstract class implementing context manager protocol
    for file readers.
    """

    def __init__(self, source, **kwargs):
        func = kwargs['parser_func']
        super(FileReader, self).__init__(*kwargs['args'], parser_func=func, **kwargs['kwargs'])
        self._pass_file = kwargs['pass_file']
        self._source_init = source
        self._mode = kwargs['mode']
        self._encoding = kwargs.get('encoding')
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


class IndexedReaderMixin(NoOpBaseReader):
    """Common interface for :py:class:`IndexedTextReader` and :py:class:`IndexedXML`."""
    @property
    def index(self):
        return self._offset_index

    @property
    def default_index(self):
        return self._offset_index

    def __len__(self):
        return len(self._offset_index)

    def __contains__(self, key):
        return key in self._offset_index

    def _item_from_offsets(self, offsets):
        raise NotImplementedError

    def get_by_id(self, elem_id):
        index = self.default_index
        if index is None:
            raise PyteomicsError('Access by ID requires building an offset index.')
        offsets = index[elem_id]
        return self._item_from_offsets(offsets)

    def get_by_ids(self, ids):
        return [self.get_by_id(key) for key in ids]

    def get_by_index(self, i):
        try:
            key = self.default_index.from_index(i, False)
        except AttributeError:
            raise PyteomicsError('Positional access requires building an offset index.')
        return self.get_by_id(key)

    def get_by_indexes(self, indexes):
        return [self.get_by_index(i) for i in indexes]

    def get_by_index_slice(self, s):
        try:
            keys = self.default_index.from_slice(s, False)
        except AttributeError:
            raise PyteomicsError('Positional access requires building an offset index.')
        return self.get_by_ids(keys)

    def get_by_key_slice(self, s):
        keys = self.default_index.between(s.start, s.stop)
        if s.step:
            keys = keys[::s.step]
        return self.get_by_ids(keys)

    def __getitem__(self, key):
        if isinstance(key, basestring):
            return self.get_by_id(key)
        if isinstance(key, int):
            return self.get_by_index(key)
        if isinstance(key, Sequence):
            if not key:
                return []
            if isinstance(key[0], int):
                return self.get_by_indexes(key)
            if isinstance(key[0], basestring):
                return self.get_by_ids(key)
        if isinstance(key, slice):
            for item in (key.start, key.stop, key.step):
                if item is not None:
                    break
            if isinstance(item, int):
                return self.get_by_index_slice(key)
            if isinstance(item, basestring):
                return self.get_by_key_slice(key)
            if item is None:
                return list(self)
        raise PyteomicsError('Unsupported query key: {}'.format(key))


class RTLocator():
    def __init__(self, reader):
        self._reader = reader

    def _get_scan_by_time(self, time):
        """Retrieve the scan object for the specified scan time.

        Parameters
        ----------
        time : float
            The time to get the nearest scan from
        Returns
        -------
        tuple: (scan_id, scan, scan_time)
        """
        if not self._reader.default_index:
            raise PyteomicsError("This method requires the index. Please pass `use_index=True` during initialization")

        scan_ids = tuple(self._reader.default_index)
        lo = 0
        hi = len(scan_ids)

        best_match = None
        best_error = float('inf')
        best_time = None
        best_id = None

        if time == float('inf'):
            scan = self._reader.get_by_id(scan_ids[-1])
            return scan_ids[-1], scan, self._reader._get_time(scan)

        while hi != lo:
            mid = (hi + lo) // 2
            sid = scan_ids[mid]
            scan = self._reader.get_by_id(sid)
            scan_time = self._reader._get_time(scan)
            err = abs(scan_time - time)
            if err < best_error:
                best_error = err
                best_match = scan
                best_time = scan_time
                best_id = sid
            if scan_time == time:
                return sid, scan, scan_time
            elif (hi - lo) == 1:
                return best_id, best_match, best_time
            elif scan_time > time:
                hi = mid
            else:
                lo = mid

    def __getitem__(self, key):
        if isinstance(key, (int, float)):
            return self._get_scan_by_time(key)[1]
        if isinstance(key, Sequence):
            return [self._get_scan_by_time(t)[1] for t in key]
        if isinstance(key, slice):
            if key.start is None:
                start_index = self._reader.default_index.from_index(0)
            else:
                start_index = self._get_scan_by_time(key.start)[0]
            if key.stop is None:
                stop_index = self._reader.default_index.from_index(-1)
            else:
                stop_index = self._get_scan_by_time(key.stop)[0]
            return self._reader[start_index:stop_index:key.step]


class TimeOrderedIndexedReaderMixin(IndexedReaderMixin):
    @property
    def time(self):
        return self._time

    def __init__(self, *args, **kwargs):
        super(TimeOrderedIndexedReaderMixin, self).__init__(*args, **kwargs)
        self._time = RTLocator(self)

    @staticmethod
    def _get_time(scan):
        raise NotImplementedError


class IndexedTextReader(IndexedReaderMixin, FileReader):
    """Abstract class for text file readers that keep an index of records for random access.
    This requires reading the file in binary mode."""

    delimiter = None
    label = None
    block_size = 1000000
    label_group = 1
    _kw_keys = ['delimiter', 'label', 'block_size', 'label_group']
    _warn_if_empty = True

    def __init__(self, source, **kwargs):
        # the underlying _file_obj gets None as encoding
        # to avoid transparent decoding of StreamReader on read() calls
        encoding = kwargs.pop('encoding', 'utf-8')
        if 'warn_if_empty' in kwargs:
            self._warn_if_empty = kwargs.pop('warn_if_empty')
        super(IndexedTextReader, self).__init__(source, mode='rb', encoding=None, **kwargs)
        self.encoding = encoding
        for attr in self._kw_keys:
            if attr in kwargs:
                setattr(self, attr, kwargs.pop(attr))
        self._offset_index = None
        if not kwargs.pop('_skip_index', False):
            self._offset_index = self.build_byte_index()

    def __getstate__(self):
        state = super(IndexedTextReader, self).__getstate__()
        state['offset_index'] = self._offset_index
        for key in self._kw_keys:
            state[key] = getattr(self, key)
        return state

    def __setstate__(self, state):
        super(IndexedTextReader, self).__setstate__(state)
        self._offset_index = state['offset_index']
        for key in self._kw_keys:
            if key in state:
                setattr(self, key, state[key])

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
        index = OffsetIndex()
        g = self._generate_offsets()
        last_offset = 0
        last_label = None
        for offset, label, keyline in g:
            if last_label is not None:
                index[last_label] = (last_offset, offset)
            last_label = label
            last_offset = offset
        assert last_label is None
        if self._warn_if_empty and not index:
            self._warn_empty()
        return index

    def _read_lines_from_offsets(self, start, end):
        self._source.seek(start)
        lines = self._source.read(end - start).decode(self.encoding).split('\n')
        return lines

    def _warn_empty(self):
        warnings.warn("{} object has an empty index for file {}. If this is unexpected, consider adjusting `label`.".format(
            self.__class__.__name__, getattr(self._source, 'name', self._source_init)))


class IndexSavingMixin(NoOpBaseReader):
    """Common interface for :py:class:`IndexSavingXML` and :py:class:`IndexSavingTextReader`."""
    _index_class = NotImplemented

    @property
    def _byte_offset_filename(self):
        try:
            path = self._source.name
        except AttributeError:
            return None
        name, ext = os.path.splitext(path)
        byte_offset_filename = '{}-{}-byte-offsets.json'.format(name, ext[1:])
        return byte_offset_filename

    def _check_has_byte_offset_file(self):
        """Check if the file at :attr:`_byte_offset_filename` exists

        Returns
        -------
        bool
            Whether the file exists
        """
        path = self._byte_offset_filename
        if path is None:
            return False
        return os.path.exists(path)

    @classmethod
    def prebuild_byte_offset_file(cls, path):
        """Construct a new XML reader, build its byte offset index and
        write it to file

        Parameters
        ----------
        path : str
            The path to the file to parse
        """
        with cls(path) as inst:
            inst.write_byte_offsets()

    def write_byte_offsets(self):
        """Write the byte offsets in :attr:`_offset_index` to the file
        at :attr:`_byte_offset_filename`
        """
        with open(self._byte_offset_filename, 'w') as f:
            self._offset_index.save(f)

    @_keepstate_method
    def build_byte_index(self):
        """Build the byte offset index by either reading these offsets
        from the file at :attr:`_byte_offset_filename`, or falling back
        to the method used by :class:`IndexedXML` or :class:`IndexedTextReader` if this operation fails
        due to an IOError
        """
        if not getattr(self, '_use_index', True): return  # indexed text readers do not have `_use_index`
        try:
            return self._read_byte_offsets()
        except (IOError, AttributeError, TypeError):
            return super(IndexSavingMixin, self).build_byte_index()

    def _read_byte_offsets(self):
        """Read the byte offset index JSON file at :attr:`_byte_offset_filename`
        and populate :attr:`_offset_index`
        """
        with open(self._byte_offset_filename, 'r') as f:
            return self._index_class.load(f)


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
                return FileReader(args[0], mode=_mode, parser_func=_func, pass_file=True, args=args[1:], kwargs=kwargs,
                    encoding=kwargs.pop('encoding', None))
            source = kwargs.pop('source', None)
            return FileReader(source, mode=_mode, parser_func=_func, pass_file=True, args=(), kwargs=kwargs, encoding=kwargs.pop('encoding', None))
        return helper
    return decorator


def _file_writer(_mode='w'):
    def decorator(_func):
        """A decorator that opens output files for writer functions.
        """
        @wraps(_func)
        def helper(*args, **kwargs):
            m = kwargs.pop('file_mode', _mode)
            enc = kwargs.pop('encoding', None)
            if len(args) > 1:
                out_arg = args[1]
            else:
                out_arg = kwargs.pop('output', None)

            with _file_obj(out_arg, m, encoding=enc) as out:
                if len(args) > 1:
                    call_args = (args[0], out) + args[2:]
                    call_kwargs = kwargs
                else:
                    call_args = args
                    call_kwargs = dict(output=out, **kwargs)
                return _func(*call_args, **call_kwargs)
        return helper
    return decorator


class WritableIndex(object):
    schema_version = (1, 0, 0)
    _schema_version_tag_key = "@pyteomics_schema_version"

    def _serializable_container(self):
        container = {'index': list(self.items())}
        return container

    def save(self, fp):
        container = self._serializable_container()
        container[self._schema_version_tag_key] = self.schema_version
        json.dump(container, fp)

    @classmethod
    def load(cls, fp):
        container = json.load(fp, object_hook=OrderedDict)
        version_tag = container.get(cls._schema_version_tag_key)
        if version_tag is None:
            # The legacy case, no special processing yet
            inst = cls()
            inst.schema_version = None
            return inst
        version_tag = tuple(version_tag)
        index = container.get("index")
        if version_tag < cls.schema_version:
            # schema upgrade case, no special processing yet
            inst = cls(index)
            inst.schema_version = version_tag
            return inst
        # no need to upgrade
        return cls(index)


class OffsetIndex(OrderedDict, WritableIndex):
    '''An augmented OrderedDict that formally wraps getting items by index
    '''

    def __init__(self, *args, **kwargs):
        super(OffsetIndex, self).__init__(*args, **kwargs)
        self._index_sequence = None

    def _invalidate(self):
        self._index_sequence = None

    @property
    def index_sequence(self):
        """Keeps a cached copy of the :meth:`items` sequence
        stored as a :class:`tuple` to avoid repeatedly copying
        the sequence over many method calls.

        Returns
        -------
        :class:`tuple`
        """
        if self._index_sequence is None:
            self._index_sequence = tuple(self.items())
        return self._index_sequence

    def __setitem__(self, key, value):
        self._invalidate()
        return super(OffsetIndex, self).__setitem__(key, value)

    def pop(self, *args, **kwargs):
        self._invalidate()
        return super(OffsetIndex, self).pop(*args, **kwargs)

    def find(self, key, *args, **kwargs):
        return self[key]

    def from_index(self, index, include_value=False):
        '''Get an entry by its integer index in the ordered sequence
        of this mapping.

        Parameters
        ----------
        index: int
            The index to retrieve.
        include_value: bool
            Whether to return both the key and the value or just the key.
            Defaults to :const:`False`.

        Returns
        -------
        object:
            If ``include_value`` is :const:`True`, a tuple of (key, value) at ``index``
            else just the key at ``index``.
        '''
        items = self.index_sequence
        if include_value:
            return items[index]
        else:
            return items[index][0]

    def from_slice(self, spec, include_value=False):
        '''Get a slice along index in the ordered sequence
        of this mapping.

        Parameters
        ----------
        spec: slice
            The slice over the range of indices to retrieve
        include_value: bool
            Whether to return both the key and the value or just the key.
            Defaults to :const:`False`

        Returns
        -------
        list:
            If ``include_value`` is :const:`True`, a tuple of (key, value) at ``index``
            else just the key at ``index`` for each ``index`` in ``spec``
        '''
        items = self.index_sequence
        return [(k, v) if include_value else k for k, v in items[spec]]

    def between(self, start, stop, include_value=False):
        keys = list(self)
        if start is not None:
            try:
                start_index = keys.index(start)
            except ValueError:
                raise KeyError(start)
        else:
            start_index = 0
        if stop is not None:
            try:
                stop_index = keys.index(stop)
            except ValueError:
                raise KeyError(stop)
        else:
            stop_index = len(keys) - 1
        if start is None or stop is None:
            pass  # won't switch indices
        else:
            start_index, stop_index = min(start_index, stop_index), max(start_index, stop_index)

        if include_value:
            return [(k, self[k]) for k in keys[start_index:stop_index + 1]]
        return keys[start_index:stop_index + 1]

    def __repr__(self):
        template = "{self.__class__.__name__}({items})"
        return template.format(self=self, items=list(self.items()))

    def _integrity_check(self):
        indices = list(self.values())
        sorted_indices = sorted(self.values())
        return indices == sorted_indices

    def sort(self):
        sorted_pairs = sorted(self.items(), key=lambda x: x[1])
        self.clear()
        self._invalidate()
        for key, value in sorted_pairs:
            self[key] = value
        return self


class IndexSavingTextReader(IndexSavingMixin, IndexedTextReader):
    _index_class = OffsetIndex


class HierarchicalOffsetIndex(WritableIndex):
    _inner_type = OffsetIndex

    def __init__(self, base=None):
        self.mapping = defaultdict(self._inner_type)
        for key, value in (base or {}).items():
            self.mapping[key] = self._inner_type(value)

    def _integrity_check(self):
        for key, value in self.items():
            if not value._integrity_check():
                return False
        return True

    def sort(self):
        for key, value in self.items():
            value.sort()
        return self

    def __getitem__(self, key):
        return self.mapping[key]

    def __setitem__(self, key, value):
        self.mapping[key] = value

    def __iter__(self):
        return iter(self.mapping)

    def __len__(self):
        return sum(len(group) for key, group in self.items())

    def __contains__(self, key):
        return key in self.mapping

    def find(self, key, element_type=None):
        if element_type is None:
            for element_type in self.keys():
                try:
                    return self.find(key, element_type)
                except KeyError:
                    continue
            raise KeyError(key)
        else:
            return self[element_type][key]

    def find_no_type(self, key):
        """Try to find `key` in each of the lower-level indexes, returning both
        value and the element type that match the key."""
        for element_type in self.keys():
            try:
                return self.find(key, element_type), element_type
            except KeyError:
                continue
        raise KeyError(key)

    def update(self, *args, **kwargs):
        self.mapping.update(*args, **kwargs)

    def pop(self, key, default=None):
        return self.mapping.pop(key, default)

    def keys(self):
        return self.mapping.keys()

    def values(self):
        return self.mapping.values()

    def items(self):
        return self.mapping.items()

    def _serializable_container(self):
        encoded_index = {}
        container = {
            'keys': list(self.keys())
        }
        for key, offset in self.items():
            encoded_index[key] = list(offset.items())
        container['index'] = encoded_index
        return container


def _make_chain(reader, readername, full_output=False):

    def concat_results(*args, **kwargs):
        results = [reader(arg, **kwargs) for arg in args]
        if pd is not None and all(isinstance(a, pd.DataFrame) for a in args):
            return pd.concat(results)
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
    try:
        if use_index is not None:
            use_index = bool(use_index)

        # if a file name is given, do not override anything; short-circuit
        if isinstance(source, basestring):
            return use_index if use_index is not None else default

        # collect information on source
        if hasattr(source, 'seekable'):
            seekable = source.seekable()
        else:
            seekable = None

        if hasattr(source, 'mode'):
            binary = 'b' in source.mode
        else:
            binary = None

        # now check for conflicts
        if seekable is False:
            if binary:
                raise PyteomicsError('Cannot work with non-seekable file in binary mode: {}.'.format(source))
            if use_index:
                warnings.warn('Cannot use indexing as {} is not seekable. Setting `use_index` to False.'.format(source))
                use_index = False
        elif binary is not None:
            if use_index is not None and binary != use_index:
                warnings.warn('use_index is {}, but the file mode is {}. '
                    'Setting `use_index` to {}'.format(use_index, source.mode, binary))
            use_index = binary
        elif use_index is None:
            warnings.warn('Could not check mode on {}. Specify `use_index` explicitly to avoid errors.'.format(source))

        if use_index is not None:
            return use_index

        return default

    except PyteomicsError:
        raise
    except Exception as e:
        if use_index is None:
            warnings.warn('Could not check mode on {}. Reason: {!r}. '
                'Specify `use_index` explicitly to avoid errors.'.format(source, e))
            return default
        return use_index


class FileReadingProcess(mp.Process):
    """Process that does a share of distributed work on entries read from file.
    Reconstructs a reader object, parses an entries from given indexes,
    optionally does additional processing, sends results back.

    The reader class must support the :py:meth:`__getitem__` dict-like lookup.
    """

    def __init__(self, reader_spec, target_spec, qin, qout, args_spec, kwargs_spec):
        super(FileReadingProcess, self).__init__(name='pyteomics-map-worker')
        self.reader_spec = reader_spec
        self.target_spec = target_spec
        self.args_spec = args_spec
        self.kwargs_spec = kwargs_spec
        self._qin = qin
        self._qout = qout
        # self._in_flag = in_flag
        self._done_flag = mp.Event()
        self.daemon = True

    def run(self):
        reader = serializer.loads(self.reader_spec)
        target = serializer.loads(self.target_spec)
        args = serializer.loads(self.args_spec)
        kwargs = serializer.loads(self.kwargs_spec)
        for key in iter(self._qin.get, None):
            item = reader[key]
            if target is not None:
                result = target(item, *args, **kwargs)
            else:
                result = item
            self._qout.put(result)
        self._done_flag.set()

    def is_done(self):
        return self._done_flag.is_set()


try:
    _NPROC = mp.cpu_count()
except NotImplementedError:
    _NPROC = 4
_QUEUE_TIMEOUT = 4
_QUEUE_SIZE = int(1e7)


class TaskMappingMixin(NoOpBaseReader):
    def __init__(self, *args, **kwargs):
        '''
        Instantiate a :py:class:`TaskMappingMixin` object, set default parameters for IPC.

        Parameters
        ----------

        queue_timeout : float, keyword only, optional
            The number of seconds to block, waiting for a result before checking to see if
            all workers are done.
        queue_size : int, keyword only, optional
            The length of IPC queue used.
        processes : int, keyword only, optional
            Number of worker processes to spawn when :py:meth:`map` is called. This can also be
            specified in the :py:meth:`map` call.
        '''
        self._queue_size = kwargs.pop('queue_size', _QUEUE_SIZE)
        self._queue_timeout = kwargs.pop('timeout', _QUEUE_TIMEOUT)
        self._nproc = kwargs.pop('processes', _NPROC)
        super(TaskMappingMixin, self).__init__(*args, **kwargs)

    def _get_reader_for_worker_spec(self):
        return self

    def _build_worker_spec(self, target, args, kwargs):
        serialized = []
        for obj, objname in [(self._get_reader_for_worker_spec(), 'reader'), (target, 'target'), (args, 'args'),
                             (kwargs, 'kwargs')]:
            try:
                serialized.append(serializer.dumps(obj))
            except serializer.PicklingError:
                msg = 'Could not serialize {0} {1} with {2.__name__}.'.format(objname, obj, serializer)
                if serializer is not dill:
                    msg += ' Try installing `dill`.'
                raise PyteomicsError(msg)
        return serialized

    def _spawn_workers(self, specifications, in_queue, out_queue, processes):
        reader_spec, target_spec, args_spec, kwargs_spec = specifications
        workers = []
        for _ in range(processes):
            worker = FileReadingProcess(
                reader_spec, target_spec, in_queue, out_queue, args_spec, kwargs_spec)
            workers.append(worker)
        return workers

    def _spawn_feeder_thread(self, in_queue, iterator, processes):
        def feeder():
            for key in iterator:
                in_queue.put(key)
            for _ in range(processes):
                in_queue.put(None)

        feeder_thread = threading.Thread(target=feeder)
        feeder_thread.daemon = True
        feeder_thread.start()
        return feeder_thread

    def map(self, target=None, processes=-1, args=None, kwargs=None, **_kwargs):
        """Execute the ``target`` function over entries of this object across up to ``processes``
        processes.

        Results will be returned out of order.

        Parameters
        ----------
        target : :class:`Callable`, optional
            The function to execute over each entry. It will be given a single object yielded by
            the wrapped iterator as well as all of the values in ``args`` and ``kwargs``
        processes : int, optional
            The number of worker processes to use. If 0 or negative,
            defaults to the number of available CPUs.
            This parameter can also be set at reader creation.
        args : :class:`Sequence`, optional
            Additional positional arguments to be passed to the target function
        kwargs : :class:`Mapping`, optional
            Additional keyword arguments to be passed to the target function
        **_kwargs
            Additional keyword arguments to be passed to the target function

        Yields
        ------
        object
            The work item returned by the target function.
        """
        if self._offset_index is None:
            raise PyteomicsError('The reader needs an index for map() calls. Create the reader with `use_index=True`.')

        if processes < 1:
            processes = self._nproc
        iterator = self._task_map_iterator()

        if args is None:
            args = tuple()
        else:
            args = tuple(args)
        if kwargs is None:
            kwargs = dict()
        else:
            kwargs = dict(kwargs)
        kwargs.update(_kwargs)

        serialized = self._build_worker_spec(target, args, kwargs)

        in_queue = mp.Queue(self._queue_size)
        out_queue = mp.Queue(self._queue_size)

        workers = self._spawn_workers(serialized, in_queue, out_queue, processes)
        feeder_thread = self._spawn_feeder_thread(in_queue, iterator, processes)
        for worker in workers:
            worker.start()

        def iterate():
            while True:
                try:
                    result = out_queue.get(True, self._queue_timeout)
                    yield result
                except Empty:
                    if all(w.is_done() for w in workers):
                        break
                    else:
                        continue

            feeder_thread.join()
            for worker in workers:
                worker.join()
        return iterate()

    def _task_map_iterator(self):
        """Returns the :class:`Iteratable` to use when dealing work items onto the input IPC
        queue used by :meth:`map`

        Returns
        -------
        :class:`Iteratable`
        """

        return iter(self._offset_index.keys())


class ChainBase(object):
    """Chain :meth:`sequence_maker` for several sources into a
    single iterable. Positional arguments should be sources like
    file names or file objects. Keyword arguments are passed to
    the :meth:`sequence_maker` function.

    Parameters
    ----------
    sources : :class:`Iterable`
        Sources for creating new sequences from, such as paths or
        file-like objects
    kwargs : :class:`Mapping`
        Additional arguments used to instantiate each sequence
    """

    def __init__(self, *sources, **kwargs):
        self.sources = sources
        self.kwargs = kwargs
        self._iterator = None

    @classmethod
    def from_iterable(cls, sources, **kwargs):
        return cls(*sources, **kwargs)

    @classmethod
    def _make_chain(cls, sequence_maker):
        if isinstance(sequence_maker, type):
            tp = type('%sChain' % sequence_maker.__class__.__name__, (cls,), {
                'sequence_maker': sequence_maker,
                '__doc__': cls.__doc__.replace(':meth:`sequence_maker`', ':class:`{}`'.format(sequence_maker.__name__))
            })
        else:
            tp = type('FunctionChain', (cls,), {
                'sequence_maker': staticmethod(sequence_maker),
                '__doc__': cls.__doc__.replace(':meth:`sequence_maker`', ':func:`{}`'.format(sequence_maker.__name__))
            })
        return tp

    def sequence_maker(self, file):
        raise NotImplementedError()

    def _create_sequence(self, file):
        return self.sequence_maker(file, **self.kwargs)

    def _iterate_over_series(self):
        for f in self.sources:
            with self._create_sequence(f) as r:
                for item in r:
                    yield item

    def __enter__(self):
        self._iterator = iter(self._iterate_over_series())
        return self

    def __exit__(self, *args, **kwargs):
        self._iterator = None

    def __iter__(self):
        return self

    def __next__(self):
        if self._iterator is None:
            self._iterator = self._iterate_over_series()
        return next(self._iterator)

    def next(self):
        return self.__next__()

    def map(self, target=None, processes=-1, queue_timeout=_QUEUE_TIMEOUT, args=None, kwargs=None, **_kwargs):
        """Execute the ``target`` function over entries of this object across up to ``processes``
        processes.

        Results will be returned out of order.

        Parameters
        ----------
        target : :class:`Callable`, optional
            The function to execute over each entry. It will be given a single object yielded by
            the wrapped iterator as well as all of the values in ``args`` and ``kwargs``
        processes : int, optional
            The number of worker processes to use. If negative, the number of processes
            will match the number of available CPUs.
        queue_timeout : float, optional
            The number of seconds to block, waiting for a result before checking to see if
            all workers are done.
        args : :class:`Sequence`, optional
            Additional positional arguments to be passed to the target function
        kwargs : :class:`Mapping`, optional
            Additional keyword arguments to be passed to the target function
        **_kwargs
            Additional keyword arguments to be passed to the target function

        Yields
        ------
        object
            The work item returned by the target function.
        """
        for f in self.sources:
            with self._create_sequence(f) as r:
                for result in r.map(target, processes, queue_timeout, args, kwargs, **_kwargs):
                    yield result


class TableJoiner(ChainBase):
    def concatenate(self, results):
        if pd is not None and all(isinstance(a, pd.DataFrame) for a in results):
            return pd.concat(results)
        if isinstance(results[0], np.ndarray):
            return np.concatenate(results)
        else:
            return np.array([b for a in results for b in a])

    def _iterate_over_series(self):
        results = [self._create_sequence(f) for f in self.sources]
        return self.concatenate(results)
