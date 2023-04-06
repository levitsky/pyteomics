"""
peff - PSI Extended FASTA Format
================================

PEFF is a forth-coming standard from PSI-HUPO formalizing and extending the
encoding of protein features and annotations for building search spaces for
proteomics. See `The PEFF specification <http://www.psidev.info/peff>`_ for
more up-to-date information on the standard.

Data manipulation
-----------------

Classes
.......

The PEFF parser inherits several properties from implementation in the :mod:`~.fasta` module,
building on top of the :class:`~.TwoLayerIndexedFASTA` reader.

Available classes:

  :py:class:`IndexedPEFF` - Parse a PEFF format file in binary-mode, supporting
  direct indexing by header string or by tag.

"""

#   Copyright 2018 Joshua Klein, Lev Levitsky
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
try:
    from collections.abc import Sequence as SequenceABC, Mapping
except ImportError:
    from collections import Sequence as SequenceABC, Mapping
from collections import OrderedDict, defaultdict

from .fasta import TwoLayerIndexedFASTA


class Header(Mapping):
    """Hold parsed properties of a key-value pair like a sequence's
    definition line.

    This object supports the :class:`Mapping` interface, and
    keys may be accessed by attribute access notation.
    """
    def __init__(self, mapping, original=None):
        self._mapping = mapping

    def __getitem__(self, key):
        return self._mapping[key]

    def __iter__(self):
        return iter(self._mapping)

    def items(self):
        return self._mapping.items()

    def keys(self):
        return self._mapping.keys()

    def values(self):
        return self._mapping.values()

    def __len__(self):
        return len(self._mapping)

    def __contains__(self, key):
        return key in self._mapping

    def __getattr__(self, key):
        if key == "_mapping":
            raise AttributeError(key)
        try:
            return self._mapping[key]
        except KeyError:
            raise AttributeError(key)

    def __repr__(self):
        return "{self.__class__.__name__}({mapping})".format(self=self, mapping=dict(self._mapping))

    def __hash__(self):
        return hash(self.defline)

    def __eq__(self, other):
        try:
            return self._mapping == other._mapping
        except AttributeError:
            return str(self) == str(other)

    def __ne__(self, other):
        return not (self == other)

    def __dir__(self):
        base = set(dir(super(Header, self)))
        keys = set(self._mapping.keys())
        return list(base | keys)


class IndexedPEFF(TwoLayerIndexedFASTA):
    """Creates an :py:class:`IndexedPEFF` object.

    Parameters
    ----------
    source : str or file
        The file to read. If a file object, it needs to be in *rb* mode.
    parse : bool, optional
        Defines whether the descriptions should be parsed in the produced tuples.
        Default is :py:const:`True`.
    kwargs : passed to the :py:class:`TwoLayerIndexedFASTA` constructor.
    """

    kv_pattern = re.compile(r"\\(?P<key>\S+)=(?P<value>.+?)(?:\s(?=\\)|$)")
    header_pattern = re.compile(r"^>?(\S+):(\S+)")
    has_feature_index = re.compile(r"^\(?(\d+):")
    header_group = 2

    class _PEFFFeature(SequenceABC):
        def __init__(self, *fields, **kwargs):
            self.fields = tuple(fields)
            self.id = kwargs.get('id')
            self.feature_type = kwargs.get("feature_type")

        def __eq__(self, other):
            return tuple(self) == tuple(other)

        def __ne__(self, other):
            return not (self == other)

        def __getitem__(self, i):
            return self.fields[i]

        def __len__(self):
            return len(self.fields)

        def __repr__(self):
            return repr(tuple(self))

        def __str__(self):
            return "(%s%s)" % (
                '%r:' % self.id if self.id is not None else '',
                '|'.join(map(str, self)), )

    def __init__(self, source, ignore_comments=False, **kwargs):
        super(IndexedPEFF, self).__init__(
            source, ignore_comments=ignore_comments, parser=self.parser,
            header_pattern=self.header_pattern, **kwargs)
        self.header_blocks = []
        self.comments = []
        self.version = None
        self.number_of_entries = 0
        self._parse_header()

    def _parse_header(self):
        self.seek(0)
        line = self.readline().decode("ascii")
        if not line.startswith("# PEFF"):
            raise ValueError("Not a PEFF File")
        self.version = tuple(map(int, line.strip()[7:].split(".")))
        current_block = defaultdict(list)
        in_header = True
        while in_header:
            line = self.readline().decode("ascii")
            if not line.startswith("#"):
                in_header = False
            line = line.strip()[2:]
            if '=' in line:
                key, value = line.split("=", 1)
                if key == "GeneralComment":
                    self.comments.append(value)
                else:
                    current_block[key].append(value)
            if line.startswith("//"):
                if current_block:
                    self.header_blocks.append(
                        Header(OrderedDict((k, v if len(v) > 1 else v[0])
                                           for k, v in current_block.items())))
                current_block = defaultdict(list)
        number_of_entries = 0
        for block in self.header_blocks:
            try:
                number_of_entries += int(block['NumberOfEntries'])
            except KeyError:
                pass
        self.number_of_entries = number_of_entries

    def _extract_parenthesis_list(self, text):
        chunks = []
        chunk = []
        paren_level = 0
        i = 0
        n = len(text)
        while i < n:
            c = text[i]
            i += 1
            if c == "(":
                if paren_level > 0:
                    chunk.append(c)
                paren_level += 1
            elif c == ")":
                if paren_level > 1:
                    chunk.append(c)
                paren_level -= 1
                if paren_level == 0:
                    if chunk:
                        chunks.append(chunk)
                    chunk = []
            else:
                chunk.append(c)
        chunks = list(map(''.join, chunks))
        return chunks

    def _split_pipe_separated_tuple(self, text):
        parts = text.split("|")
        return parts

    def _coerce_types(self, key, value):
        value = value.strip()
        feature_id_match = self.has_feature_index.search(value)
        if feature_id_match:
            feature_id = int(feature_id_match.group(1))
            value = self.has_feature_index.sub('', value)
        else:
            feature_id = None
        if "|" in value:
            value = self._split_pipe_separated_tuple(value)
            result = []
            for i, v in enumerate(value):
                result.append(self._coerce_value(key, v, i))
            return self._PEFFFeature(*result, feature_type=key, id=feature_id)
        else:
            return self._coerce_value(key, value, 0)

    def _coerce_value(self, key, value, index):
        try:
            return int(value)
        except ValueError:
            pass
        try:
            return float(value)
        except ValueError:
            pass
        return str(value)

    def parser(self, line):
        match = self.header_pattern.match(line)
        if not match:
            raise ValueError(
                "Failed to parse {!r} using {!r}".format(
                    line, self))
        storage = OrderedDict()
        prefix = None
        db_uid = None
        if line.startswith(">"):
            line = line[1:]
        prefix, line = line.split(":", 1)
        db_uid, line = line.split(" ", 1)
        storage['Prefix'] = prefix
        storage['Tag'] = db_uid
        kv_pattern = re.compile(r"\\(?P<key>\S+)=(?P<value>.+?)(?:\s(?=\\)|$)")
        for key, value in kv_pattern.findall(line):
            if not (value.startswith("(") or " (" in value):
                storage[key] = self._coerce_types(key, value)
            else:
                # multi-value
                storage[key] = [self._coerce_types(key, v) for v in self._extract_parenthesis_list(value)]
        return Header(storage)
