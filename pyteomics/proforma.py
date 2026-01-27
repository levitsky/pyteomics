'''
proforma - Proteoform and Peptidoform Notation
==============================================

ProForma is a notation for defining modified amino acid sequences using
a set of controlled vocabularies, as well as encoding uncertain or partial
information about localization. See `ProForma specification <https://www.psidev.info/proforma>`_
for more up-to-date information.

For more details, see the :mod:`pyteomics.proforma` online.
'''

import itertools
import re
import warnings
from typing import Any, Callable, Dict, Iterable, List, Optional, ClassVar, Sequence, Tuple, Type, Union, Generic, TypeVar, NamedTuple
from collections import deque, namedtuple
from functools import partial
from itertools import chain
from array import array as _array
from enum import Enum
from numbers import Integral

from .mass import Composition, std_aa_mass, Unimod, nist_mass, calculate_mass, std_ion_comp, mass_charge_ratio, std_aa_comp
from .auxiliary import PyteomicsError, BasicComposition
from .auxiliary.utils import add_metaclass
from .auxiliary.psims_util import load_psimod, load_xlmod, load_gno, obo_cache, _has_psims

try:
    import numpy as np
except ImportError:
    np = None


_WATER_MASS = calculate_mass(formula="H2O")

std_aa_mass = std_aa_mass.copy()
std_aa_mass['X'] = 0

element_symbols = set(nist_mass)

T = TypeVar('T')


class ProFormaError(PyteomicsError):
    def __init__(self, message, index=None, parser_state=None, **kwargs):
        super(ProFormaError, self).__init__(PyteomicsError, message, index, parser_state)
        self.message = message
        self.index = index
        self.parser_state = parser_state


class PrefixSavingMeta(type):
    '''A subclass-registering-metaclass that provides easy
    lookup of subclasses by prefix attributes.
    '''

    def __new__(mcs, name, parents, attrs):
        new_type = type.__new__(mcs, name, parents, attrs)
        prefix = attrs.get("prefix_name")
        if prefix:
            new_type.prefix_map[prefix.lower()] = new_type
        short = attrs.get("short_prefix")
        if short:
            new_type.prefix_map[short.lower()] = new_type
        return new_type

    def find_by_tag(self, tag_name):
        if tag_name is None:
            raise ValueError("tag_name cannot be None!")
        tag_name = tag_name.lower()
        return self.prefix_map[tag_name]


class TagTypeEnum(Enum):
    unimod = 0
    psimod = 1
    massmod = 2
    generic = 3
    info = 4
    gnome = 5
    xlmod = 6

    formula = 7
    glycan = 8

    localization_marker = 9
    position_label = 10

    position_modifier = 11
    comup = 12
    comkp = 13
    limit = 14
    custom = 15

    group_placeholder = 999


class ModificationTagStyle(Enum):
    Unset = 0
    ShortId = 1
    LongId = 2
    ShortName = 3
    LongName = 4


class ModificationSourceType(Enum):
    """
    Whether a tag was generated from explicit user input (``Explicit``), a constant
    modification rule (``Constant``), or from a variable expansion (``Generated``).

    Used to track sources in :class:`ProteoformCombinator` machinery.
    """
    Explicit = 0
    Constant = 1
    Generated = 2


_sentinel = object()


class ModificationMassNotFoundError(ProFormaError):
    pass


class CompositionNotFoundError(ProFormaError):
    pass


class MissingChargeStateError(ProFormaError):
    pass


class UnknownMonosaccharideError(ProFormaError):
    pass


@add_metaclass(PrefixSavingMeta)
class TagBase(object):
    '''A base class for all tag types.

    Attributes
    ----------
    type: Enum
        An element of :class:`TagTypeEnum` saying what kind of tag this is.
    value: object
        The data stored in this tag, usually an externally controlled name
    extra: list
        Any extra tags that were nested within this tag. Usually limited to INFO
        tags but may be other synonymous controlled vocabulary terms.
    group_id: str or None
        A short label denoting which group, if any, this tag belongs to
    '''
    __slots__ = ("type", "value", "extra", "group_id", )

    type: TagTypeEnum
    value: Any
    extra: List["TagBase"]
    group_id: Optional[str]

    prefix_name: ClassVar[Optional[str]] = None
    short_prefix: ClassVar[Optional[str]] = None
    prefix_map: ClassVar[Dict[str, Type['TagBase']]] = {}

    def __init__(self, type, value, extra=None, group_id=None):
        self.type = type
        self.value = value
        self.extra = extra or []
        self.group_id = group_id

    def copy(self):
        return self.__class__(self.value, [e.copy() for e in self.extra], self.group_id)

    def __str__(self):
        part = self._format_main()
        had_marker = False
        if self.extra:
            rest = []
            for e in self.extra:
                rest.append(str(e))
                had_marker |= isinstance(e, GroupLabelBase) and e.group_id == self.group_id
            label = '|'.join([part] + rest)
        else:
            label = part
        if self.group_id and not had_marker:
            label = '%s%s' % (label, self.group_id)
        return '%s' % label

    def __repr__(self):
        template = "{self.__class__.__name__}({self.value!r}, {self.extra!r}, {self.group_id!r})"
        return template.format(self=self)

    def __eq__(self, other):
        if other is None:
            return False
        if isinstance(other, str):
            return str(self) == other
        return (self.type == other.type) and (self.value == other.value) and (self.extra == other.extra) \
            and (self.group_id == other.group_id)

    def __ne__(self, other):
        return not self == other

    def is_modification(self) -> bool:
        return self.type in (
            TagTypeEnum.formula,
            TagTypeEnum.generic,
            TagTypeEnum.glycan,
            TagTypeEnum.gnome,
            TagTypeEnum.unimod,
            TagTypeEnum.massmod,
            TagTypeEnum.psimod,
            TagTypeEnum.custom,
        )

    def find_modification(self) -> Optional["TagBase"]:
        if self.is_modification():
            return self
        for tag in self.extra:
            if tag.is_modification:
                return tag
        return None

    def find_tag_type(self, tag_type: TagTypeEnum) -> List['TagBase']:
        '''Search this tag or tag collection for elements with a particular
        tag type and return them.

        Parameters
        ----------
        tag_type : TagTypeEnum
            A label from :class:`TagTypeEnum`, or an equivalent type.

        Returns
        -------
        matches : list
            The list of all tags in this object which match the requested tag type.
        '''
        out = []
        if self.type == tag_type:
            out.append(self)
        if not self.extra:
            return out
        for e in self.extra:
            if e.type == tag_type:
                out.append(e)
        return out

    @classmethod
    def parse(cls, buffer) -> 'TagBase':
        return process_tag_tokens(buffer)

    def has_mass(self) -> bool:
        """
        Check if this tag carries a mass value.

        Returns
        -------
        bool
        """
        return False

    def has_composition(self) -> bool:
        return False


class GroupLabelBase(TagBase):
    __slots__ = ()

    def __str__(self):
        part = self._format_main()
        if self.extra:
            rest = [str(e) for e in self.extra]
            label = '|'.join([part] + rest)
        else:
            label = part
        return '%s' % label


class PositionLabelTag(GroupLabelBase):
    '''A tag to mark that a position is involved in a group in some way, but does
    not imply any specific semantics.
    '''
    __slots__ = ()

    def __init__(self, value=None, extra=None, group_id=None):
        assert group_id is not None
        value = group_id
        super(PositionLabelTag, self).__init__(
            TagTypeEnum.position_label, value, extra, group_id)

    def _format_main(self):
        return "{self.group_id}".format(self=self)


class LocalizationMarker(GroupLabelBase):
    '''A tag to mark a particular localization site
    '''
    __slots__ = ()

    def __init__(self, value, extra=None, group_id=None):
        assert group_id is not None
        super(LocalizationMarker, self).__init__(
            TagTypeEnum.localization_marker, float(value), extra, group_id)

    def _format_main(self):
        return "{self.group_id}({self.value:.4g})".format(self=self)


class InformationTag(TagBase):
    '''A tag carrying free text describing the location
    '''
    __slots__ = ()

    prefix_name = "INFO"

    def __init__(self, value, extra=None, group_id=None):
        super(InformationTag, self).__init__(
            TagTypeEnum.info, str(value), extra, group_id)

    def _format_main(self):
        return f"INFO:{self.value}"


class PositionModifierTag(TagBase):
    __slots__ = ()

    prefix_name = "Position"

    def __init__(self, value, extra=None, group_id=None):
        super().__init__(TagTypeEnum.position_modifier, value, extra, group_id)

    def _format_main(self):
        return f"{self.prefix_name}:{self.value}"


class LimitModifierTag(TagBase):
    __slots__ = ()

    prefix_name = "Limit"

    def __init__(self, value, extra=None, group_id=None):
        if not isinstance(value, (int, float)):
            try:
                value = int(value)
            except (ValueError, TypeError):
                pass
        super().__init__(TagTypeEnum.limit, value, extra, group_id)

    def _format_main(self):
        return f"{self.prefix_name}:{self.value}"


class ColocaliseModificationsOfKnownPostionTag(TagBase):
    __slots__ = ()

    prefix_name = "ColocaliseModificationsOfKnownPosition"
    short_prefix = "CoMKP"

    def __init__(self, extra=None, group_id=None):
        super().__init__(TagTypeEnum.comkp, None, extra, group_id)

    def copy(self):
        return self.__class__([e.copy() for e in (self.extra or [])], self.group_id)

    def _format_main(self):
        return self.short_prefix


class ColocaliseModificationsOfUnknownPostionTag(TagBase):
    __slots__ = ()

    prefix_name = "ColocaliseModificationsOfUnknownPosition"
    short_prefix = "CoMUP"

    def __init__(self, extra=None, group_id=None):
        super().__init__(TagTypeEnum.comup, None, extra, group_id)

    def copy(self):
        return self.__class__([e.copy() for e in self.extra or []], self.group_id)

    def _format_main(self):
        return self.short_prefix


class ModificationResolver(object):
    name: str
    symbol: str

    _database: Optional[Any]
    _cache: Optional[Dict[Tuple[Optional[str], Optional[int], frozenset], Any]]

    def __init__(self, name, **kwargs):
        self.name = name.lower()
        self.symbol = self.name[0]
        self._database = None
        self._cache = {}

    def clear_cache(self):
        """Clear the modification definition cache"""
        self._cache.clear()

    def enable_caching(self, flag: bool=True):
        """
        Enable or disable caching of modification definitions.

        If `flag` is :const:`False`, this will also dispose of any
        existing cached values.

        Parameters
        ----------
        flag : :class:`bool`
            Whether or not to disable the cache
        """
        if flag:
            if not self._cache:
                self._cache = {}
        else:
            self._cache = None

    def load_database(self):
        raise NotImplementedError()

    @property
    def database(self):
        if not self._database:
            self._database = self.load_database()
        return self._database

    @database.setter
    def database(self, database):
        self._database = database

    def parse_identifier(self, identifier: str) -> Tuple[Optional[str], Optional[int]]:
        """Parse a string that is either a CV prefixed identifier or name.

        Parameters
        ----------
        identifier : str
            The identifier string to parse, removing CV prefix as needed.

        Returns
        -------
        name : str, optional
            A textual identifier embedded in the qualified identifier, if any, otherwise
            :const:`None`.
        id : int, optional
            An integer ID embedded in the qualified identifier, if any, otherwise
            :const:`None`.
        """
        tokens = identifier.split(":", 1)
        if len(tokens) > 1:
            prefix = tokens[0].lower()
            if prefix == self.name or prefix == self.symbol:
                identifier = tokens[1]

        if identifier.isdigit():
            id = int(identifier)
            name = None
        else:
            name = identifier
            id = None
        return name, id

    def _resolve_impl(self, name: str=None, id: int=None, **kwargs) -> Dict[str, Any]:
        raise NotImplementedError()

    def resolve(self, name: str=None, id: int=None, **kwargs):
        if self._cache is None:
            return self._resolve_impl(name, id, **kwargs)
        cache_key = (name, id, frozenset(kwargs.items()))
        if cache_key in self._cache:
            return self._cache[cache_key].copy()
        try:
            value = self._resolve_impl(name, id, **kwargs)
        except KeyError:
            if name.startswith(("+", "-")):
                value = {
                    "composition": None,
                    "mass": float(name),
                    "name": name,
                    "id": None,
                    "provider": self.name,
                    "source": self,
                }
            else:
                raise
        self._cache[cache_key] = value
        return value.copy()

    def __call__(self, name=None, id=None, **kwargs):
        return self.resolve(name, id, **kwargs)

    def __eq__(self, other):
        return self.name == other.name

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(self.name)



class UnimodResolver(ModificationResolver):
    def __init__(self, **kwargs):
        super(UnimodResolver, self).__init__("unimod", **kwargs)
        self._database = kwargs.get("database")
        self.strict = kwargs.get("strict", True)

    def load_database(self):
        if _has_psims:
            return obo_cache.resolve("http://www.unimod.org/obo/unimod.obo")
        return Unimod()

    def _resolve_impl(self, name=None, id=None, **kwargs):
        strict = kwargs.get("strict", self.strict)
        exhaustive = kwargs.get("exhaustive", True)
        if name is not None:
            defn = self.database.by_title(name, strict=strict)
            if not defn:
                defn = self.database.by_name(name, strict=strict)
            if not defn and exhaustive and strict:
                defn = self.database.by_title(name, strict=False)
                if not defn:
                    defn = self.database.by_name(name, strict=False)
            if defn and isinstance(defn, list):
                warnings.warn(
                    "Multiple matches found for {!r} in Unimod, taking the first, {}.".format(
                        name, defn[0]['record_id']))
                defn = defn[0]
            if not defn:
                raise KeyError(name)
        elif id is not None:
            defn = self.database[id]
            if not defn:
                raise KeyError(id)
        else:
            raise ValueError("Must provide one of `name` or `id`")
        if isinstance(defn, dict):
            return {
                'composition': defn['composition'],
                'name': defn['title'],
                'id': defn['record_id'],
                'mass': defn['mono_mass'],
                'provider': self.name,
                "source": self
            }
        else:
            name = defn.ex_code_name
            if not name:
                name = defn.code_name
            return {
                "composition": defn.composition,
                "name": name,
                "id": defn.id,
                "mass": defn.monoisotopic_mass,
                "provider": self.name,
                "source": self
            }


class PSIModResolver(ModificationResolver):
    def __init__(self, **kwargs):
        super(PSIModResolver, self).__init__('psimod', **kwargs)
        self._database = kwargs.get("database")

    def load_database(self):
        return load_psimod()

    def _resolve_impl(self, name=None, id=None, **kwargs):
        if name is not None:
            defn = self.database[name]
        elif id is not None:
            defn = self.database['MOD:{:05d}'.format(id)]
        else:
            raise ValueError("Must provide one of `name` or `id`")

        # Non-standard amino acids are listed with `DiffMono` = `none`
        # but have a valid `MassMono` definition. Normally, `MassMono` is
        # the full mass of the residue plus the modification so it'd double count the
        # amino acid to use that value. Non-standard amino acids are a special case
        # because they *should* only be used with the amino acid X
        mass = None
        for key in ["DiffMono", "MassMono"]:
            if key in defn:
                try:
                    mass = float(defn[key])
                    break
                except (KeyError, TypeError, ValueError):
                    continue
        else:
            raise ModificationMassNotFoundError(
                "Could not resolve the mass of %r from %r" % ((name, id), defn)
            )

        # As with `DiffMono` for non-standard amino acids, but for chemical formulas -> Compositions
        for key in ["DiffFormula", "Formula"]:
            if key in defn and defn[key] is not None:
                composition = Composition()
                diff_formula_tokens = defn[key].strip().split(" ")
                for i in range(0, len(diff_formula_tokens), 2):
                    element = diff_formula_tokens[i]
                    count = diff_formula_tokens[i + 1]
                    if count:
                        count = int(count)
                    if element.startswith("("):
                        j = element.index(")")
                        isotope = element[1:j]
                        element = "%s[%s]" % (element[j + 1:], isotope)
                    composition[element] += count
                break
        else:
            composition = None
            warnings.warn("No formula was found for %r in PSI-MOD, composition will be missing" % ((name, id), ))
        return {
            'mass': mass,
            'composition': composition,
            'name': defn.name,
            'id': defn.id,
            'provider': self.name,
            "source": self
        }


class XLMODResolver(ModificationResolver):
    def __init__(self, **kwargs):
        super(XLMODResolver, self).__init__('xlmod', **kwargs)
        self._database = kwargs.get("database")

    def load_database(self):
        return load_xlmod()

    def _resolve_impl(self, name=None, id=None, **kwargs):
        if name is not None:
            defn = self.database[name]
        elif id is not None:
            defn = self.database['XLMOD:{:05d}'.format(id)]
        else:
            raise ValueError("Must provide one of `name` or `id`")
        try:
            mass = float(defn['monoIsotopicMass'])
        except (KeyError, TypeError, ValueError):
            raise ModificationMassNotFoundError("Could not resolve the mass of %r from %r" % ((name, id), defn))
        if 'deadEndFormula' in defn:
            composition = Composition(defn['deadEndFormula'].replace(" ", '').replace("D", "H[2]"))
        elif 'bridgeFormula' in defn:
            composition = Composition(
                defn['bridgeFormula'].replace(" ", '').replace("D", "H[2]"))
        return {
            'mass': mass,
            'composition': composition,
            'name': defn.name,
            'id': defn.id,
            'provider': self.name,
            "source": self
        }


# TODO: Implement resolve walking up the graph to get the mass. Can't really
# get any more information without glypy/glyspace interaction
class GNOResolver(ModificationResolver):
    mass_pattern = re.compile(r"(\d+(:?\.\d+)) Da")

    def __init__(self, **kwargs):
        super(GNOResolver, self).__init__('gnome', **kwargs)
        self._database = kwargs.get("database")

    def load_database(self):
        return load_gno()

    def get_mass_from_glycan_composition(self, term):
        '''Parse the Byonic-style glycan composition from property GNO:00000202
        to get the counts of each monosaccharide and use that to calculate mass.

        The mass computed here is exact and dehydrated, distinct from the rounded-off
        mass that :meth:`get_mass_from_term` will produce by walking up the CV term
        hierarchy. However, not all glycan compositions are representable in GNO:00000202
        format, so this may silently be absent or incomplete, hence the double-check in
        :meth:`get_mass_from_term`.

        Parameters
        ----------
        term : psims.controlled_vocabulary.Entity
            The CV entity being parsed.

        Returns
        -------
        mass : float or :const:`None`
            If a glycan composition is found on the term, the computed
            mass will be returned. Otherwise the :const:`None` is returned
        '''
        val = term.get('GNO:00000202')
        monosaccharides = BasicComposition()
        composition = Composition()
        if val:
            tokens = re.findall(r"([A-Za-z0-9]+)\((\d+)\)", val)
            mass = 0.0
            for symbol, count in tokens:
                count = int(count)
                try:
                    mono_mass, mono_comp, symbol = GlycanModification.valid_monosaccharides[symbol]
                    mass += mono_mass * count
                    composition += mono_comp * count
                    monosaccharides[symbol] += count
                except KeyError:
                    continue
            return mass, monosaccharides, composition
        return None, None, None

    def get_mass_from_term(self, term, raw_mass):
        '''Walk up the term hierarchy and find the mass group
        term near the root of the tree, and return the most accurate
        mass available for the provided term.

        The mass group term's mass is rounded to two decimal places, leading
        to relatively large errors.

        Parameters
        ----------
        term : psims.controlled_vocabulary.Entity
            The CV entity being parsed.

        Returns
        -------
        mass : float or :const:`None`
            If a root node is found along the term's lineage, computed
            mass will be returned. Otherwise the :const:`None` is returned.
            The mass may be
        '''
        root_id = 'GNO:00000001'
        parent = term.parent()
        if isinstance(parent, list):
            parent = parent[0]
        while parent.id != root_id:
            next_parent = parent.parent()
            if isinstance(next_parent, list):
                next_parent = next_parent[0]
            if next_parent.id == root_id:
                break
            parent = next_parent
        match = self.mass_pattern.search(parent.name)
        if not match:
            return None
        # This will have a small mass error.
        rough_mass = float(match.group(1)) - _WATER_MASS
        if raw_mass is not None and abs(rough_mass - raw_mass) < 1:
            return raw_mass
        warnings.warn(
            ("An accurate glycan composition could not be inferred from %s. "
             "Only a rough approximation is available.") % (term, ))
        return rough_mass

    def _resolve_impl(self, name=None, id=None, **kwargs):
        if name is not None:
            term = self.database[name]
        elif id is not None:
            term = self.database[id]
        else:
            raise ValueError("Must provide one of `name` or `id`")
        raw_mass, monosaccharides, composition = self.get_mass_from_glycan_composition(term)

        rec = {
            "name": term.name,
            "id": term.id,
            "provider": self.name,
            "composition": composition,
            "monosaccharides": monosaccharides,
            "mass": self.get_mass_from_term(term, raw_mass),
            "source": self
        }
        return rec


class GenericResolver(ModificationResolver):

    def __init__(self, resolvers, **kwargs):
        super(GenericResolver, self).__init__('generic', **kwargs)
        self.resolvers = list(resolvers)

    def load_database(self):
        return None

    def parse_identifier(self, identifier):
        """Parse a string that is either a CV prefixed identifier or name.

        Does no parsing as a :class:`GenericModification` is never qualified.

        Parameters
        ----------
        identifier : str
            The identifier string to parse, removing CV prefix as needed.

        Returns
        -------
        name : str, optional
            A textual identifier embedded in the qualified identifier, if any, otherwise
            :const:`None`.
        id : int, optional
            An integer ID embedded in the qualified identifier, if any, otherwise
            :const:`None`.
        """
        return identifier, None

    def _resolve_impl(self, name=None, id=None, **kwargs):
        defn = None
        for resolver in self.resolvers:
            try:
                defn = resolver(name=name, id=id, **kwargs)
                break
            except KeyError:
                continue
            except ModificationMassNotFoundError:
                warnings.warn("Could not resolve the mass for %r in %r" % ((name, id), resolver))
                continue
        if defn is None:
            if name is None:
                raise KeyError(id)
            elif id is None:
                raise KeyError(name)
            else:
                raise ValueError("Must provide one of `name` or `id`")
        return defn


class CustomResolver(ModificationResolver):
    store: Dict[str, Dict[str, Any]]

    def __init__(self, store: Dict[str, Dict[str, Any]]=None, **kwargs):
        if store is None:
            store = {}
        super().__init__("custom", **kwargs)
        self.store = store

    def _resolve_impl(self, name = None, id = None, **kwargs):
        if name is not None:
            return self.store[name]
        elif id is not None:
            return self.store[id]
        else:
            raise ValueError("Must provide one of `name` or `id`")

    def register(self, name, state: Dict[str, Any], **kwargs):
        state = state.copy()
        state.update(kwargs)
        state['id'] = name
        no_mass = "mass" not in state
        no_comp = "composition" not in state
        if no_mass and no_comp:
            raise ValueError("A custom modification definition *must* include at least one of `mass` or `composition`")
        self.store[name] = state


class ModificationBase(TagBase):
    '''A base class for all modification tags with marked prefixes.

    While :class:`ModificationBase` is hashable, its equality testing
    brings in additional tag-related information. For pure modification
    identity comparison, use :attr:`key` to get a :class:`ModificationToken`
    free of these concerns.
    '''

    _tag_type = None
    __slots__ = ('_definition', 'style', '_generated')

    _generated: ModificationSourceType

    def __init__(self, value, extra=None, group_id=None, style=None):
        if style is None:
            style = ModificationTagStyle.Unset
        super(ModificationBase, self).__init__(
            self._tag_type, value, extra, group_id)
        self._definition = None
        self._generated = ModificationSourceType.Explicit
        self.style = style

    def copy(self):
        return self.__class__(self.value, [e.copy() for e in self.extra], self.group_id, self.style)

    def __reduce__(self):
        return self.__class__, (self.value, self.extra, self.group_id, self.style), self.__getstate__()

    def __getstate__(self):
        if self._definition is None:
            return None
        state = self._definition.copy()
        state['source'] = None
        return state

    def __setstate__(self, state):
        self._definition = state

    def __eq__(self, other):
        if isinstance(other, ModificationToken):
            return other == self
        return super(ModificationBase, self).__eq__(other)

    def __hash__(self):
        return hash((self.id, self.provider))

    @property
    def key(self) -> 'ModificationToken':
        '''Get a safe-to-hash-and-compare :class:`ModificationToken`
        representing this modification without tag-like properties.

        Returns
        --------
        ModificationToken
        '''
        return ModificationToken(self.value, self.id, self.provider, self.__class__)

    @property
    def definition(self) -> Dict[str, Any]:
        '''A :class:`dict` of properties describing this modification, given
        by the providing controlled vocabulary. This value is cached, and
        should not be modified.

        Returns
        -------
        dict
        '''
        if self._definition is None:
            self._definition = self.resolve()
        return self._definition

    @property
    def mass(self):
        '''The monoisotopic mass shift this modification applies

        Returns
        -------
        float
        '''
        return self.definition['mass']

    def has_mass(self):
        """
        Check if this tag carries a mass value.

        Returns
        -------
        bool
        """
        return True

    def has_composition(self):
        return True

    @property
    def composition(self) -> Optional[Composition]:
        '''The chemical composition shift this modification applies'''
        return self.definition.get('composition')

    @property
    def charge(self) -> Optional[int]:
        return self.definition.get('charge')

    @property
    def id(self) -> Optional[int]:
        '''The unique identifier given to this modification by its provider

        Returns
        -------
        str or int
        '''
        return self.definition.get('id')

    @property
    def name(self):
        '''The primary name of this modification from its provider.

        Returns
        -------
        str
        '''
        return self.definition.get('name')

    @property
    def provider(self):
        '''The name of the controlled vocabulary that provided this
        modification.

        Returns
        -------
        str
        '''
        return self.definition.get('provider')

    def _populate_from_definition(self, definition):
        self._definition = definition

    def _format_main(self) -> str:
        if self.style == ModificationTagStyle.Unset or self.style is None:
            return "{self.prefix_name}:{self.value}".format(self=self)
        elif self.style == ModificationTagStyle.LongId:
            return "{self.prefix_name}:{self.id}".format(self=self)
        elif self.style == ModificationTagStyle.ShortId:
            return "{self.short_prefix}:{self.id}".format(self=self)
        elif self.style == ModificationTagStyle.LongName:
            return "{self.prefix_name}:{self.name}".format(self=self)
        elif self.style == ModificationTagStyle.ShortName:
            return "{self.short_prefix}:{self.name}".format(self=self)
        else:
            warnings.warn("Unknown formatting style {!r}".format(self.style))
            return "{self.prefix_name}:{self.value}".format(self=self)

    def resolve(self):
        '''Find the term and return it's properties
        '''
        keys = self.resolver.parse_identifier(self.value)
        return self.resolver(*keys)


class MassModification(TagBase):
    '''A modification defined purely by a signed mass shift in Daltons.

    The value of a :class:`MassModification` is always a :class:`float`
    '''
    __slots__ = ('_significant_figures', '_generated')

    prefix_name = "Obs"
    _generated: ModificationSourceType

    def __init__(self, value, extra=None, group_id=None):
        if isinstance(value, str):
            sigfigs = len(value.split('.')[-1].rstrip('0'))
        else:
            sigfigs = 4
        self._significant_figures = sigfigs
        self._generated = ModificationSourceType.Explicit
        super(MassModification, self).__init__(
            TagTypeEnum.massmod, float(value), extra, group_id)

    def copy(self):
        return self.__class__(self.value, [e.copy() for e in self.extra] if self.extra else [], self.group_id)

    def _format_main(self):
        if self.value >= 0:
            return ('+{0:0.{1}f}'.format(self.value, self._significant_figures)).rstrip('0').rstrip('.')
        else:
            return ('{0:0.{1}f}'.format(self.value, self._significant_figures)).rstrip('0').rstrip('.')

    @property
    def provider(self):
        return None

    @property
    def id(self):
        return self._format_main()

    @property
    def key(self) -> "ModificationToken":
        '''Get a safe-to-hash-and-compare :class:`ModificationToken`
        representing this modification without tag-like properties.

        Returns
        --------
        ModificationToken
        '''
        return ModificationToken(self.value, self.id, self.provider, self.__class__)

    @property
    def mass(self) -> float:
        return self.value

    def has_mass(self) -> bool:
        """
        Check if this tag carries a mass value.

        Returns
        -------
        bool
        """
        return True

    def has_composition(self) -> bool:
        return False

    def __eq__(self, other):
        if isinstance(other, ModificationToken):
            return other == self
        return super(MassModification, self).__eq__(other)

    def __hash__(self):
        return hash((self.id, self.provider))


class FormulaModification(ModificationBase):
    prefix_name = "Formula"

    charge_carrier_pattern: re.Pattern = re.compile(r':z((?:-|\+)?\d*)$')
    isotope_pattern: re.Pattern = re.compile(r'\[(?P<isotope>\d+)(?P<element>[A-Z][a-z]*)(?P<quantity>[\-+]?\d+)\]')
    _tag_type = TagTypeEnum.formula

    @staticmethod
    def _normalize_isotope_notation(match):
        '''Rewrite ProForma isotope notation to Pyteomics-compatible
        isotope notation.

        Parameters
        ----------
        match : Match
            The matched isotope notation string parsed by the regular expression.

        Returns
        reformatted : str
            The re-written isotope notation
        '''
        parts = match.groupdict()
        return "{element}[{isotope}]{quantity}".format(**parts)

    @classmethod
    def parse(cls, value: str):
        normalized = value.replace(" ", "")
        # If there is a [ character in the formula, we know there are isotopes which
        # need to be normalized.
        if "[" in normalized:
            normalized = cls.isotope_pattern.sub(
                cls._normalize_isotope_notation, normalized
            )
        if ":z" in normalized:
            matched = cls.charge_carrier_pattern.search(normalized)
            if not matched:
                raise ProFormaError(
                    "{normalized!r} is a malformed formula".format(
                        normalized=normalized
                    ),
                    None,
                    None,
                )
            charge = matched.group(1)
            charge = int(charge)
            normalized = cls.charge_carrier_pattern.sub("", normalized)
        else:
            charge = None
        composition = Composition(formula=normalized)
        return composition, charge

    def resolve(self):
        composition, charge = self.parse(self.value)
        return {
            "mass": composition.mass(),
            "composition": composition,
            "name": self.value,
            "charge": charge
        }


monosaccharide_description = namedtuple('monosaccharide_description', ('mass', 'composition', "symbol"))


class GlycanModification(ModificationBase):
    prefix_name = "Glycan"

    _tag_type = TagTypeEnum.glycan

    valid_monosaccharides = {
        "Hex": monosaccharide_description(162.0528, Composition("C6H10O5"), 'Hex'),
        "HexNAc": monosaccharide_description(203.0793, Composition("C8H13N1O5"), 'HexNAc'),
        "HexS": monosaccharide_description(242.009, Composition("C6H10O8S1"), 'HexS'),
        "HexP": monosaccharide_description(242.0191, Composition("C6H11O8P1"), 'HexP'),
        "HexNAcS": monosaccharide_description(283.0361, Composition("C8H13N1O8S1"), 'HexNAcS'),
        "dHex": monosaccharide_description(146.0579, Composition("C6H10O4"), 'dHex'),
        "NeuAc": monosaccharide_description(291.0954, Composition("C11H17N1O8"), 'NeuAc'),
        "NeuGc": monosaccharide_description(307.0903, Composition("C11H17N1O9"), 'NeuGc'),
        "Pen": monosaccharide_description(132.0422, Composition("C5H8O4"), 'Pen'),
        "Fuc": monosaccharide_description(146.0579, Composition("C6H10O4"), 'Fuc')
    }

    valid_monosaccharides['Neu5Ac'] = valid_monosaccharides['NeuAc']
    valid_monosaccharides['Neu5Gc'] = valid_monosaccharides['NeuGc']
    valid_monosaccharides['Pent'] = valid_monosaccharides['Pen']
    valid_monosaccharides['d-Hex'] = valid_monosaccharides['dHex']

    monomer_tokenizer = re.compile(
        r"|".join(sorted(valid_monosaccharides.keys(), key=len, reverse=True)))
    tokenizer = re.compile(r"(%s|[A-Za-z]+)\s*(\d*)\s*" % monomer_tokenizer.pattern)

    @property
    def monosaccharides(self):
        return self.definition.get('monosaccharides')

    def resolve(self):
        composite = BasicComposition()
        for tok, cnt in self.tokenizer.findall(self.value):
            if cnt:
                cnt = int(cnt)
            else:
                cnt = 1
            if tok not in self.valid_monosaccharides:
                parts = self.monomer_tokenizer.findall(tok)
                t = 0
                for p in parts:
                    if p not in self.valid_monosaccharides:
                        break
                    t += len(p)
                if t != len(tok):
                    raise ValueError("{tok!r} is not a valid monosaccharide name".format(tok=tok))
                else:
                    for p in parts[:-1]:
                        sym = self.valid_monosaccharides[p].symbol
                        composite[sym] += 1
                    sym = self.valid_monosaccharides[parts[-1]].symbol
                    composite[sym] += cnt
            else:
                sym = self.valid_monosaccharides[tok].symbol
                composite[sym] += cnt
        mass = 0
        chemcomp = Composition()
        for key, cnt in composite.items():
            try:
                m, c, sym = self.valid_monosaccharides[key]
            except KeyError:
                raise UnknownMonosaccharideError(key)
            mass += m * cnt
            chemcomp += c * cnt
        return {
            "mass": mass,
            "composition": chemcomp,
            "name": self.value,
            "monosaccharides": composite
        }


class UnimodModification(ModificationBase):
    __slots__ = ()

    resolver = UnimodResolver()

    prefix_name = "UNIMOD"
    short_prefix = "U"
    _tag_type = TagTypeEnum.unimod


class PSIModModification(ModificationBase):
    __slots__ = ()

    resolver = PSIModResolver()

    prefix_name = "MOD"
    short_prefix = 'M'
    _tag_type = TagTypeEnum.psimod


class GNOmeModification(ModificationBase):
    __slots__ = ()

    resolver = GNOResolver()

    prefix_name = "GNO"
    short_prefix = 'G'
    _tag_type = TagTypeEnum.gnome

    @property
    def monosaccharides(self):
        return self.definition.get('monosaccharides')


class XLMODModification(ModificationBase):
    __slots__ = ()

    resolver = XLMODResolver()

    prefix_name = "XLMOD"
    # short_prefix = 'XL'
    _tag_type = TagTypeEnum.xlmod


class CustomModification(ModificationBase):
    __slots__ = ()

    resolver = CustomResolver()

    prefix_name = 'Custom'
    short_prefix = 'C'

    @classmethod
    def register(cls, name, state: Dict[str, Any], **kwargs):
        return cls.resolver.register(name, state, **kwargs)


class GenericModification(ModificationBase):
    __slots__ = ()
    _tag_type = TagTypeEnum.generic
    resolver = GenericResolver([
        # Do exact matching here first. Then default to non-strict matching as a final
        # correction effort.
        partial(UnimodModification.resolver, exhaustive=False),
        PSIModModification.resolver,
        XLMODModification.resolver,
        GNOmeModification.resolver,
        # Some really common names aren't actually found in the XML exactly, so default
        # to non-strict matching now to avoid masking other sources here.
        partial(UnimodModification.resolver, strict=False)
    ])

    def __init__(self, value, extra=None, group_id=None, style=None):
        super(GenericModification, self).__init__(
            value, extra, group_id, style)

    def _format_main(self):
        return self.value

    def resolve(self):
        '''Find the term, searching through all available vocabularies and
        return the first match's properties
        '''
        keys = self.resolver.parse_identifier(self.value)
        defn = self.resolver(*keys)
        if defn is not None:
            return defn
        raise KeyError(keys)


def set_unimod_path(path):
    '''Set the path to load the Unimod database from for resolving
    ProForma Unimod modifications.

    .. note::

        This method ensures that the Unimod modification database loads
        quickly from a local database file instead of downloading a new
        copy from the internet.

    Parameters
    ----------
    path : str or file-like object
        A path to or file-like object for the "unimod.xml" file.

    Returns
    -------
    :class:`~pyteomics.mass.mass.Unimod`
    '''
    db = Unimod(path)
    UnimodModification.resolver.database = db
    return db


class ModificationToken(object):
    '''Describes a particular modification from a particular provider, independent
    of a :class:`TagBase`'s state.

    This class is meant to be used in place of a :class:`ModificationBase` object
    when equality testing and hashing is desired, but do not want extra properties
    to be involved.

    :class:`ModificationToken` is comparable and hashable, and can be compared with
    :class:`ModificationBase` subclass instances safely. It can be called to create
    a new instance of the :class:`ModificationBase` it is equal to.

    Attributes
    ----------
    name : str
        The name of the modification being represented, as the user specified it.
    id : int or str
        Whatever unique identifier the providing controlled vocabulary gave to this
        modification
    provider : str
        The name of the providing controlled vocabulary.
    source_cls : type
        A sub-class of :class:`ModificationBase` that will be used to fulfill this
        token if requested, providing it a resolver.
    '''
    __slots__ = ('name', 'id', 'provider', 'source_cls')

    name: str
    id: int
    provider: Callable
    source_cls: Union[Type[ModificationBase], Type[MassModification], Type['ModificationToken']]

    def __init__(self, name: str, id: int, provider: Callable, source_cls: Type):
        self.name = name
        self.id = id
        self.provider = provider
        self.source_cls = source_cls

    def __eq__(self, other):
        if other is None:
            return False
        if isinstance(other, (ModificationToken, ModificationBase, MassModification)):
            return self.id == other.id and self.provider == other.provider
        return False

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash((self.id, self.provider))

    def __call__(self):
        '''Create a new :class:`ModificationBase`
        instance from the provided :attr:`name`
        against :attr:`source_cls`'s resolver.

        Returns
        -------
        ModificationBase
        '''
        return self.source_cls(self.name)

    def __repr__(self):
        template = "{self.__class__.__name__}({self.name!r}, {self.id!r}, {self.provider!r}, {self.source_cls})"
        return template.format(self=self)


def split_tags(tokens: List[str]) -> List[List[str]]:
    '''Split a token array into discrete sets of tag
    tokens.

    Parameters
    ----------
    tokens: list
        The characters of the tag token buffer

    Returns
    -------
    list of list:
        The tokens for each contained tag
    '''
    starts = [0]
    ends = []
    for i, c in enumerate(tokens):
        if c == '|':
            ends.append(i)
            starts.append(i + 1)
        elif (i != 0 and c == '#'):
            ends.append(i)
            starts.append(i)
    ends.append(len(tokens))
    out = []
    for i, start in enumerate(starts):
        end = ends[i]
        tag = tokens[start:end]
        if len(tag) == 0:
            continue
        # Short circuit on INFO tags which can't be broken
        # if (tag[0] == 'i' and tag[:5] == ['i', 'n', 'f', 'o', ':']) or (tag[0] == 'I' and tag[:5] == ['I', 'N', 'F', 'O', ':']):
        #     tag = tokens[start:]
        #     out.append(tag)
        #     break
        out.append(tag)
    return out


def find_prefix(tokens: List[str]) -> Tuple[str, str]:
    '''Find the prefix, if any of the tag defined by `tokens`
    delimited by ":".

    Parameters
    ----------
    tokens: list
        The tag tokens to search

    Returns
    -------
    prefix: str or None
        The prefix string, if found
    rest: str
        The rest of the tokens, merged as a string
    '''
    for i, c in enumerate(tokens):
        if c == ':':
            return ''.join(tokens[:i]), ''.join(tokens[i + 1:])
    return None, ''.join(tokens)


def process_marker(tokens: Sequence[str]) -> Union[PositionLabelTag, LocalizationMarker]:
    '''Process a marker, which is a tag whose value starts with #.

    Parameters
    ----------
    tokens: list or str
        The tag tokens to parse

    Returns
    -------
    PositionLabelTag or LocalizationMarker
    '''
    if tokens[1:3] == 'XL':
        return PositionLabelTag(None, group_id=''.join(tokens))
    else:
        group_id = None
        value = None
        for i, c in enumerate(tokens):
            if c == '(':
                group_id = ''.join(tokens[:i])
                if tokens[-1] != ')':
                    raise Exception(
                        "Localization marker with score missing closing parenthesis")
                value = float(''.join(tokens[i + 1:-1]))
                return LocalizationMarker(value, group_id=group_id)
        else:
            group_id = ''.join(tokens)
            return PositionLabelTag(group_id=group_id)


def process_tag_tokens(tokens: List[str]) -> TagBase:
    '''Convert a tag token buffer into a parsed :class:`TagBase` instance
    of the appropriate sub-type with zero or more sub-tags.

    Parameters
    ----------
    tokens: list
        The tokens to parse

    Returns
    -------
    TagBase:
        The parsed tag
    '''
    parts = split_tags(tokens)
    main_tag = parts[0]
    if main_tag[0] in ('+', '-'):
        main_tag = ''.join(main_tag)
        main_tag = MassModification(main_tag)
    elif main_tag[0] == '#':
        main_tag = process_marker(main_tag)
    else:
        prefix, value = find_prefix(main_tag)
        if prefix is None:
            value = ''.join(value)
            if value.lower() in TagBase.prefix_map:
                main_tag = TagBase.prefix_map[value.lower()]()
            else:
                main_tag = GenericModification(''.join(value))
        else:
            try:
                tag_type = TagBase.find_by_tag(prefix)
                main_tag = tag_type(value)
            except KeyError:
                main_tag_str = ''.join(main_tag)
                main_tag = GenericModification(main_tag_str)

    if len(parts) > 1:
        extras = []
        for part in parts[1:]:
            prefix, value = find_prefix(part)
            if prefix is None:
                if value[0] == "#":
                    marker = process_marker(value)
                    if isinstance(marker, PositionLabelTag):
                        main_tag.group_id = ''.join(value)
                    else:
                        main_tag.group_id = marker.group_id
                        extras.append(marker)
                else:
                    value = ''.join(value)
                    if value.lower() in TagBase.prefix_map:
                        extra_tag = TagBase.prefix_map[value.lower()]()
                    else:
                        extra_tag = GenericModification("".join(value))
                    extras.append(extra_tag)
            else:
                try:
                    tag_type = TagBase.find_by_tag(prefix)
                    extra_tag = tag_type(value)
                except KeyError:
                    part_str = ''.join(part)
                    extra_tag = GenericModification(part_str)
                extras.append(extra_tag)
        main_tag.extra = extras
    return main_tag


class ModificationTarget(object):
    aa: str
    n_term: bool
    c_term: bool

    def __init__(self, aa, n_term=False, c_term=False):
        self.aa = aa
        self.n_term = n_term
        self.c_term = c_term

    def __eq__(self, other):
        if isinstance(other, str):
            return str(self) == other
        else:
            return (
                self.aa == other.aa
                and self.n_term == other.n_term
                and self.c_term == other.c_term
            )

    def __ne__(self, other):
        if isinstance(other, str):
            return str(self) != other
        else:
            return (
                self.aa != other.aa
                or self.n_term != other.n_term
                or self.c_term != other.c_term
            )

    def __hash__(self):
        return hash(str(self))

    def __str__(self):
        buffer = []
        if self.n_term:
            buffer.append('N-term')
        if self.c_term:
            buffer.append('C-term')
        if self.aa:
            buffer.append(self.aa)
        return ':'.join(buffer)

    def __repr__(self):
        return str(self)

    def is_valid(self, aa: str, n_term: bool, c_term: bool) -> bool:
        if (n_term and self.n_term) or (c_term and self.c_term):
            if (self.aa and aa.upper() == self.aa.upper()) or self.aa is None:
                return True
            return False
        return self.aa.upper() == aa.upper() or self.aa is None


class ModificationRule(object):
    '''Define a fixed modification rule which dictates a modification tag is
    always applied at one or more amino acid residues.

    Attributes
    ----------
    modification_tag: TagBase
        The modification to apply
    targets: list
        The list of amino acids this applies to
    '''
    __slots__ = ('modification_tag', 'targets')

    modification_tag: TagBase
    targets: List[ModificationTarget]

    def __init__(self, modification_tag: TagBase, targets: Union[ModificationTarget, List[ModificationTarget], None]=None):
        self.modification_tag = modification_tag
        self.targets = targets
        self._validate_targets()

    def is_valid(self, aa: str, n_term: bool, c_term: bool) -> bool:
        return any(target.is_valid(aa, n_term, c_term) for target in self.targets)

    def _validate_targets(self):
        validated_targets = []
        if self.targets is None:
            self.targets = []
        elif not isinstance(self.targets, list):
            self.targets = [self.targets]
        for target in self.targets:
            if isinstance(target, ModificationTarget):
                validated_targets.append(target)
            elif target in VALID_AA:
                validated_targets.append(ModificationTarget(target, False, False))
            elif target in ("N-term", "C-term"):
                n_term = target == "N-term"
                c_term = target == "C-term"
                validated_targets.append(ModificationTarget(None, n_term, c_term))
            elif target.startswith(("N-term:", "C-term:")):
                tokens = target.split(":")
                if len(tokens) == 2:
                    if tokens[1] in VALID_AA:
                        n_term = tokens[0] == "N-term"
                        c_term = tokens[0] == "C-term"
                        validated_targets.append(ModificationTarget(tokens[1], n_term, c_term))
                    else:
                        raise PyteomicsError(
                            "Modification rule {0} has an invalid amino acid specific terminal target {2} in {1}".format(
                                self,
                                target,
                                tokens[1]
                            )
                        )
                else:
                    raise PyteomicsError(
                        "Modification rule {0} has an empty amino acid specific terminal target {1}".format(
                            self, target
                        )
                    )
            else:
                raise PyteomicsError(
                    "Modification rule {0} has an invalid target {1}".format(
                        self, target
                    )
                )

        self.targets = validated_targets

    def __eq__(self, other):
        if other is None:
            return False
        return self.modification_tag == other.modification_tag and self.targets == other.targets

    def __ne__(self, other):
        return not self == other

    def __str__(self):
        targets = ','.join(map(str, self.targets))
        return "<[{self.modification_tag}]@{targets}>".format(self=self, targets=targets)

    def __repr__(self):
        return "{self.__class__.__name__}({self.modification_tag!r}, {self.targets})".format(self=self)


class StableIsotope(object):
    '''
    Define a fixed isotope that is applied globally to all amino acids.

    Attributes
    ----------
    isotope: str
        The stable isotope string, of the form [<isotope-number>]<element> or a special
        isotopoform's name.
    '''
    __slots__ = ('isotope', )
    isotope: str

    def __init__(self, isotope):
        self.isotope = isotope

    def copy(self):
        return self.__class__(self.isotope)

    def __eq__(self, other):
        if other is None:
            return False
        return self.isotope == other.isotope

    def __ne__(self, other):
        return not self == other

    def __str__(self):
        return "<{self.isotope}>".format(self=self)

    def __repr__(self):
        return "{self.__class__.__name__}({self.isotope})".format(self=self)


class IntersectionEnum(Enum):
    no_overlap = 0
    full_contains_interval = 1
    full_contained_in_interval = 2
    start_overlap = 3
    end_overlap = 4


class TaggedInterval(object):
    '''Define a fixed interval over the associated sequence which contains the localization
    of the associated tag or denotes a region of general sequence order ambiguity.

    Attributes
    ----------
    start: int
        The starting position (inclusive) of the interval along the primary sequence
    end: int
        The ending position (exclusive) of the interval along the primary sequence
    tags: list[TagBase]
        The tags being localized
    ambiguous : bool
        Whether the interval is ambiguous or not
    '''
    __slots__ = ('start', 'end', 'tags', 'ambiguous')

    start: int
    end: Optional[int]
    tags: Optional[List[TagBase]]
    ambiguous: bool

    def __init__(self, start, end=None, tags=None, ambiguous=False):
        self.start = start
        self.end = end
        self.tags = tags
        self.ambiguous = ambiguous

    def copy(self):
        return self.__class__(
            self.start,
            self.end,
            [v.copy() for v in self.tags] if self.tags else [],
            self.ambiguous
        )

    def __eq__(self, other):
        if other is None:
            return False
        return self.start == other.start and self.end == other.end and self.tags == other.tags

    def __ne__(self, other):
        return not self == other

    def __str__(self):
        return f"({'?' if self.ambiguous else ''}{self.start}-{self.end}){self.tags!r}"

    def __repr__(self):
        return f"{self.__class__.__name__}({self.start}, {self.end}, {self.tags}, ambiguous={self.ambiguous})"

    def as_slice(self):
        return slice(self.start, self.end)

    def contains(self, i):
        return self.start <= i < self.end

    def __contains__(self, i):
        return self.contains(i)

    def _check_slice(self, qstart, qend, warn_ambiguous):
        # Fully contained interval
        valid = qstart <= self.start and qend >= self.end
        case = IntersectionEnum.full_contained_in_interval if valid else IntersectionEnum.no_overlap
        if not valid:
            # Spans the beginning but not the end
            valid = qstart <= self.start and qend > self.start
            if valid:
                case = IntersectionEnum.start_overlap
                if warn_ambiguous:
                    warnings.warn("Slice bisecting interval %s" % (self, ))

        if not valid:
            # Spans the end but not the beginning
            valid = qstart < self.end and qend > self.end
            if valid:
                case = IntersectionEnum.end_overlap
                if warn_ambiguous:
                    warnings.warn("Slice bisecting interval %s" % (self, ))

        if not valid:
            # Contained interval
            valid = qstart >= self.start and qend < self.end
            if valid:
                case = IntersectionEnum.full_contains_interval
                if warn_ambiguous:
                    warnings.warn("Slice bisecting interval %s" % (self, ))
        return valid, case

    def _update_coordinates_sliced(self, start=None, end=None, warn_ambiguous=True):
        if end is None:
            qend = self.end + 1
        else:
            qend = end
        if start is None:
            qstart = self.start - 1
        else:
            qstart = start

        valid, intersection_type = self._check_slice(qstart, qend, warn_ambiguous)
        if self.ambiguous and intersection_type not in (IntersectionEnum.full_contained_in_interval, IntersectionEnum.no_overlap):
            raise ValueError("Cannot bisect an ambiguous interval")
        if not valid:
            return None
        new = self.copy()
        if start is not None:
            diff = self.start - start
            if diff < 0:
                diff = 0
            new.start = diff
        if end is not None:
            width = min(new.end, end) - self.start
        else:
            width = self.end - max(start, self.start)
        new.end = new.start + width
        return new


class Adduct(NamedTuple):
    name: str
    charge: int
    count: int

    def __str__(self):
        base = f"{self.name}:z{'+' if self.charge > 0 else ''}{self.charge}"
        if self.count > 1:
            return base + f"^{self.count}"
        return base

    def composition(self) -> Composition:
        if self.name == 'e-':
            return Composition({"e-": self.count})
        return Composition(formula=self.name) * self.count

    def mass(self) -> float:
        return Composition(formula=self.name).mass * self.count

    def total_charge(self) -> int:
        return self.charge * self.count


class ChargeState(object):
    """Describes the charge and adduct types of the structure.

    Attributes
    ----------
    charge : int
        The total charge state as a signed number.
    adducts : list[Adduct]
        Each charge carrier associated with the molecule.
    """
    __slots__ = ("charge", "adducts")

    charge: int
    adducts: List[Adduct]

    @classmethod
    def from_adducts(cls, adducts: List[Adduct]):
        acc = 0
        for a in adducts:
            acc += a.charge * a.count
        return cls(acc, adducts)

    def __init__(self, charge, adducts=None):
        if adducts is None:
            adducts = [Adduct("H", 1, charge)]
        self.charge = charge
        self.adducts = adducts

    def for_mz_calculation(self) -> Tuple[float, int]:
        """
        Get the total mass of the charge carrier(s) and their collective charge
        to plug into the formula for mass-to-charge-ratio, ``(mass of molecule + mass of charge carrier) / charge``

        Returns
        -------
        charge_carrier_mass : float
            The total mass of the charge carriers(s) in the adducting group(s)
        charge : int
            The total charge contributed by all the charge carriers
            in the adducting group(s)
        """
        mass = 0.0
        for a in self.adducts:
            mass += a.mass()
        return (mass, self.charge)

    def composition(self):
        comp = Composition()
        for a in self.adducts:
            comp += a.composition()
        return comp

    def __str__(self):
        if len(self.adducts) > 1 or self.adducts[0].name != 'H':
            tokens = []
            tokens.append("[")
            tokens.append(','.join((map(str, self.adducts))))
            tokens.append("]")
            return ''.join(tokens)
        else:
            return f'{self.charge:d}'

    def __repr__(self):
        template = "{self.__class__.__name__}({self.charge}, {self.adducts})"
        return template.format(self=self)


class TokenBuffer(Generic[T]):
    '''A token buffer that wraps the accumulation and reset logic
    of a list of :class:`str` objects.

    Implements a subset of the Sequence protocol.

    Attributes
    ----------
    buffer: list
        The list of tokens accumulated since the last parsing.
    '''
    buffer: List[str]
    boundaries: List[int]

    def __init__(self, initial=None):
        self.buffer = list(initial or [])
        self.boundaries = []

    def append(self, c: str):
        '''
        Append a new character to the buffer.

        Parameters
        ----------
        c: str
            The character appended
        '''
        self.buffer.append(c)

    def extend(self, cs: str):
        '''
        Extend the buffer with additional characters

        Parameters
        ----------
        cs: str
            The chracters to append
        '''
        self.buffer.extend(cs)

    def reset(self):
        '''Discard the content of the current buffer.
        '''
        if self.buffer:
            self.buffer = []
        if self.boundaries:
            self.boundaries = []

    def __bool__(self):
        return bool(self.buffer)

    def __iter__(self):
        return iter(self.buffer)

    def __getitem__(self, i):
        return self.buffer[i]

    def __len__(self):
        return len(self.buffer)

    def tokenize(self) -> List[str]:
        i = 0
        pieces = []
        for k in self.boundaries + [len(self)]:
            piece = self.buffer[i:k]
            i = k
            pieces.append(piece)
        return pieces

    def _transform(self, value: T) -> T:
        return value

    def process(self) -> Union[T, List[T]]:
        if self.boundaries:
            value = [self._transform(v) for v in self.tokenize()]
        else:
            value = self._transform(self.buffer)
        self.reset()
        return value

    def bound(self) -> int:
        k = len(self)
        self.boundaries.append(k)
        return k

    def __call__(self) -> Union[T, List[T]]:
        return self.process()


class NumberParser(TokenBuffer[int]):
    '''A buffer which accumulates tokens until it is asked to parse them into
    :class:`int` instances.
    '''

    def _transform(self, value) -> int:
        return int(''.join(value))


class StringParser(TokenBuffer[str]):
    '''A buffer which accumulates tokens until it is asked to parse them into
    :class:`str` instances.
    '''

    def _transform(self, value) -> str:
        return ''.join(value)


class AdductParser(StringParser):
    '''A buffer which accumulates tokens related to adducts until it is asked to parse them into
    a list of [(str, int)] tuples, where the first element is the adduct name
    and the second element is the number of adducts of that type.
    '''
    token_pattern = re.compile(r'(?P<number>[+-]?\d*)(?P<adduct>[A-Za-z]+)(?P<charge>\d*[+-])')
    token_pattern2 = re.compile(r'(?P<adduct>[A-Za-z]+):[zZ](?P<charge>(-|\+)\d+)(?:\^(?P<number>\d+))?')

    def parse_form1(self, token: str) -> Optional[Tuple[str, int, int]]:
        parsed = self.token_pattern.match(token)
        if not parsed:
            return None
        gdict = parsed.groupdict()
        if gdict['adduct'] == 'e':
            adduct = 'e-'
        else:
            adduct = gdict['adduct']
        if gdict['number'] == '+' or gdict['number'] == '':
            number = 1
        elif gdict['number'] == '-':
            number = -1
        else:
            number = int(gdict['number'])
        charge = int(gdict['charge'][:-1]) if gdict['charge'][:-1] else 1
        if gdict['charge'][-1] == '-':
            charge = -charge
        return (adduct, charge, number)

    def parse_form2(self, token: str):
        parsed = self.token_pattern2.match(token)
        if not parsed:
            return None
        gdict = parsed.groupdict()
        if gdict['adduct'] == 'e':
            adduct = 'e-'
        else:
            adduct = gdict['adduct']
        if gdict['number'] == '+' or gdict['number'] == '':
            number = 1
        elif gdict['number'] == '-':
            number = -1
        elif gdict['number'] is not None:
            number = int(gdict['number'])
        else:
            number = 1
        charge = int(gdict['charge']) if gdict['charge'] else 1
        return (adduct, charge, number)

    def process(self):
        value = []
        for token in self.tokenize():
            if not isinstance(token, str):
                token = ''.join(token)
            try:
                adduct = self.parse_form1(token)
                if not adduct:
                    adduct = self.parse_form2(token)
                if not adduct:
                    raise ProFormaError(
                        "Invalid adduct token {!r} in {!r}".format(token, self.buffer)
                    )
                value.append(Adduct(*adduct))
            except AttributeError:
                raise ProFormaError("Invalid adduct token {!r} in {!r}".format(token, self.buffer))
        return value


class TagParser(TokenBuffer[TagBase]):
    '''A buffer which accumulates tokens until it is asked to parse them into
    :class:`TagBase` instances.

    Implements a subset of the Sequence protocol.

    Attributes
    ----------
    buffer: list
        The list of tokens accumulated since the last parsing.
    group_ids: set
        The set of all group IDs that have been produced so far.
    '''

    def __init__(self, initial=None, group_ids=None):
        super(TagParser, self).__init__(initial)
        if group_ids:
            self.group_ids = set(group_ids)
        else:
            self.group_ids = set()

    def _transform(self, value):
        tag = process_tag_tokens(value)
        if tag.group_id:
            self.group_ids.add(tag.group_id)
        return tag

    def process(self):
        value = super(TagParser, self).process()
        if not isinstance(value, list):
            value = [value]
        return value


class ParserStateEnum(Enum):
    before_sequence = 0
    tag_before_sequence = 1
    global_tag = 2
    fixed_spec = 3
    labile_tag = 4
    sequence = 5
    tag_in_sequence = 6
    interval_tag = 7
    tag_after_sequence = 8
    stable_isotope = 9
    post_tag_before = 10
    unlocalized_count = 11
    post_global = 12
    post_global_aa = 13
    post_interval_tag = 14
    post_tag_after = 15
    charge_state_start = 16
    charge_state_number = 17
    charge_state_adduct_start = 18
    charge_state_adduct_end = 19
    inter_chain_cross_link_start = 20
    chimeric_start = 21
    interval_initial = 22
    post_global_terminal = 23

    peptidoform_name_start = 24
    peptidoform_name_level = 25
    peptidoform_name_text = 26
    peptidoform_name_close = 27

    done = 999


BEFORE = ParserStateEnum.before_sequence
TAG_BEFORE = ParserStateEnum.tag_before_sequence
FIXED = ParserStateEnum.fixed_spec
GLOBAL = ParserStateEnum.global_tag
ISOTOPE = ParserStateEnum.stable_isotope
LABILE = ParserStateEnum.labile_tag
SEQ = ParserStateEnum.sequence
TAG = ParserStateEnum.tag_in_sequence
INTERVAL_TAG = ParserStateEnum.interval_tag
INTERVAL_INIT = ParserStateEnum.interval_initial
TAG_AFTER = ParserStateEnum.tag_after_sequence
POST_TAG_BEFORE = ParserStateEnum.post_tag_before
POST_TAG_AFTER = ParserStateEnum.post_tag_after
UNLOCALIZED_COUNT = ParserStateEnum.unlocalized_count

POST_GLOBAL = ParserStateEnum.post_global
POST_GLOBAL_AA = ParserStateEnum.post_global_aa
POST_GLOBAL_TERM = ParserStateEnum.post_global_terminal
POST_INTERVAL_TAG = ParserStateEnum.post_interval_tag

CHARGE_START = ParserStateEnum.charge_state_start
CHARGE_NUMBER = ParserStateEnum.charge_state_number

ADDUCT_START = ParserStateEnum.charge_state_adduct_start
ADDUCT_END = ParserStateEnum.charge_state_adduct_end

PEPTIDOFORM_NAME_START = ParserStateEnum.peptidoform_name_start
PEPTIDOFORM_NAME_LEVEL = ParserStateEnum.peptidoform_name_level
PEPTIDOFORM_NAME_TEXT = ParserStateEnum.peptidoform_name_text
PEPTIDOFORM_NAME_CLOSE = ParserStateEnum.peptidoform_name_close

DONE = ParserStateEnum.done

VALID_AA_UPPER = set("QWERTYIPASDFGHKLCVNMXUOJZB")
VALID_AA = {s.lower() for s in VALID_AA_UPPER} | VALID_AA_UPPER
TERMINAL_SPEC_CHARS = set('N-term') | set('C-term') | set("ncT: ")


class Parser:
    """
    A parser for the ProForma 2 syntax.

    Attributes
    ----------
    sequence : str
        The sequence to be parsed
    index : int
        The current index parsing from
    depth : int
        The current depth of the brace type being parsed within
    length : int
        The total length in characters of the sequence to parse
    state : ParserStateEnum
        The state of the parser is currently in, dictating how it will interpret
        the next token read.
    """
    sequence: str
    index: int
    depth: int
    length: int
    state: ParserStateEnum

    labile_modifications: List[TagBase]
    fixed_modifications: List[TagBase]
    unlocalized_modifications: List[TagBase]
    intervals: List[TaggedInterval]
    isotopes: List[StableIsotope]
    n_term: List[TagBase]
    c_term: List[TagBase]
    positions: List
    current_aa: str
    current_interval: Optional[TaggedInterval]
    current_tag: TagParser
    current_unlocalized_count: NumberParser
    current_aa_targets: StringParser
    charge_buffer: Optional[NumberParser]
    adduct_buffer: Optional[AdductParser]

    def __init__(self, sequence: str, case_sensitive_aa: bool=False):
        """
        Instantiate a ProForma 2 parser for the specified sequence.

        Parameters
        ----------
        sequence : str
            The sequence to parse
        case_sensitive_aa : bool
            Whether to treat amino acids as case sensitive (older behavior) while the specification
            states they should be handled insensitively.
        """
        self.sequence = sequence
        self.index = 0
        self.depth = 0
        self.length = len(sequence)
        self.state = ParserStateEnum.before_sequence
        self._VALID_AA = VALID_AA if not case_sensitive_aa else VALID_AA_UPPER

        self.n_term = []
        self.c_term = []
        self.intervals = []
        self.positions = []

        self.adduct_buffer = None
        self.charge_buffer = None
        self.current_aa = None
        self.current_interval = None
        self.current_tag = TagParser()
        self.current_aa_targets = StringParser()
        self.current_unlocalized_count = NumberParser()

        self.unlocalized_modifications = []
        self.labile_modifications = []
        self.fixed_modifications = []
        self.unlocalized_modifications = []
        self.intervals = []
        self.isotopes = []
        self.name_level = None
        self.name_buffer = StringParser()
        self.names = {}

    @property
    def i(self) -> int:
        return self.index

    @i.setter
    def i(self, value: int):
        self.index = value

    @property
    def n(self) -> int:
        return self.length

    def pack_sequence_position(self):
        self.positions.append(
            (
                self.current_aa,
                self.current_tag() if self.current_tag else None,
            )
        )
        self.current_aa = None

    def handle_before(self, c: str):
        if c == '[':
            self.state = TAG_BEFORE
            self.depth = 1
        elif c == '{':
            self.state = LABILE
            self.depth = 1
        elif c == '<':
            self.state = FIXED
        elif c in self._VALID_AA:
            self.current_aa = c
            self.state = SEQ
        elif c == '(':
            if (self.index + 1) < self.length:
                if self.sequence[self.index + 1] == '>':
                    self.index += 1
                    self.state = PEPTIDOFORM_NAME_LEVEL
                    self.depth = 1
                    self.name_level = 1
                else:
                    self.state = INTERVAL_INIT
                    self.current_interval = TaggedInterval(len(self.positions) + 1)
        else:
            raise ProFormaError(
                f"Error In State {self.state}, unexpected {c} found at index {self.i}",
                self.i,
                self.state,
            )

    def handle_seq(self, c: str):
        state = self.state
        if state == INTERVAL_INIT:
            self.state = SEQ
            if c == '?':
                if self.current_interval is not None:
                    self.current_interval.ambiguous = True
                # continue
                return True
        if c in self._VALID_AA:
            if self.current_aa is not None:
                self.pack_sequence_position()
            self.current_aa = c
        elif c == '[':
            self.state = TAG
            if self.current_tag:
                self.current_tag.bound()
            self.depth = 1
        elif c == '(':
            if self.current_interval is not None:
                raise ProFormaError(
                    (
                        f"Error In State {self.state}, nested range found at index {self.index}. "
                        "Nested ranges are not yet supported by ProForma."
                    ),
                    self.index,
                    self.state,
                )
            self.current_interval = TaggedInterval(len(self.positions) + 1)
            self.state = INTERVAL_INIT
        elif c == ')':
            self.pack_sequence_position()
            if self.current_interval is None:
                raise ProFormaError(
                    f"Error In State {self.state}, unexpected {c} found at index {self.index}",
                    self.index,
                    self.state,
                )
            else:
                self.current_interval.end = len(self.positions)
                if self.i + 1 < self.n and self.sequence[self.i + 1] == "[":
                    self.i += 1
                    self.depth = 1
                    self.state = INTERVAL_TAG
                else:
                    self.intervals.append(self.current_interval)
                    self.current_interval = None
        elif c == '-':
            if self.current_aa:
                self.pack_sequence_position()
            self.state = TAG_AFTER
            if self.i >= self.n or self.sequence[self.i + 1] != "[":
                raise ProFormaError("Missing Opening Tag", self.i, self.state)
            self.i += 1
            self.depth = 1
        elif c == '/':
            self.state = CHARGE_START
            self.charge_buffer = NumberParser()
            self.adduct_buffer = AdductParser()
        elif c == '+':
            raise ProFormaError(
                f"Error In State {self.state}, {c} found at index {self.i}. Chimeric representation not supported",
                self.i,
                self.state,
            )
        else:
            raise ProFormaError(
                f"Error In State {self.state}, unexpected {c} found at index {self.i}",
                self.i,
                self.state,
            )

    def handle_tag(self, c: str):
        if c == "[":
            self.depth += 1
            self.current_tag.append(c)
        elif c == "]":
            self.depth -= 1
            if self.depth <= 0:
                self.depth = 0
                if self.state == TAG:
                    self.state = SEQ
                elif self.state == TAG_BEFORE:
                    self.state = POST_TAG_BEFORE
                elif self.state == TAG_AFTER:
                    self.c_term = self.current_tag()
                    self.state = POST_TAG_AFTER
                elif self.state == GLOBAL:
                    self.state = POST_GLOBAL
                elif self.state == INTERVAL_TAG:
                    self.state = POST_INTERVAL_TAG
                    # self.current_interval.tags.append(self.current_tag())
                    self.depth = 0
            else:
                self.current_tag.append(c)
        else:
            self.current_tag.append(c)

    def handle_fixed(self, c: str):
        if c == '[':
            self.state = GLOBAL
        else:
            # Do validation here
            self.state = ISOTOPE
            self.current_tag.reset()
            self.current_tag.append(c)

    def handle_isotope(self, c: str):
        if c != ">":
            self.current_tag.append(c)
        else:
            # Not technically a tag, but exploits the current buffer
            self.isotopes.append(StableIsotope("".join(self.current_tag)))
            self.current_tag.reset()
            self.state = BEFORE

    def handle_labile(self, c: str):
        if c == "{":
            self.depth += 1
        elif c == "}":
            self.depth -= 1
            if self.depth <= 0:
                self.depth = 0
                self.labile_modifications.append(self.current_tag()[0])
                self.state = BEFORE
        else:
            self.current_tag.append(c)

    def handle_post_interval_tag(self, c: str):
        if c == "[":
            self.current_tag.bound()
            self.state = INTERVAL_TAG
        elif c in self._VALID_AA:
            self.current_aa = c
            self.current_interval.tags = self.current_tag()
            self.intervals.append(self.current_interval)
            self.current_interval = None
            self.state = SEQ
        elif c == "-":
            self.state = TAG_AFTER
            # Unroll next state to immediately fall into a tag parsing state instead of
            # including a separate post-dash state
            if self.i >= self.n or self.sequence[self.i] != "[":
                raise ProFormaError("Missing Closing Tag", self.i, self.state)
            self.i += 1
            self.depth = 1
        elif c == "/":
            self.state = CHARGE_START
            self.charge_buffer = NumberParser()
        elif c == "+":
            raise ProFormaError(
                f"Error In State {self.state}, {self.c} found at index {self.i}. Chimeric representation not supported",
                self.i,
                self.state,
            )
        else:
            raise ProFormaError(
                f"Error In State {self.state}, unexpected {self.c} found at index {self.i}",
                self.i,
                self.state,
            )

    def handle_post_tag_before(self, c: str):
        if c == "?":
            self.unlocalized_modifications.extend(self.current_tag())
            self.state = BEFORE
        elif c == "-":
            self.n_term = self.current_tag()
            self.state = BEFORE
        elif c == "^":
            self.state = UNLOCALIZED_COUNT
        elif c == "[":
            self.current_tag.bound()
            self.state = TAG_BEFORE
        else:
            raise ProFormaError(
                f"Error In State {self.state}, unexpected {self.c} found at index {self.i}",
                self.i,
                self.state,
            )

    def handle_unlocalized_count(self, c: str):
        if c.isdigit():
            self.current_unlocalized_count.append(c)
        elif c == "[":
            self.state = TAG_BEFORE
            self.depth = 1
            tags = self.current_tag()
            tags, tag = tags[:-1], tags[-1]
            self.unlocalized_modifications.extend(tags)
            multiplicity = self.current_unlocalized_count()
            for _ in range(multiplicity):
                self.unlocalized_modifications.append(tag)
        elif c == "?":
            self.state = BEFORE
            tags = self.current_tag()
            tags, tag = tags[:-1], tags[-1]
            self.unlocalized_modifications.extend(tags)
            multiplicity = self.current_unlocalized_count()
            for _ in range(multiplicity):
                self.unlocalized_modifications.append(tag)
        else:
            raise ProFormaError(
                f"Error In State {self.state}, unexpected {c} found at index {self.i}",
                self.i,
                self.state,
            )

    def handle_post_global(self, c: str):
        if c == "@":
            self.state = POST_GLOBAL_AA
        else:
            raise ProFormaError(
                (
                    f"Error In State {self.state}, fixed modification detected without "
                    f"target amino acids found at index {self.i}"
                ),
                self.i,
                self.state,
            )

    def handle_post_global_aa(self, c: str):
        if c in self._VALID_AA or c in TERMINAL_SPEC_CHARS:
            self.current_aa_targets.append(c)
        elif c == ",":
            # the next character should be another amino acid
            self.current_aa_targets.bound()
        elif c == ">":
            try:
                v = self.current_aa_targets()
                self.fixed_modifications.append(
                    ModificationRule(self.current_tag()[0], v)
                )
            except PyteomicsError as err:
                raise ProFormaError(
                    (
                        f"Error In State {self.state}, fixed modification detected invalid "
                        f"target found at index {self.i}: {err}"
                    ),
                    self.i,
                    self.state,
                )
            self.state = BEFORE
        else:
            raise ProFormaError(
                f"Error In State {self.state}, unclosed fixed modification rule",
                self.i,
                self.state,
            )

    def handle_post_tag_after(self, c: str):
        if c == "/":
            self.state = CHARGE_START
            self.charge_buffer = NumberParser()
        elif c == "+":
            raise ProFormaError(
                f"Error In State {self.state}, {c} found at index {self.i}. Chimeric representation not supported",
                self.i,
                self.state,
            )

    def handle_charge_start(self, c: str):
        if c in "+-":
            self.charge_buffer.append(c)
            self.state = CHARGE_NUMBER
        elif c.isdigit():
            self.charge_buffer.append(c)
            self.state = CHARGE_NUMBER
        elif c == "/":
            self.state = ParserStateEnum.inter_chain_cross_link_start
            raise ProFormaError(
                "Inter-chain cross-linked peptides are not yet supported",
                self.i,
                self.state,
            )
        elif c == '[':
            self.state = ParserStateEnum.charge_state_adduct_start
        else:
            raise ProFormaError(
                f"Error In State {self.state}, unexpected {c} found at index {self.i}",
                self.i,
                self.state,
            )

    def handle_charge_number(self, c: str):
        if c.isdigit():
            self.charge_buffer.append(c)
        elif c == "[":
            self.state = ADDUCT_START
            self.adduct_buffer = AdductParser()
        else:
            raise ProFormaError(
                f"Error In State {self.state}, unexpected {c} found at index {self.i}",
                self.i,
                self.state,
            )

    def handle_adduct_start(self, c: str):
        if c.isdigit() or c in "+:-" or c.isalpha():
            self.adduct_buffer.append(c)
        elif c == ",":
            self.adduct_buffer.bound()
        elif c == "]":
            self.state = ADDUCT_END

    def handle_adduct_end(self, c: str):
        if c == "+":
            raise ProFormaError(
                f"Error In State {self.state}, {c} found at index {self.i}. Chimeric representation not supported",
                self.i,
                self.state,
            )

    def handle_name_level(self, c: str):
        if c == '>' and self.name_level < 3:
            self.name_level += 1
        elif c == ')':
            self.names[self.name_level] = self.name_buffer()
            self.name_level = 0
            self.state = BEFORE
            self.depth = 0
        else:
            self.name_buffer.append(c)
            self.state = PEPTIDOFORM_NAME_TEXT

    def handle_name_text(self, c: str):
        if c == ')':
            self.depth -= 1
            if self.depth <= 0:
                name = self.name_buffer()
                self.names[self.name_level] = name
                self.state = BEFORE
            else:
                self.name_buffer.append(c)
        elif c == '(':
            self.depth += 1
            self.name_buffer.append(c)
        else:
            self.name_buffer.append(c)

    def step(self) -> bool:
        if self.index < self.length:
            c = self.sequence[self.index]

            # Initial state prior to sequence content
            if self.state == BEFORE:
                self.handle_before(c)
            # The body of the amino acid sequence.
            elif self.state == SEQ or self.state == INTERVAL_INIT:
                self.handle_seq(c)

            # Tag parsing which rely on `current_tag` to buffer tokens.
            elif (
                self.state == TAG
                or self.state == TAG_BEFORE
                or self.state == TAG_AFTER
                or self.state == GLOBAL
                or self.state == INTERVAL_TAG
            ):
                self.handle_tag(c)

            # Handle transition to fixed modifications or isotope labeling from opening signal.
            elif self.state == FIXED:
                self.handle_fixed(c)
            # Handle fixed isotope rules, which rely on `current_tag` to buffer tokens
            elif self.state == ISOTOPE:
                self.handle_isotope(c)
            # Handle labile modifications, which rely on `current_tag` to buffer tokens
            elif self.state == LABILE:
                self.handle_labile(c)
            # The intermediate state between an interval tag and returning to sequence parsing.
            # A new tag may start immediately, leading to it being appended to the interval instead
            # instead of returning to the primary sequence. Because this state may also occur at the
            # end of a sequence, it must also handle sequence-terminal transitions like C-terminal tags,
            # charge states, and the like.
            elif self.state == POST_INTERVAL_TAG:
                self.handle_post_interval_tag(c)
            # An intermediate state for discriminating which type of tag-before-sequence type
            # we just finished parsing.
            elif self.state == POST_TAG_BEFORE:
                self.handle_post_tag_before(c)
            elif self.state == UNLOCALIZED_COUNT:
                self.handle_unlocalized_count(c)
            elif self.state == POST_GLOBAL:
                self.handle_post_global(c)
            elif self.state == POST_GLOBAL_AA:
                self.handle_post_global_aa(c)
            elif self.state == POST_TAG_AFTER:
                self.handle_post_tag_after(c)
            elif self.state == CHARGE_START:
                self.handle_charge_start(c)
            elif self.state == CHARGE_NUMBER:
                self.handle_charge_number(c)
            elif self.state == ADDUCT_START:
                self.handle_adduct_start(c)
            elif self.state == ADDUCT_END:
                self.handle_adduct_end(c)
            elif self.state == PEPTIDOFORM_NAME_LEVEL:
                self.handle_name_level(c)
            elif self.state == PEPTIDOFORM_NAME_TEXT:
                self.handle_name_text(c)
            else:
                raise ProFormaError(
                    f"Error In State {self.state}, unexpected {c} found at index {self.i}",
                    self.i,
                    self.state,
                )
            self.index += 1
        return self.index < self.length

    def finish(
        self,
    ) -> Tuple[List[Tuple[str, Optional[List[TagBase]]]], Dict[str, Any]]:
        """
        Post-process the parser's accumulated parsed token data and return the parsed
        sequence and metadata.

        Returns
        -------
        sequence : List[Tuple[str, Optional[List[TagBase]]]]
            The primary amino acid sequence of the ProForma string
        metadata : Dict[str, Any]
            All other information outside the main sequence, including unlocalized, labile, or global modifications,
            names, charge states, and more.
        """
        if self.charge_buffer:
            charge_number = self.charge_buffer()
            if self.adduct_buffer:
                adducts = self.adduct_buffer()
            else:
                adducts = None
            charge_state = ChargeState(charge_number, adducts)
        elif self.adduct_buffer:
            adducts = self.adduct_buffer()
            charge_state = ChargeState.from_adducts(adducts)
        else:
            charge_state = None
        if self.current_aa:
            self.pack_sequence_position()
        if self.state in (
            ISOTOPE,
            TAG,
            TAG_AFTER,
            TAG_BEFORE,
            LABILE,
        ):
            raise ProFormaError(
                f"Error In State {self.state}, unclosed group reached end of string!",
                self.i,
                self.state,
            )
        return self.positions, {
            "n_term": self.n_term,
            "c_term": self.c_term,
            "unlocalized_modifications": self.unlocalized_modifications,
            "labile_modifications": self.labile_modifications,
            "fixed_modifications": self.fixed_modifications,
            "intervals": self.intervals,
            "isotopes": self.isotopes,
            "group_ids": sorted(self.current_tag.group_ids),
            "charge_state": charge_state,
            "names": self.names
        }

    def parse(self):
        while self.step():
            pass
        return self.finish()

    def __call__(self, *args, **kwds):
        return self.parse()


def parse(sequence: str, **kwargs) -> Tuple[List[Tuple[str, Optional[List[TagBase]]]], Dict[str, Any]]:
    """
    Tokenize a ProForma sequence into a sequence of amino acid+tag positions, and a
    mapping of sequence-spanning modifiers.

    .. note::
        This is a state machine parser, but with certain sub-state paths
        unrolled to avoid an explosion of formal intermediary states.

    Parameters
    ----------
    sequence: str
        The sequence to parse
    **kwargs :
        Forwarded to :class:`Parser`

    Returns
    -------
    parsed_sequence: list[tuple[str, list[TagBase]]]
        The (amino acid: str, TagBase or None) pairs denoting the positions along the primary sequence
    modifiers: dict
        A mapping listing the labile modifications, fixed modifications, stable isotopes, unlocalized
        modifications, tagged intervals, and group IDs
    """
    parser = Parser(sequence, **kwargs)
    return parser.parse()


def _parse(sequence):
    '''Tokenize a ProForma sequence into a sequence of amino acid+tag positions, and a
    mapping of sequence-spanning modifiers.

    .. warning::
        This is the older parser which was designed for ProForma v2.0. There are some syntactic
        constructs that were introduced in v2.1 that are not compatible. This function is retained
        for handling sequences in the older format only.

    .. note::
        This is a state machine parser, but with certain sub-state paths
        unrolled to avoid an explosion of formal intermediary states.

    Parameters
    ----------
    sequence: str
        The sequence to parse

    Returns
    -------
    parsed_sequence: list[tuple[str, list[TagBase]]]
        The (amino acid: str, TagBase or None) pairs denoting the positions along the primary sequence
    modifiers: dict
        A mapping listing the labile modifications, fixed modifications, stable isotopes, unlocalized
        modifications, tagged intervals, and group IDs
    '''
    labile_modifications = []
    fixed_modifications = []
    unlocalized_modifications = []
    intervals = []
    isotopes = []

    n_term = None
    c_term = None

    i = 0
    n = len(sequence)

    positions = []
    state = BEFORE
    depth = 0

    current_aa = None
    current_tag = TagParser()
    current_interval = None
    current_unlocalized_count = NumberParser()
    current_aa_targets = StringParser()

    charge_buffer = None
    adduct_buffer = None

    # A mostly context free finite state machine unrolled
    # by hand.
    while i < n:
        c = sequence[i]
        i += 1
        # Initial state prior to sequence content
        if state == BEFORE:
            if c == '[':
                state = TAG_BEFORE
                depth = 1
            elif c == '{':
                state = LABILE
                depth = 1
            elif c == '<':
                state = FIXED
            elif c in VALID_AA:
                current_aa = c
                state = SEQ
            else:
                raise ProFormaError(
                    f"Error In State {state}, unexpected {c} found at index {i}", i, state)
        # The body of the amino acid sequence.
        elif state == SEQ or state == INTERVAL_INIT:
            if state == INTERVAL_INIT:
                state = SEQ
                if c == '?':
                    if current_interval is not None:
                        current_interval.ambiguous = True
                    continue
            if c in VALID_AA:
                if current_aa is not None:
                    positions.append((current_aa, current_tag() if current_tag else None))
                current_aa = c
            elif c == '[':
                state = TAG
                if current_tag:
                    current_tag.bound()
                depth = 1
            elif c == '(':
                if current_interval is not None:
                    raise ProFormaError(
                        ("Error In State {state}, nested range found at index {i}. "
                         "Nested ranges are not yet supported by ProForma.").format(
                            **locals()), i, state)
                current_interval = TaggedInterval(len(positions) + 1)
                state = INTERVAL_INIT
            elif c == ')':
                positions.append(
                    (current_aa, current_tag() if current_tag else None))
                current_aa = None
                if current_interval is None:
                    raise ProFormaError("Error In State {state}, unexpected {c} found at index {i}".format(**locals()), i, state)
                else:
                    current_interval.end = len(positions)
                    if i < n and sequence[i] == '[':
                        i += 1
                        depth = 1
                        state = INTERVAL_TAG
                    else:
                        intervals.append(current_interval)
                        current_interval = None
            elif c == '-':
                if current_aa:
                    positions.append((current_aa, current_tag() if current_tag else None))
                    current_aa = None
                state = TAG_AFTER
                if i >= n or sequence[i] != '[':
                    raise ProFormaError("Missing Closing Tag", i, state)
                i += 1
                depth = 1
            elif c == '/':
                state = CHARGE_START
                charge_buffer = NumberParser()
            elif c == '+':
                raise ProFormaError(
                    "Error In State {state}, {c} found at index {i}. Chimeric representation not supported".format(**locals()), i, state)
            else:
                raise ProFormaError("Error In State {state}, unexpected {c} found at index {i}".format(**locals()), i, state)
        # Tag parsing which rely on `current_tag` to buffer tokens.
        elif state == TAG or state == TAG_BEFORE or state == TAG_AFTER or state == GLOBAL or state == INTERVAL_TAG:
            if c == '[':
                depth += 1
                current_tag.append(c)
            elif c == ']':
                depth -= 1
                if depth <= 0:
                    depth = 0
                    if state == TAG:
                        state = SEQ
                    elif state == TAG_BEFORE:
                        state = POST_TAG_BEFORE
                    elif state == TAG_AFTER:
                        c_term = current_tag()
                        state = POST_TAG_AFTER
                    elif state == GLOBAL:
                        state = POST_GLOBAL
                    elif state == INTERVAL_TAG:
                        state = POST_INTERVAL_TAG
                        depth = 0
                else:
                    current_tag.append(c)
            else:
                current_tag.append(c)
        # Handle transition to fixed modifications or isotope labeling from opening signal.
        elif state == FIXED:
            if c == '[':
                state = GLOBAL
            else:
                # Do validation here
                state = ISOTOPE
                current_tag.reset()
                current_tag.append(c)
        # Handle fixed isotope rules, which rely on `current_tag` to buffer tokens
        elif state == ISOTOPE:
            if c != '>':
                current_tag.append(c)
            else:
                # Not technically a tag, but exploits the current buffer
                isotopes.append(StableIsotope(''.join(current_tag)))
                current_tag.reset()
                state = BEFORE
        # Handle labile modifications, which rely on `current_tag` to buffer tokens
        elif state == LABILE:
            if c == '{':
                depth += 1
            elif c == '}':
                depth -= 1
                if depth <= 0:
                    depth = 0
                    labile_modifications.append(current_tag()[0])
                    state = BEFORE
            else:
                current_tag.append(c)
        # The intermediate state between an interval tag and returning to sequence parsing.
        # A new tag may start immediately, leading to it being appended to the interval instead
        # instead of returning to the primary sequence. Because this state may also occur at the
        # end of a sequence, it must also handle sequence-terminal transitions like C-terminal tags,
        # charge states, and the like.
        elif state == POST_INTERVAL_TAG:
            if c == '[':
                current_tag.bound()
                state = INTERVAL_TAG
            elif c in VALID_AA:
                current_aa = c
                current_interval.tags = current_tag()
                intervals.append(current_interval)
                current_interval = None
                state = SEQ
            elif c == '-':
                state = TAG_AFTER
                if i >= n or sequence[i] != '[':
                    raise ProFormaError("Missing Closing Tag", i, state)
                i += 1
                depth = 1
            elif c == '/':
                state = CHARGE_START
                charge_buffer = NumberParser()
            elif c == '+':
                raise ProFormaError(
                    "Error In State {state}, {c} found at index {i}. Chimeric representation not supported".format(**locals()), i, state)
            else:
                raise ProFormaError(
                    "Error In State {state}, unexpected {c} found at index {i}".format(**locals()), i, state)
        # An intermediate state for discriminating which type of tag-before-sequence type
        # we just finished parsing.
        elif state == POST_TAG_BEFORE:
            if c == '?':
                unlocalized_modifications.append(current_tag()[0])
                state = BEFORE
            elif c == '-':
                n_term = current_tag()
                state = BEFORE
            elif c == '^':
                state = UNLOCALIZED_COUNT
            else:
                raise ProFormaError(
                    "Error In State {state}, unexpected {c} found at index {i}".format(**locals()), i, state)
        elif state == UNLOCALIZED_COUNT:
            if c.isdigit():
                current_unlocalized_count.append(c)
            elif c == '[':
                state = TAG_BEFORE
                depth = 1
                tag = current_tag()[0]
                multiplicity = current_unlocalized_count()
                for _ in range(multiplicity):
                    unlocalized_modifications.append(tag)
            elif c == '?':
                state = BEFORE
                tag = current_tag()[0]
                multiplicity = current_unlocalized_count()
                for _ in range(multiplicity):
                    unlocalized_modifications.append(tag)
            else:
                raise ProFormaError(
                    "Error In State {state}, unexpected {c} found at index {i}".format(**locals()), i, state)
        elif state == POST_GLOBAL:
            if c == '@':
                state = POST_GLOBAL_AA
            else:
                raise ProFormaError(
                    ("Error In State {state}, fixed modification detected without "
                     "target amino acids found at index {i}").format(**locals()), i, state)
        elif state == POST_GLOBAL_AA:
            if c in VALID_AA or c in TERMINAL_SPEC_CHARS:
                current_aa_targets.append(c)
            elif c == ',':
                # the next character should be another amino acid
                current_aa_targets.bound()
            elif c == '>':
                try:
                    v = current_aa_targets()
                    fixed_modifications.append(
                        ModificationRule(current_tag()[0], v))
                except PyteomicsError as err:
                    raise ProFormaError(
                        (
                            "Error In State {state}, fixed modification detected invalid "
                            "target found at index {i}: {err}"
                        ).format(state=state, i=i, err=err),
                        i,
                        state,
                    )
                state = BEFORE
            else:
                raise ProFormaError(
                    ("Error In State {state}, unclosed fixed modification rule").format(**locals()), i, state)
        elif state == POST_TAG_AFTER:
            if c == '/':
                state = CHARGE_START
                charge_buffer = NumberParser()
            elif c == '+':
                raise ProFormaError(
                    "Error In State {state}, {c} found at index {i}. Chimeric representation not supported".format(**locals()), i, state)
        elif state == CHARGE_START:
            if c in '+-':
                charge_buffer.append(c)
                state = CHARGE_NUMBER
            elif c.isdigit():
                charge_buffer.append(c)
                state = CHARGE_NUMBER
            elif c == '/':
                state = ParserStateEnum.inter_chain_cross_link_start
                raise ProFormaError("Inter-chain cross-linked peptides are not yet supported", i, state)
            else:
                raise ProFormaError(
                    "Error In State {state}, unexpected {c} found at index {i}".format(**locals()), i, state)
        elif state == CHARGE_NUMBER:
            if c.isdigit():
                charge_buffer.append(c)
            elif c == "[":
                state = ADDUCT_START
                adduct_buffer = AdductParser()
            else:
                raise ProFormaError(
                    "Error In State {state}, unexpected {c} found at index {i}".format(**locals()), i, state)
        elif state == ADDUCT_START:
            if c.isdigit() or c in "+-" or c.isalpha():
                adduct_buffer.append(c)
            elif c == ',':
                adduct_buffer.bound()
            elif c == ']':
                state = ADDUCT_END
        elif state == ADDUCT_END:
            if c == '+':
                raise ProFormaError(
                    "Error In State {state}, {c} found at index {i}. Chimeric representation not supported".format(**locals()), i, state)
        else:
            raise ProFormaError("Error In State {state}, unexpected {c} found at index {i}".format(**locals()), i, state)
    if charge_buffer:
        charge_number = charge_buffer()
        if adduct_buffer:
            adducts = adduct_buffer()
        else:
            adducts = None
        charge_state = ChargeState(charge_number, adducts)
    else:
        charge_state = None
    if current_aa:
        positions.append((current_aa, current_tag() if current_tag else None))
    if state in (ISOTOPE, TAG, TAG_AFTER, TAG_BEFORE, LABILE, ):
        raise ProFormaError("Error In State {state}, unclosed group reached end of string!".format(**locals()), i, state)
    return positions, {
        'n_term': n_term,
        'c_term': c_term,
        'unlocalized_modifications': unlocalized_modifications,
        'labile_modifications': labile_modifications,
        'fixed_modifications': fixed_modifications,
        'intervals': intervals,
        'isotopes': isotopes,
        'group_ids': sorted(current_tag.group_ids),
        'charge_state': charge_state,
    }


def to_proforma(
    sequence,
    n_term: Optional[List[TagBase]] = None,
    c_term: Optional[List[TagBase]] = None,
    unlocalized_modifications: Optional[List[TagBase]] = None,
    labile_modifications: Optional[List[TagBase]] = None,
    fixed_modifications: Optional[List[TagBase]] = None,
    intervals: Optional[List[TaggedInterval]]=None,
    isotopes: Optional[List[StableIsotope]] = None,
    charge_state: Optional[ChargeState]=None,
    group_ids: Iterable[str]=None,
    names: Optional[Dict[int, str]] = None,
):
    '''Convert a sequence plus modifiers into formatted text following the
    ProForma specification.

    Parameters
    ----------
    sequence : list[tuple[str, TagBase]]
        The primary sequence of the peptidoform/proteoform to render
    n_term : Optional[TagBase]
        The N-terminal modification, if any.
    c_term : Optional[TagBase]
        The C-terminal modification, if any.
    unlocalized_modifications : Optional[list[TagBase]]
        Any modifications which aren't assigned to a specific location.
    labile_modifications : Optional[list[TagBase]]
        Any labile modifications
    fixed_modifications : Optional[list[ModificationRule]]
        Any fixed modifications
    intervals : Optional[list[TaggedInterval]]
        A list of modified intervals, if any
    isotopes : Optional[list[StableIsotope]]
        Any global stable isotope labels applied
    charge_state : Optional[ChargeState]
        An optional charge state value
    group_ids : Optional[list[str]]
        Any group identifiers. This parameter is currently not used.

    Returns
    -------
    str
    '''
    if names is None:
        names = {}
    buffer = []
    if 3 in names:
        buffer.append("(>>>")
        buffer.append(names[3])
        buffer.append(")")
    if isotopes:
        for iso in isotopes:
            buffer.append(str(iso))
    if fixed_modifications:
        for rule in fixed_modifications:
            buffer.append(str(rule))
    if 2 in names:
        buffer.append("(>>")
        buffer.append(names[2])
        buffer.append(")")
    if 1 in names:
        buffer.append("(>")
        buffer.append(names[1])
        buffer.append(")")
    primary = deque()
    for aa, tags in sequence:
        if not tags:
            primary.append(str(aa))
        else:
            primary.append(str(aa) + "".join(["[{0!s}]".format(t) for t in tags]))
    if intervals:
        for iv in sorted(intervals, key=lambda x: x.start):
            if iv.ambiguous:
                primary[iv.start] = "(?" + primary[iv.start]
            else:
                primary[iv.start] = "(" + primary[iv.start]

            terminator = "{0!s})".format(primary[iv.end - 1])
            if iv.tags:
                terminator += "".join("[{!s}]".format(t) for t in iv.tags)
            primary[iv.end - 1] = terminator
    if n_term:
        primary.appendleft("".join("[{!s}]".format(t) for t in n_term) + "-")
    if c_term:
        primary.append("-" + "".join("[{!s}]".format(t) for t in c_term))
    if charge_state:
        primary.append("/{!s}".format(charge_state))
    if labile_modifications:
        primary.extendleft(["{{{!s}}}".format(m) for m in labile_modifications])
    if unlocalized_modifications:
        primary.appendleft("?")
        primary.extendleft(["[{!s}]".format(m) for m in unlocalized_modifications])
    primary.appendleft("".join(buffer))
    return "".join(primary)


class _ProFormaProperty(Generic[T]):
    def __init__(self, name):
        self.name = name

    def __get__(self, obj, cls) -> T:
        return obj.properties[self.name]

    def __set__(self, obj, value: T):
        obj.properties[self.name] = value

    def __repr__(self):
        template = "{self.__class__.__name__}({self.name!r})"
        return template.format(self=self)


class ProForma(object):
    '''Represent a parsed ProForma sequence.

    The preferred way to instantiate this class is via the :meth:`parse`
    method.

    Attributes
    ----------
    sequence : list[tuple[str, List[TagBase]]]
        The list of (amino acid, tag collection) pairs making up the primary sequence of the
        peptide.
    isotopes : list[StableIsotope]
        A list of any stable isotope rules that apply to this peptide
    charge_state : int, optional
        An optional charge state that may have been provided
    intervals : list[Interval]
        Any annotated intervals that contain either sequence ambiguity or a
        tag over that interval.
    labile_modifications : list[ModificationBase]
        Any modifications that were parsed as labile, and may not appear at
        any location on the peptide primary sequence.
    unlocalized_modifications : list[ModificationBase]
        Any modifications that were not localized but may be attached to peptide
        sequence evidence.
    n_term : list[ModificationBase]
        Any modifications on the N-terminus of the peptide
    c_term : list[ModificationBase]
        Any modifications on the C-terminus of the peptide
    group_ids : set
        The collection of all groupd identifiers on this sequence.
    mass : float
        The computed mass for the fully modified peptide, including labile
        and unlocalized modifications. **Does not include stable isotopes at this time**
    '''

    sequence: List[Tuple[str, Optional[List[TagBase]]]]
    properties: Dict[str, Any]

    def __init__(self, sequence, properties):
        """
        Initialize a :class:`ProForma` instance from a parse tree.

        To construct an instance from a string directly, see :meth:`ProForma.parse`.

        See Also
        --------
        :meth:`ProForma.parse`
        """
        self.sequence = sequence
        self.properties = properties

    isotopes = _ProFormaProperty[List[StableIsotope]]('isotopes')
    _charge_state = _ProFormaProperty('charge_state')

    intervals = _ProFormaProperty[List[TaggedInterval]]('intervals')
    fixed_modifications = _ProFormaProperty[List[ModificationRule]]("fixed_modifications")
    labile_modifications = _ProFormaProperty[List[TagBase]]('labile_modifications')
    unlocalized_modifications = _ProFormaProperty[List[TagBase]]("unlocalized_modifications")

    n_term = _ProFormaProperty[List[TagBase]]("n_term")
    c_term = _ProFormaProperty[List[TagBase]]("c_term")

    group_ids = _ProFormaProperty('group_ids')

    def __str__(self):
        return to_proforma(self.sequence, **self.properties)

    def __repr__(self):
        return "{self.__class__.__name__}({self.sequence}, {self.properties})".format(self=self)

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, i):
        if isinstance(i, slice):
            props = self.properties.copy()
            ivs = []
            for iv in props['intervals']:
                iv = iv._update_coordinates_sliced(
                    i.start, i.stop)
                if iv is None:
                    continue
                ivs.append(iv)
            props['intervals'] = ivs

            if not (i.start is None or i.start == 0):
                props['n_term'] = None
            n = len(self)
            if not (i.stop is None or i.stop >= n):
                props['c_term'] = None

            subseq = self.__class__(self.sequence[i], props)
            if subseq.group_ids:
                kept_group_ids = []
                for group_id in subseq.group_ids:
                    tag_hits = subseq.find_tags_by_id(group_id, include_position=True)
                    if not tag_hits:
                        continue
                    kept_group_ids.append(group_id)
                    # We sliced the sequence, but only the localization markers were captured,
                    # not the actual modification definition. Update the first occurrence of the
                    # localization marker with a group id marked modification tag.
                    if all(not isinstance(v, LocalizationMarker) for _, v in tag_hits):
                        i = tag_hits[0]
                        val: TagBase
                        for val in self.find_tags_by_id(group_id, include_position=False):
                            if not isinstance(val, LocalizationMarker):
                                val = val.copy()
                                for j, tag in enumerate(subseq[i][1]):
                                    if tag.group_id == group_id:
                                        subseq[i][1][j] = val
            return subseq
        else:
            return self.sequence[i]

    def __setitem__(self, i, val: Tuple[str, Optional[List[TagBase]]]):
        self.sequence[i] = val

    def __eq__(self, other):
        if isinstance(other, str):
            return str(self) == other
        elif other is None:
            return False
        else:
            return self.sequence == other.sequence and self.properties == other.properties

    def __ne__(self, other):
        return not self == other

    def __iter__(self):
        return iter(self.sequence)

    @property
    def charge_state(self) -> Optional[ChargeState]:
        z = self._charge_state
        return z

    @classmethod
    def parse(cls, string, **kwargs):
        '''Parse a ProForma string.

        Parameters
        ----------
        string : str
            The string to parse
        **kwargs :
            Forwarded to :class:`Parser`
        Returns
        -------
        ProForma
        '''
        return cls(*parse(string, **kwargs))

    @property
    def mass(self):
        '''
        Compute the *total* monoisotopic neutral mass of the peptidoform.

        This does not include the adduct.
        '''
        mass = 0.0

        fixed_modifications = self.properties['fixed_modifications']

        n_term_v = 0
        c_term_v = len(self) - 1
        for i, position in enumerate(self.sequence):
            aa = position[0].upper()
            try:
                mass += std_aa_mass[aa]
            except KeyError:
                warnings.warn("%r does not have an exact mass" % (aa, ))
            n_term = i == n_term_v
            c_term = i == c_term_v
            for rule in fixed_modifications:
                if rule.is_valid(aa, n_term, c_term):
                    if rule.modification_tag.has_mass():
                        mass += rule.modification_tag.mass
            tags = position[1]
            if tags:
                for tag in tags:
                    if tag.has_mass():
                        mass += tag.mass
        for mod in self.properties['labile_modifications']:
            mass += mod.mass
        for mod in self.properties['unlocalized_modifications']:
            mass += mod.mass
        if self.properties.get('n_term'):
            for mod in self.properties['n_term']:
                if mod.has_mass():
                    mass += mod.mass
        mass += calculate_mass(formula="H")
        if self.properties.get('c_term'):
            for mod in self.properties['c_term']:
                if mod.has_mass():
                    mass += mod.mass

        mass += calculate_mass(formula="OH")
        for iv in self.properties['intervals']:
            for tag in iv.tags:
                if tag.has_mass():
                    mass += tag.mass
        return mass

    def mz(self, charge_state: Union[int, ChargeState, None] = None, **kwargs) -> float:
        """
        Compute the *total* m/z of the peptidoform in the specified charge state, or fall back
        to the peptidoform ion's defined charge state and adduction.

        This method first tries to get the composition of the peptidoform ion with :meth:`composition`
        and then forwards ``kwargs`` to :meth:`Composition.mass` to compute m/z with full flexibility,
        but if that fails due to missing modification compositions, this method falls back to directly
        computing monoisotopic mass and uses the charge state to get the m/z.

        .. warning::
            If no charge state of any kind is available, this will raise a :py:class:`MissingChargeStateError`.

        Parameters
        ----------
        charge_state : int or :class:`ChargeState`, optional
            The charge state either as in integer number of protons gained/lost,
            or a :class:`ChargeState` instance. If not provided, :attr:`charge_state`
            will be used.
        **kwargs :
            Forwarded to :meth:`Composition.mass`

        Returns
        -------
        float
        """
        if charge_state is None:
            charge_state = self.charge_state
        elif isinstance(charge_state, Integral):
            charge_state = ChargeState(int(charge_state))
        elif not isinstance(charge_state, ChargeState):
            raise TypeError(
                f"Expected a charge state-like type, got {type(charge_state)}"
            )
        if charge_state is None:
            raise MissingChargeStateError(
                f"Requested an m/z value without providing a charge state and the peptidoform {self!r} does "
                "not have a charge state itself."
            )
        try:
            composition = self.composition(include_charge=charge_state, ignore_missing=False)
            return composition.mass(**kwargs)
        except ProFormaError:
            charge_carrier_mass, charge = charge_state.for_mz_calculation()
            return (self.mass + charge_carrier_mass) / abs(charge)

    def fragments(self, ion_shift, charge=1, reverse=None, include_labile=True, include_unlocalized=True):
        """
        The function generates all possible fragments of the requested
        series type.

        Parameters
        ----------
        ion_shift : float or str
            The mass shift of the ion series, or the name of the ion series
        charge : int
            The charge state of the theoretical fragment masses to generate.
            Defaults to 1+. If 0 is passed, neutral masses will be returned.
        reverse : bool, optional
            Whether to fragment from the N-terminus (``False``) or C-terminus (``True``).
            If ``ion_shift`` is a :class:`str`, the terminal will be inferred from
            the series name. Otherwise, defaults to ``False``.
        include_labile : bool, optional
            Whether or not to include dissociated modification masses.
            Defaults to ``True``
        include_unlocalized : bool, optional
            Whether or not to include unlocalized modification masses.
            Defaults to ``True``

        Returns
        -------
        np.ndarray

        Examples
        --------

        >>> p = proforma.ProForma.parse("PEPTIDE")
        >>> p.fragments('b', charge=1)
        array([ 98.06004032, 227.1026334 , 324.15539725, 425.20307572,
                538.2871397 , 653.31408272])
        >>> p.fragments('y', charge=1)
        array([148.06043424, 263.08737726, 376.17144124, 477.21911971,
               574.27188356, 703.31447664])

        """
        if isinstance(ion_shift, str):
            if ion_shift[0] in 'xyz':
                reverse = True
            ion_shift = std_ion_comp[ion_shift].mass(absolute=False)

        n = len(self.sequence)
        masses = _array('d')

        mass = 0
        mass += ion_shift

        fixed_modifications = self.properties['fixed_modifications']

        intervals = self.intervals
        if intervals:
            intervals = sorted(intervals, key=lambda x: x.start, reverse=reverse)
        intervals = deque(intervals)

        if not include_labile:
            for mod in self.properties['labile_modifications']:
                mass += mod.mass

        if not reverse:
            if self.properties.get('n_term'):
                for mod in self.properties['n_term']:
                    if mod.has_mass():
                        mass += mod.mass
        else:
            if self.properties.get('c_term'):
                for mod in self.properties['c_term']:
                    if mod.has_mass():
                        mass += mod.mass

        if include_unlocalized:
            for mod in self.properties['unlocalized_modifications']:
                if mod.has_mass():
                    mass += mod.mass

        mass += _WATER_MASS

        if not reverse:
            iterator = (iter(range(0, n - 1)))
            n_term_v = 0
            c_term_v = n - 1
        else:
            iterator = (reversed(range(1, n)))
            n_term_v = n - 1
            c_term_v = 0

        for i in iterator:
            position = self.sequence[i]

            aa = position[0].upper()
            if aa != 'X':
                try:
                    mass += std_aa_mass[aa]
                except KeyError:
                    warnings.warn("%r does not have an exact mass" % (aa, ))

            n_term = i == n_term_v
            c_term = i == c_term_v
            for rule in fixed_modifications:
                if rule.is_valid(aa, n_term, c_term):
                    if rule.modification_tag.has_mass():
                        mass += rule.modification_tag.mass

            tags = position[1]
            if tags:
                for tag in tags:
                    if tag.has_mass():
                        mass += tag.mass

            while intervals and intervals[0].contains(i):
                iv = intervals.popleft()
                for tag in iv.tags:
                    if tag.has_mass():
                        mass += tag.mass

            masses.append(mass)

        if np is not None:
            masses = np.asarray(masses)
            if charge != 0:
                return mass_charge_ratio(masses, charge)
            return masses
        if charge != 0:
            for i, mass in enumerate(masses):
                masses[i] = mass_charge_ratio(mass, charge)
        return masses

    def find_tags_by_id(self, tag_id, include_position=True):
        '''Find all occurrences of a particular tag ID

        Parameters
        ----------
        tag_id : str
            The tag ID to search for
        include_position : bool
            Whether or not to return the locations for matched
            tag positions

        Returns
        -------
        list[tuple[Any, TagBase]] or list[TagBase]
        '''
        if not tag_id.startswith("#"):
            tag_id = "#" + tag_id
        matches = []
        for i, (_token, tags) in enumerate(self.sequence):
            if tags:
                for tag in tags:
                    if tag.group_id == tag_id:
                        if include_position:
                            matches.append((i, tag))
                        else:
                            matches.append(tag)
        for iv in self.properties['intervals']:
            if iv.tag.group_id == tag_id:
                matches.append((iv, iv.tag) if include_position else iv.tag)
        for ulmod in self.properties['unlocalized_modifications']:
            if ulmod.group_id == tag_id:
                matches.append(('unlocalized_modifications', ulmod)
                               if include_position else ulmod)
        for lamod in self.properties['labile_modifications']:
            if lamod.group_id == tag_id:
                matches.append(('labile_modifications', lamod)
                               if include_position else lamod)
        return matches

    @property
    def tags(self):
        return [tag for tags_at in [pos[1] for pos in self if pos[1]] for tag in tags_at]

    def generate_proteoforms(self):
        return iter(ProteoformCombinator(self))

    def copy(self):
        sequence = []
        for (aa, tags) in self:
            if tags:
                tags = [t.copy() for t in tags]
            sequence.append((aa, tags))
        properties = self.properties.copy()
        for k in [
            "labile_modifications",
            "n_term",
            "c_term",
            "unlocalized_modifications",
            "fixed_modifications",
            "isotopes",
        ]:
            properties[k] = [v.copy() for v in properties[k]]
        properties['names'] = properties['names'].copy()
        return self.__class__(sequence, properties)

    def composition(self, include_charge: Union[bool, ChargeState]=False, aa_comp=None, ignore_missing=False) -> Composition:
        '''
        Calculate the elemental composition of the ProForma sequence.

        Parameters
        ----------
        include_charge : bool or :class:`ChargeState`, optional
            If True, then :attr:`charge_state` will be included in the composition.
            If a :class:`ChargeState` instance is passed, this charge and adduction
            will be included instead. Otherwise, composition of the neutral molecule
            will be returned. Defaults to False.
        aa_comp : dict, optional
            A dictionary mapping amino acid symbols to their respective
            compositions. If not provided, the standard amino acid composition
            will be used. ``X`` *always* has a mass of 0.0, regardless of this
            argument.
        ignore_missing : bool, optional
            If True, tags with missing composition will be silently ignored. If False (default),
            a :py:class:`CompositionNotFoundError` will be raised.

            .. note::
                Amino acids not found in `aa_mass` will result in errors even with `ignore_missing=True`.

        Returns
        -------
        Composition
            :py:class:`Composition` object representing the composition of the ProForma sequence.
        '''
        if ignore_missing:
            def get_comp(tag):
                try:
                    return tag.composition or Composition({})
                except AttributeError:
                    return Composition({})
        else:
            def get_comp(tag):
                try:
                    comp = tag.composition
                except AttributeError as e:
                    raise CompositionNotFoundError(f'No composition found for tag {tag}') from e
                if comp is None:
                    raise CompositionNotFoundError(f'No composition found for tag {tag}')
                return comp

        comp = Composition()
        if aa_comp is None:
            aa_comp = std_aa_comp
        try:
            for i, (aa, tags) in enumerate(self.sequence):
                if aa != 'X':
                    comp += aa_comp[aa]
                for tag in tags or []:
                    if not tag.has_composition():
                        if isinstance(tag, MassModification) and not ignore_missing:
                            raise CompositionNotFoundError(f"No composition found for tag {tag}")
                        continue
                    comp += get_comp(tag)
                for rule in self.fixed_modifications:
                    if rule.is_valid(aa, i == 0, i == len(self.sequence) - 1):
                        comp += get_comp(rule.modification_tag)
            for tag in chain(self.labile_modifications, self.unlocalized_modifications,
                             chain.from_iterable(interval.tags for interval in self.intervals)):
                comp += get_comp(tag)
        except KeyError as e:
            raise CompositionNotFoundError(f'No composition found for amino acid {aa}') from e
        if self.n_term:
            comp += get_comp(self.n_term)
        else:
            comp['H'] += 1  # Add hydrogen for N-terminus
        if self.c_term:
            comp += get_comp(self.c_term)
        else:
            comp += Composition({'O': 1, 'H': 1})  # Add -OH for C-terminus
        if include_charge and self.charge_state:
            comp += self.charge_state.composition()
        return comp


class GeneratorModificationRuleDirective:
    """
    A helper for :class:`ProteoformCombinator` that maps modification rules to sequence locations.

    This type probably shouldn't be created directly.

    TODO: limit directives, position group_id tags
    """
    rule: ModificationRule
    region: Optional[TaggedInterval] = None
    colocal_known: bool = False
    colocal_unknown: bool = False
    limit: int = 1

    def __init__(self, rule, region=None, colocal_known: bool = False, colocal_unknown: bool = False, limit: int = 1):
        self.rule = rule
        self.region = region
        self.colocal_known = colocal_known
        self.colocal_unknown = colocal_unknown
        self.limit = limit

    def create(self) -> TagBase:
        return self.rule.modification_tag.copy()

    def __repr__(self):
        return f"{self.__class__.__name__}({self.rule}, {self.region}, {self.colocal_known}, {self.colocal_unknown})"

    def find_positions(self, sequence: ProForma) -> List[int]:
        n = len(sequence) - 1
        positions = []
        group_id = self.rule.modification_tag.group_id
        for i, (aa, tags) in enumerate(sequence):
            if self.region and not self.region.contains(i):
                continue
            elif group_id is not None:
                if not tags:
                    continue
                for tag in tags:
                    # TODO: Implement combinatoric limits here
                    if tag.group_id == group_id:
                        positions.append(i)
            elif self.rule.is_valid(aa, i == 0, i == n):
                if not tags:
                    positions.append(i)
                else:
                    known = []
                    unknown = []
                    for tag in tags:
                        if isinstance(tag, (ModificationBase, MassModification)):
                            if tag._generated in (ModificationSourceType.Explicit, ModificationSourceType.Constant):
                                known.append(tag)
                            elif tag._generated == ModificationSourceType.Generated:
                                unknown.append(tag)
                    total_at = len(known) + len(unknown)
                    can_known = (known and self.colocal_known) or not known
                    can_unknown = (unknown and self.colocal_unknown) or not unknown
                    if ((can_known and can_unknown) or (not can_known and not can_unknown)) and total_at < self.limit:
                        positions.append(i)
        return positions

    @classmethod
    def from_unlocalized_rule(cls, tag: TagBase) -> "GeneratorModificationRuleDirective":
        mod = tag.find_modification()
        if not mod:
            return
        position_constraints = tag.find_tag_type(TagTypeEnum.position_modifier)
        targets = [ModificationTarget(v.value) for v in position_constraints]
        colocal_known = bool(tag.find_tag_type(TagTypeEnum.comkp))
        colocal_unknown = bool(tag.find_tag_type(TagTypeEnum.comup))
        rule = ModificationRule(modification_tag=mod, targets=targets)
        limit = max([t.value for t in tag.find_tag_type(TagTypeEnum.limit)] + [1])
        return cls(rule, None, colocal_known, colocal_unknown, limit)

    @classmethod
    def from_region_rule(cls, region: TaggedInterval) -> List['GeneratorModificationRuleDirective']:
        rules = []
        for tag in (region.tags or []):
            mod = tag.find_modification()
            if not mod:
                continue
            position_constraints = tag.find_tag_type(TagTypeEnum.position_modifier)
            targets = [v.value for v in position_constraints]
            colocal_known = bool(tag.find_tag_type(TagTypeEnum.comkp))
            colocal_unknown = bool(tag.find_tag_type(TagTypeEnum.comup))
            rule = ModificationRule(modification_tag=mod, targets=targets)
            limit = max([t.value for t in tag.find_tag_type(TagTypeEnum.limit)] + [1])
            rules.append(cls(rule, region, colocal_known, colocal_unknown, limit))
        return rules


class ProteoformCombinator:
    """
    Generate combinations of modification (co)localizations for
    modifications that aren't at a fixed position specified in
    the original sequence.

    Attributes
    ----------
    template: :class:`ProForma`
        The template sequence to apply any combination of rules to
    variable_rules: list[:class:`GeneratorModificationRuleDirective`]
        The rules to apply in combinations to the template sequence
    """
    template: ProForma
    variable_rules: List[GeneratorModificationRuleDirective]

    def __init__(self, base_proteoform: ProForma):
        self.template = base_proteoform.copy()
        self.variable_rules = []
        self._extract_rules()
        self._apply_fixed_modifications()
        self._iter = self.generate()

    def _apply_fixed_modifications(self):
        for c in self.template.fixed_modifications:
            rule = GeneratorModificationRuleDirective(c)
            positions = rule.find_positions(self.template)
            for i in positions:
                (aa, tags) = self.template[i]
                if not tags:
                    tags = []
                tag = rule.create()
                if isinstance(tag, (MassModification, ModificationBase)):
                    tag._generated = ModificationSourceType.Constant
                tags.append(tag)
                self.template[i] = (aa, tags)
        self.template.fixed_modifications.clear()

    def _extract_rules(self) -> List[GeneratorModificationRuleDirective]:
        rules = []
        remains = []
        for iv in self.template.intervals:
            block = GeneratorModificationRuleDirective.from_region_rule(iv)
            if block:
                rules.extend(block)
                iv = iv.copy()
                iv.tags = [t for t in iv.tags if not t.is_modification()]
                remains.append(iv)
            else:
                remains.append(iv)
        self.template.intervals = remains

        remains = []
        for rule in self.template.unlocalized_modifications:
            rule_ = GeneratorModificationRuleDirective.from_unlocalized_rule(rule)
            if rule_:
                rules.append(rule_)
            else:
                remains.append(rule)
        self.template.unlocalized_modifications = remains
        self.variable_rules = rules

    def __iter__(self):
        return self

    def __next__(self):
        return next(self._iter)

    def generate(self):
        position_choices = []
        for rule in self.variable_rules:
            positions_for = rule.find_positions(self.template)
            position_choices.append([None] + positions_for)

        for slots in itertools.product(*position_choices):
            template = self.template.copy()
            valid = True
            for rule, idx in zip(self.variable_rules, slots):
                if idx is None:
                    continue
                if idx not in rule.find_positions(template):
                    valid = False
                    break
                (aa, tags) = template[idx]
                if not tags:
                    tags = []
                tag = rule.create()
                tag._generated = ModificationSourceType.Generated
                tags.append(tag)
                template[idx] = (aa, tags)
            if valid:
                yield template
