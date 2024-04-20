'''
proforma - Proteoform and Peptidoform Notation
==============================================

ProForma is a notation for defining modified amino acid sequences using
a set of controlled vocabularies, as well as encoding uncertain or partial
information about localization. See `ProForma specification <https://www.psidev.info/proforma>`_
for more up-to-date information.

For more details, see the :mod:`pyteomics.proforma` online.
'''

import re
import warnings
from collections import deque, namedtuple
from functools import partial
from array import array as _array

try:
    from enum import Enum
except ImportError:
    # Python 2 doesn't have a builtin Enum type
    Enum = object

from .mass import Composition, std_aa_mass, Unimod, nist_mass, calculate_mass, std_ion_comp, mass_charge_ratio
from .auxiliary import PyteomicsError, BasicComposition
from .auxiliary.utils import add_metaclass

try:
    import numpy as np
except ImportError:
    np = None

try:
    from psims.controlled_vocabulary.controlled_vocabulary import (load_psimod, load_xlmod, load_gno, obo_cache, load_unimod)
    _has_psims = True
except ImportError:
    def _needs_psims(name):
        raise ImportError("Loading %s requires the `psims` library. To access it, please install `psims`" % name)

    load_psimod = partial(_needs_psims, 'PSIMOD')
    load_xlmod = partial(_needs_psims, 'XLMOD')
    load_gno = partial(_needs_psims, 'GNO')
    load_unimod = partial(_needs_psims, 'UNIMOD')
    obo_cache = None
    _has_psims = False

_WATER_MASS = calculate_mass(formula="H2O")

std_aa_mass = std_aa_mass.copy()
std_aa_mass['X'] = 0

element_symbols = set(nist_mass)
element_symbols.remove("e*")
element_symbols.add('e')


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
    group_placeholder = 999


class ModificationTagStyle(Enum):
    Unset = 0
    ShortId = 1
    LongId = 2
    ShortName = 3
    LongName = 4


_sentinel = object()


class ModificationMassNotFoundError(ProFormaError):
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
    __slots__ = ("type", "value", "extra", "group_id")

    prefix_name = None
    short_prefix = None
    prefix_map = {}

    def __init__(self, type, value, extra=None, group_id=None):
        self.type = type
        self.value = value
        self.extra = extra
        self.group_id = group_id

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

    def find_tag_type(self, tag_type):
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
    def parse(cls, buffer):
        return process_tag_tokens(buffer)


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
        return str(self.value)


class ModificationResolver(object):
    def __init__(self, name, **kwargs):
        self.name = name.lower()
        self.symbol = self.name[0]
        self._database = None
        self._cache = {}

    def clear_cache(self):
        """Clear the modification definition cache"""
        self._cache.clear()

    def enable_caching(self, flag=True):
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

    def parse_identifier(self, identifier):
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

    def _resolve_impl(self, name=None, id=None, **kwargs):
        raise NotImplementedError()

    def resolve(self, name=None, id=None, **kwargs):
        if self._cache is None:
            return self._resolve_impl(name, id, **kwargs)
        cache_key = (name, id, frozenset(kwargs.items()))
        if cache_key in self._cache:
            return self._cache[cache_key].copy()
        value = self._resolve_impl(name, id, **kwargs)
        self._cache[cache_key] = value
        return  value.copy()

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
        try:
            mass = float(defn.DiffMono)
        except (KeyError, TypeError, ValueError):
            raise ModificationMassNotFoundError("Could not resolve the mass of %r from %r" % ((name, id), defn))
        if defn.DiffFormula is not None:
            composition = Composition()
            diff_formula_tokens = defn.DiffFormula.strip().split(" ")
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
            "name":term.name,
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


class ModificationBase(TagBase):
    '''A base class for all modification tags with marked prefixes.

    While :class:`ModificationBase` is hashable, its equality testing
    brings in additional tag-related information. For pure modification
    identity comparison, use :attr:`key` to get a :class:`ModificationToken`
    free of these concerns..
    '''

    _tag_type = None
    __slots__ = ('_definition', 'style')

    def __init__(self, value, extra=None, group_id=None, style=None):
        if style is None:
            style = ModificationTagStyle.Unset
        super(ModificationBase, self).__init__(
            self._tag_type, value, extra, group_id)
        self._definition = None
        self.style = style

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
    def key(self):
        '''Get a safe-to-hash-and-compare :class:`ModificationToken`
        representing this modification without tag-like properties.

        Returns
        --------
        ModificationToken
        '''
        return ModificationToken(self.value, self.id, self.provider, self.__class__)

    @property
    def definition(self):
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
        -------float
        '''
        return self.definition['mass']

    @property
    def composition(self):
        '''The chemical composition shift this modification applies'''
        return self.definition.get('composition')

    @property
    def id(self):
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

    def _format_main(self):
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
    __slots__ = ('_significant_figures', )

    prefix_name = "Obs"

    def __init__(self, value, extra=None, group_id=None):
        if isinstance(value, str):
            sigfigs = len(value.split('.')[-1].rstrip('0'))
        else:
            sigfigs = 4
        self._significant_figures = sigfigs
        super(MassModification, self).__init__(
            TagTypeEnum.massmod, float(value), extra, group_id)

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
    def key(self):
        '''Get a safe-to-hash-and-compare :class:`ModificationToken`
        representing this modification without tag-like properties.

        Returns
        --------
        ModificationToken
        '''
        return ModificationToken(self.value, self.id, self.provider, self.__class__)

    @property
    def mass(self):
        return self.value

    def __eq__(self, other):
        if isinstance(other, ModificationToken):
            return other == self
        return super(MassModification, self).__eq__(other)

    def __hash__(self):
        return hash((self.id, self.provider))


class FormulaModification(ModificationBase):
    prefix_name = "Formula"

    isotope_pattern = re.compile(r'\[(?P<isotope>\d+)(?P<element>[A-Z][a-z]*)(?P<quantity>[\-+]?\d+)\]')
    _tag_type = TagTypeEnum.formula

    def _normalize_isotope_notation(self, match):
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

    def resolve(self):
        normalized = self.value.replace(' ', '')
        # If there is a [ character in the formula, we know there are isotopes which
        # need to be normalized.
        if '[' in normalized:
            normalized = self.isotope_pattern.sub(self._normalize_isotope_notation, normalized)
        composition = Composition(formula=normalized)
        return {
            "mass": composition.mass(),
            "composition": composition,
            "name": self.value
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

    def __init__(self, name, id, provider, source_cls):
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


def split_tags(tokens):
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


def find_prefix(tokens):
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


def process_marker(tokens):
    '''Process a marker, which is a tag whose value starts with #.

    Parameters
    ----------
    tokens: list
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


def process_tag_tokens(tokens):
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
                    extras.append(GenericModification(''.join(value)))
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

    def is_valid(self, aa, n_term, c_term):
        if (n_term and self.n_term) or (c_term and self.c_term):
            if (self.aa and aa == self.aa) or self.aa is None:
                return True
            return False
        return self.aa == aa or self.aa is None


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

    def __init__(self, modification_tag, targets=None):
        self.modification_tag = modification_tag
        self.targets = targets
        self._validate_targets()

    def is_valid(self, aa, n_term, c_term):
        return any(target.is_valid(aa, n_term, c_term) for target in self.targets)

    def _validate_targets(self):
        validated_targets = []
        if self.targets is None:
            self.targets = []
        elif not isinstance(self.targets, list):
            self.targets = [self.targets]
        for target in self.targets:
            if target in VALID_AA:
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
    '''Define a fixed isotope that is applied globally to all amino acids.

    Attributes
    ----------
    isotope: str
        The stable isotope string, of the form [<isotope-number>]<element> or a special
        isotopoform's name.
    '''
    __slots__ = ('isotope', )

    def __init__(self, isotope):
        self.isotope = isotope

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

    def __init__(self, start, end=None, tags=None, ambiguous=False):
        self.start = start
        self.end = end
        self.tags = tags
        self.ambiguous = ambiguous

    def __eq__(self, other):
        if other is None:
            return False
        return self.start == other.start and self.end == other.end and self.tags == other.tags

    def __ne__(self, other):
        return not self == other

    def __str__(self):
        return "({self.start}-{self.end}){self.tags!r}".format(self=self)

    def __repr__(self):
        return "{self.__class__.__name__}({self.start}, {self.end}, {self.tags})".format(self=self)

    def as_slice(self):
        return slice(self.start, self.end)

    def contains(self, i):
        return self.start <= i < self.end

    def __contains__(self, i):
        return self.contains(i)

    def copy(self):
        return self.__class__(self.start, self.end, self.tags)

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


class ChargeState(object):
    '''Describes the charge and adduct types of the structure.

    Attributes
    ----------
    charge : int
        The total charge state as a signed number.
    adducts : list[str]
        Each charge carrier associated with the molecule.
    '''
    __slots__ = ("charge", "adducts")

    def __init__(self, charge, adducts=None):
        if adducts is None:
            adducts = []
        self.charge = charge
        self.adducts = adducts

    def __str__(self):
        tokens = [str(self.charge)]
        if self.adducts:
            tokens.append("[")
            tokens.append(','.join(str(adduct) for adduct in self.adducts))
            tokens.append("]")
        return ''.join(tokens)

    def __repr__(self):
        template = "{self.__class__.__name__}({self.charge}, {self.adducts})"
        return template.format(self=self)


class TokenBuffer(object):
    '''A token buffer that wraps the accumulation and reset logic
    of a list of :class:`str` objects.

    Implements a subset of the Sequence protocol.

    Attributes
    ----------
    buffer: list
        The list of tokens accumulated since the last parsing.
    '''
    def __init__(self, initial=None):
        self.buffer = list(initial or [])
        self.boundaries = []

    def append(self, c):
        '''Append a new character to the buffer.

        Parameters
        ----------
        c: str
            The character appended
        '''
        self.buffer.append(c)

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

    def tokenize(self):
        i = 0
        pieces = []
        for k in self.boundaries + [len(self)]:
            piece = self.buffer[i:k]
            i = k
            pieces.append(piece)
        return pieces

    def _transform(self, value):
        return value

    def process(self):
        if self.boundaries:
            value = [self._transform(v) for v in self.tokenize()]
        else:
            value = self._transform(self.buffer)
        self.reset()
        return value

    def bound(self):
        k = len(self)
        self.boundaries.append(k)
        return k

    def __call__(self):
        return self.process()


class NumberParser(TokenBuffer):
    '''A buffer which accumulates tokens until it is asked to parse them into
    :class:`int` instances.
    '''

    def _transform(self, value):
        return int(''.join(value))


class StringParser(TokenBuffer):
    '''A buffer which accumulates tokens until it is asked to parse them into
    :class:`str` instances.
    '''

    def _transform(self, value):
        return ''.join(value)


class TagParser(TokenBuffer):
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
DONE = ParserStateEnum.done

VALID_AA = set("QWERTYIPASDFGHKLCVNMXUOJZB")
TERMINAL_SPEC_CHARS = set('N-term') | set('C-term') | set("ncT: ")

def parse(sequence):
    '''Tokenize a ProForma sequence into a sequence of amino acid+tag positions, and a
    mapping of sequence-spanning modifiers.

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
                    "Error In State {state}, unexpected {c} found at index {i}".format(**locals()), i, state)
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
                adduct_buffer = StringParser()
            else:
                raise ProFormaError(
                    "Error In State {state}, unexpected {c} found at index {i}".format(**locals()), i, state)
        elif state == ADDUCT_START:
            if c.isdigit() or c in "+-" or c in element_symbols:
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


def to_proforma(sequence, n_term=None, c_term=None, unlocalized_modifications=None,
                labile_modifications=None, fixed_modifications=None, intervals=None,
                isotopes=None, charge_state=None, group_ids=None):
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
    primary = deque()
    for aa, tags in sequence:
        if not tags:
            primary.append(str(aa))
        else:
            primary.append(str(aa) + ''.join(['[{0!s}]'.format(t) for t in tags]))
    if intervals:
        for iv in sorted(intervals, key=lambda x: x.start):
            if iv.ambiguous:
                primary[iv.start] = '(?' + primary[iv.start]
            else:
                primary[iv.start] = '(' + primary[iv.start]

            terminator = '{0!s})'.format(primary[iv.end - 1])
            if iv.tags:
                terminator += ''.join('[{!s}]'.format(t) for t in iv.tags)
            primary[iv.end - 1] = terminator
    if n_term:
        primary.appendleft(''.join("[{!s}]".format(t) for t in n_term) + '-')
    if c_term:
        primary.append('-' + ''.join("[{!s}]".format(t) for t in c_term))
    if charge_state:
        primary.append("/{!s}".format(charge_state))
    if labile_modifications:
        primary.extendleft(['{{{!s}}}'.format(m) for m in labile_modifications])
    if unlocalized_modifications:
        primary.appendleft("?")
        primary.extendleft(['[{!s}]'.format(m) for m in unlocalized_modifications])
    if isotopes:
        primary.extendleft(['{!s}'.format(m) for m in isotopes])
    if fixed_modifications:
        primary.extendleft(['{!s}'.format(m) for m in fixed_modifications])
    return ''.join(primary)


class _ProFormaProperty(object):
    def __init__(self, name):
        self.name = name

    def __get__(self, obj, cls):
        return obj.properties[self.name]

    def __set__(self, obj, value):
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

    def __init__(self, sequence, properties):
        self.sequence = sequence
        self.properties = properties

    isotopes = _ProFormaProperty('isotopes')
    charge_state = _ProFormaProperty('charge_state')

    intervals = _ProFormaProperty('intervals')
    fixed_modifications = _ProFormaProperty('fixed_modifications')
    labile_modifications = _ProFormaProperty('labile_modifications')
    unlocalized_modifications = _ProFormaProperty('unlocalized_modifications')

    n_term = _ProFormaProperty('n_term')
    c_term = _ProFormaProperty('c_term')

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

            return self.__class__(self.sequence[i], props)
        else:
            return self.sequence[i]

    def __eq__(self, other):
        if isinstance(other, str):
            return str(self) == other
        elif other is None:
            return False
        else:
            return self.sequence == other.sequence and self.properties == other.properties

    def __ne__(self, other):
        return not self == other

    @classmethod
    def parse(cls, string):
        '''Parse a ProForma string.

        Parameters
        ----------
        string : str
            The string to parse

        Returns
        -------
        ProForma
        '''
        return cls(*parse(string))

    @property
    def mass(self):
        mass = 0.0

        fixed_modifications = self.properties['fixed_modifications']

        n_term_v = 0
        c_term_v = len(self) - 1
        for i, position in enumerate(self.sequence):
            aa = position[0]
            try:
                mass += std_aa_mass[aa]
            except KeyError:
                warnings.warn("%r does not have an exact mass" % (aa, ))
            n_term = i == n_term_v
            c_term = i == c_term_v
            for rule in fixed_modifications:
                if rule.is_valid(aa, n_term, c_term):
                    mass += rule.modification_tag.mass
            tags = position[1]
            if tags:
                for tag in tags:
                    try:
                        mass += tag.mass
                    except (AttributeError, KeyError):
                        continue
        for mod in self.properties['labile_modifications']:
            mass += mod.mass
        for mod in self.properties['unlocalized_modifications']:
            mass += mod.mass
        if self.properties.get('n_term'):
            for mod in self.properties['n_term']:
                try:
                    mass += mod.mass
                except (AttributeError, KeyError):
                    continue
        mass += calculate_mass(formula="H")
        if self.properties.get('c_term'):
            for mod in self.properties['c_term']:
                try:
                    mass += mod.mass
                except (AttributeError, KeyError):
                    continue

        mass += calculate_mass(formula="OH")
        for iv in self.properties['intervals']:
            try:
                mass += iv.tag.mass
            except (AttributeError, KeyError):
                continue
        return mass

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
            intervals = sorted(intervals, key=lambda x: x.start)
        intervals = deque(intervals)

        if not include_labile:
            for mod in self.properties['labile_modifications']:
                mass += mod.mass

        if not reverse:
            if self.properties.get('n_term'):
                for mod in self.properties['n_term']:
                    try:
                        mass += mod.mass
                    except (AttributeError, KeyError):
                        continue
        else:
            if self.properties.get('c_term'):
                for mod in self.properties['c_term']:
                    try:
                        mass += mod.mass
                    except (AttributeError, KeyError):
                        continue

        if include_unlocalized:
            for mod in self.properties['unlocalized_modifications']:
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

            aa = position[0]
            try:
                mass += std_aa_mass[aa]
            except KeyError:
                warnings.warn("%r does not have an exact mass" % (aa, ))

            n_term = i == n_term_v
            c_term = i == c_term_v
            for rule in fixed_modifications:
                if rule.is_valid(aa, n_term, c_term):
                    mass += rule.modification_tag.mass

            tags = position[1]
            if tags:
                for tag in tags:
                    try:
                        mass += tag.mass
                    except (AttributeError, KeyError):
                        continue

            while intervals and intervals[0].contains(i):
                iv = intervals.popleft()

                try:
                    mass += iv.tag.mass
                except (AttributeError, KeyError):
                    continue

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
