"""
mass - molecular masses and isotope distributions
=================================================

Summary
-------

This module defines general functions for mass and isotope abundance
calculations. For most of the functions, the user can define a given
substance in various formats, but all of them would be reduced to the
:py:func:`Composition <Composition.__init__>` object describing its
chemical composition.


Classes
-------

  :py:func:`Composition <Composition.__init__>` - a class storing chemical
  composition of a substance.

  :py:class:`Unimod` - a class representing a Python interface to the
  `Unimod database <http://unimod.org/>`_
  (see :py:mod:`pyteomics.mass.unimod` for a much more powerful alternative).

Mass calculations
-----------------

  :py:func:`calculate_mass` - a general routine for mass / m/z
  calculation. Can calculate mass for a polypeptide sequence, chemical
  formula or elemental composition. Supplied with an ion type and
  charge, the function would calculate m/z.

  :py:func:`fast_mass` - a less powerful but much faster function for
  polypeptide mass calculation.

  :py:func:`fast_mass2` - a version of `fast_mass` that supports *modX* notation.

Isotopic abundances
-------------------

  :py:func:`isotopic_composition_abundance` - calculate the relative
  abundance of a given isotopic composition.

  :py:func:`most_probable_isotopic_composition` - finds the most
  abundant isotopic composition for a molecule defined by a
  polypeptide sequence, chemical formula or elemental composition.

  :py:func:`isotopologues` - iterate over possible isotopic conposition of a molecule,
  possibly filtered by abundance.

Data
----

  :py:data:`nist_mass` - a dict with exact masses of the most abundant
  isotopes.

  :py:data:`std_aa_comp` - a dict with the elemental compositions
  of the standard twenty amino acid residues, selenocysteine and pyrrolysine.

  :py:data:`std_ion_comp` - a dict with the relative elemental
  compositions of the standard peptide fragment ions.

  :py:data:`std_aa_mass` - a dict with the monoisotopic masses
  of the standard twenty amino acid residues, selenocysteine and pyrrolysine.

-----------------------------------------------------------------------------
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

from __future__ import division
import math
from .. import parser
from ..auxiliary import PyteomicsError, _nist_mass, BasicComposition
from itertools import chain, product, combinations_with_replacement
from collections import defaultdict
try:
    from urllib import urlopen
except ImportError:
    from urllib.request import urlopen
from datetime import datetime
import re
import operator

nist_mass = _nist_mass
"""
A dict with the exact element masses downloaded from the NIST website:
http://www.nist.gov/pml/data/comp.cfm . There are entries for each
element containing the masses and relative abundances of several
abundant isotopes and a separate entry for undefined isotope with zero
key, mass of the most abundant isotope and 1.0 abundance.
"""

def _make_isotope_string(element_name, isotope_num):
    """Form a string label for an isotope."""
    if isotope_num == 0:
        return element_name
    else:
        return '{}[{}]'.format(element_name, isotope_num)

def _parse_isotope_string(label):
    """Parse an string with an isotope label and return the element name and
    the isotope number.

    >>> _parse_isotope_string('C')
    ('C', 0)
    >>> _parse_isotope_string('C[12]')
    ('C', 12)
    """
    element_name, num = re.match(_isotope_string, label).groups()
    isotope_num = int(num) if num else 0
    return element_name, isotope_num

# Initialize std_aa_comp and std_ion_comp before the Composition class
# description, fill it later.
std_aa_comp = {}
"""A dictionary with elemental compositions of the twenty standard
amino acid residues, selenocysteine, pyrrolysine,
and standard H- and -OH terminal groups.
"""

std_ion_comp = {}
"""A dict with relative elemental compositions of the standard peptide
fragment ions. An elemental composition of a fragment ion is calculated as a
difference between the total elemental composition of an ion
and the sum of elemental compositions of its constituting amino acid residues.
"""

_isotope_string = r'^([A-Z][a-z+]*)(?:\[(\d+)\])?$'
_atom = r'([A-Z][a-z+]*)(?:\[(\d+)\])?([+-]?\d+)?'
_formula = r'^({})*$'.format(_atom)

class Composition(BasicComposition):
    """
    A Composition object stores a chemical composition of a
    substance. Basically, it is a dict object, with the names
    of chemical elements as keys and values equal to an integer number of
    atoms of the corresponding element in a substance.

    The main improvement over dict is that Composition objects allow
    adding and subtraction.
    """
    _kw_sources = {'formula', 'sequence', 'parsed_sequence', 'split_sequence', 'composition'}

    def _from_parsed_sequence(self, parsed_sequence, aa_comp):
        self.clear()
        comp = defaultdict(int)
        for aa in parsed_sequence:
            if aa in aa_comp:
                for elem, cnt in aa_comp[aa].items():
                    comp[elem] += cnt
            else:
                try:
                    mod, aa = parser._split_label(aa)
                    for elem, cnt in chain(
                            aa_comp[mod].items(), aa_comp[aa].items()):
                        comp[elem] += cnt

                except (PyteomicsError, KeyError):
                    raise PyteomicsError(
                            'No information for %s in `aa_comp`' % aa)
        self._from_composition(comp)

    def _from_split_sequence(self, split_sequence, aa_comp):
        self.clear()
        comp = defaultdict(int)
        for group in split_sequence:
            i = 0
            while i < len(group):
                for j in range(len(group)+1, -1, -1):
                    try:
                        label = ''.join(group[i:j])
                        for elem, cnt in aa_comp[label].items():
                            comp[elem] += cnt
                    except KeyError:
                        continue
                    else:
                        i = j
                        break
                if j == 0:
                    raise PyteomicsError("Invalid group starting from position %d: %s" % (i+1, group))
        self._from_composition(comp)

    def _from_sequence(self, sequence, aa_comp):
        parsed_sequence = parser.parse(
            sequence,
            labels=aa_comp,
            show_unmodified_termini=True)
        self._from_parsed_sequence(parsed_sequence, aa_comp)

    def _from_formula(self, formula, mass_data):
        if not re.match(_formula, formula):
            raise PyteomicsError('Invalid formula: ' + formula)
        for elem, isotope, number in re.findall(_atom, formula):
            if not elem in mass_data:
                raise PyteomicsError('Unknown chemical element: ' + elem)
            self[_make_isotope_string(elem, int(isotope) if isotope else 0)
                    ] += int(number) if number else 1

    def _from_composition(self, comp):
        for isotope_string, num_atoms in comp.items():
            element_name, isotope_num = _parse_isotope_string(
                isotope_string)

            # Remove explicitly undefined isotopes (e.g. X[0]).
            self[_make_isotope_string(element_name, isotope_num)] = num_atoms

    def __init__(self, *args, **kwargs):
        """
        A Composition object stores a chemical composition of a
        substance. Basically it is a dict object, in which keys are the names
        of chemical elements and values contain integer numbers of
        corresponding atoms in a substance.

        The main improvement over dict is that Composition objects allow
        addition and subtraction.

        A Composition object can be initialized with one of the
        following arguments: formula, sequence, parsed_sequence or
        split_sequence.

        If none of these are specified, the constructor will look at the first
        positional argument and try to build the object from it. Without
        positional arguments, a Composition will be constructed directly from
        keyword arguments.

        If there's an ambiguity, i.e. the argument is both a valid sequence
        and a formula (such as 'HCN'), it will be treated as a sequence. You
        need to provide the 'formula' keyword to override this.

        .. warning::

            Be careful when supplying a list with a parsed sequence or a split
            sequence as a keyword argument. It must be
            obtained with enabled `show_unmodified_termini` option.
            When supplying it as a positional argument, the option doesn't
            matter, because the positional argument is always converted to
            a sequence prior to any processing.

        Parameters
        ----------
        formula : str, optional
            A string with a chemical formula. All elements must be present in
            `mass_data`.
        sequence : str, optional
            A polypeptide sequence string in modX notation.
        parsed_sequence : list of str, optional
            A polypeptide sequence parsed into a list of amino acids.
        split_sequence : list of tuples of str, optional
            A polypeptyde sequence parsed into a list of tuples
            (as returned be :py:func:`pyteomics.parser.parse` with
            ``split=True``).
        aa_comp : dict, optional
            A dict with the elemental composition of the amino acids (the
            default value is std_aa_comp).
        mass_data : dict, optional
            A dict with the masses of chemical elements (the default
            value is :py:data:`nist_mass`). It is used for formulae parsing only.
        charge : int, optional
            If not 0 then additional protons are added to the composition.
        ion_comp : dict, optional
            A dict with the relative elemental compositions of peptide ion
            fragments (default is :py:data:`std_ion_comp`).
        ion_type : str, optional
            If specified, then the polypeptide is considered to be in the form
            of the corresponding ion. Do not forget to specify the charge state!
        """
        defaultdict.__init__(self, int)

        aa_comp = kwargs.get('aa_comp', std_aa_comp)
        mass_data = kwargs.get('mass_data', nist_mass)


        kw_given = self._kw_sources.intersection(kwargs)
        if len(kw_given) > 1:
            raise PyteomicsError('Only one of {} can be specified!\n'
                    'Given: {}'.format(', '.join(self._kw_sources),
                        ', '.join(kw_given)))
        elif kw_given:
            kwa = kw_given.pop()
            getattr(self, '_from_' + kwa)(kwargs[kwa],
                    mass_data if kwa == 'formula' else aa_comp)

        # can't build from kwargs
        elif args:
            if isinstance(args[0], dict):
                self._from_composition(args[0])
            elif isinstance(args[0], str):
                try:
                    self._from_sequence(args[0], aa_comp)
                except PyteomicsError:
                    try:
                        self._from_formula(args[0], mass_data)
                    except PyteomicsError:
                        raise PyteomicsError(
                                'Could not create a Composition object from '
                                'string: "{}": not a valid sequence or '
                                'formula'.format(args[0]))
            else:
                try:
                    self._from_sequence(parser.tostring(args[0], True), aa_comp)
                except:
                    raise PyteomicsError('Could not create a Composition object'
                            ' from `{}`. A Composition object must be '
                            'specified by sequence, parsed or split sequence,'
                            ' formula or dict.'.format(args[0]))
        else:
            self._from_composition(kwargs)

        ion_comp = kwargs.get('ion_comp', std_ion_comp)
        if 'ion_type' in kwargs:
            self += ion_comp[kwargs['ion_type']]

        # Get charge
        charge = self['H+']
        if 'charge' in kwargs:
            if charge:
                raise PyteomicsError(
                    'Charge is specified both by the number of protons and '
                    '`charge` in kwargs')
            charge = kwargs['charge']
            self['H+'] = charge

    def mass(self, **kwargs):
        """Calculate the mass or *m/z* of a :py:class:`Composition`.

        Parameters
        ----------
        average : bool, optional
            If :py:const:`True` then the average mass is calculated. Note that mass
            is not averaged for elements with specified isotopes. Default is
            :py:const:`False`.
        charge : int, optional
            If not 0 then m/z is calculated: the mass is increased
            by the corresponding number of proton masses and divided
            by `charge`.
        mass_data : dict, optional
            A dict with the masses of the chemical elements (the default
            value is :py:data:`nist_mass`).
        ion_comp : dict, optional
            A dict with the relative elemental compositions of peptide ion
            fragments (default is :py:data:`std_ion_comp`).
        ion_type : str, optional
            If specified, then the polypeptide is considered to be in the form
            of the corresponding ion. Do not forget to specify the charge state!

        Returns
        -------
        mass : float
        """
        composition = self
        mass_data = kwargs.get('mass_data', nist_mass)

        # Calculate mass
        mass = 0.0
        average = kwargs.get('average', False)
        for isotope_string, amount in composition.items():
            element_name, isotope_num = _parse_isotope_string(isotope_string)
            # Calculate average mass if required and the isotope number is
            # not specified.
            if (not isotope_num) and average:
                for isotope, data in mass_data[element_name].items():
                    if isotope:
                        mass += (amount * data[0] * data[1])
            else:
                mass += (amount * mass_data[element_name][isotope_num][0])

        # Calculate m/z if required
        charge = kwargs.get('charge', composition['H+'])
        if charge:
            if not composition['H+']:
                mass += mass_data['H+'][0][0] * charge
            mass /= charge
        return mass

std_aa_comp.update({
    'A':   Composition({'H': 5, 'C': 3, 'O': 1, 'N': 1}),
    'C':   Composition({'H': 5, 'C': 3, 'S': 1, 'O': 1, 'N': 1}),
    'D':   Composition({'H': 5, 'C': 4, 'O': 3, 'N': 1}),
    'E':   Composition({'H': 7, 'C': 5, 'O': 3, 'N': 1}),
    'F':   Composition({'H': 9, 'C': 9, 'O': 1, 'N': 1}),
    'G':   Composition({'H': 3, 'C': 2, 'O': 1, 'N': 1}),
    'H':   Composition({'H': 7, 'C': 6, 'N': 3, 'O': 1}),
    'I':   Composition({'H': 11, 'C': 6, 'O': 1, 'N': 1}),
    'K':   Composition({'H': 12, 'C': 6, 'N': 2, 'O': 1}),
    'L':   Composition({'H': 11, 'C': 6, 'O': 1, 'N': 1}),
    'M':   Composition({'H': 9, 'C': 5, 'S': 1, 'O': 1, 'N': 1}),
    'N':   Composition({'H': 6, 'C': 4, 'O': 2, 'N': 2}),
    'P':   Composition({'H': 7, 'C': 5, 'O': 1, 'N': 1}),
    'Q':   Composition({'H': 8, 'C': 5, 'O': 2, 'N': 2}),
    'R':   Composition({'H': 12, 'C': 6, 'N': 4, 'O': 1}),
    'S':   Composition({'H': 5, 'C': 3, 'O': 2, 'N': 1}),
    'T':   Composition({'H': 7, 'C': 4, 'O': 2, 'N': 1}),
    'V':   Composition({'H': 9, 'C': 5, 'O': 1, 'N': 1}),
    'W':   Composition({'C': 11, 'H': 10, 'N': 2, 'O': 1}),
    'Y':   Composition({'H': 9, 'C': 9, 'O': 2, 'N': 1}),
    'U':   Composition({'H': 5, 'C': 3, 'O': 1, 'N': 1, 'Se' : 1}),
    'O':   Composition({'H': 19, 'C': 12, 'O': 2, 'N': 3}),
    'H-':  Composition({'H': 1}),
    '-OH': Composition({'O': 1, 'H': 1}),
    })

std_ion_comp.update({
    'M':        Composition(formula=''),
    'M-H2O':    Composition(formula='H-2O-1'),
    'M-NH3':    Composition(formula='N-1H-3'),
    'a':        Composition(formula='H-2O-1' + 'C-1O-1'),
    'a-H2O':    Composition(formula='H-2O-1' + 'C-1O-1' + 'H-2O-1'),
    'a-NH3':    Composition(formula='H-2O-1' + 'C-1O-1' + 'N-1H-3'),
    'b':        Composition(formula='H-2O-1'),
    'b-H2O':    Composition(formula='H-2O-1' + 'H-2O-1'),
    'b-NH3':    Composition(formula='H-2O-1' + 'N-1H-3'),
    'c':        Composition(formula='H-2O-1' + 'NH3'),
    'c-1':      Composition(formula='H-2O-1' + 'NH3' + 'H-1'),
    'c-dot':    Composition(formula='H-2O-1' + 'NH3' + 'H1'),
    'c+1':      Composition(formula='H-2O-1' + 'NH3' + 'H1'),
    'c+2':      Composition(formula='H-2O-1' + 'NH3' + 'H2'),
    'c-H2O':    Composition(formula='H-2O-1' + 'NH3' + 'H-2O-1'),
    'c-NH3':    Composition(formula='H-2O-1'),
    'x':        Composition(formula='H-2O-1' + 'CO2'),
    'x-H2O':    Composition(formula='H-2O-1' + 'CO2' + 'H-2O-1'),
    'x-NH3':    Composition(formula='H-2O-1' + 'CO2' + 'N-1H-3'),
    'y':        Composition(formula=''),
    'y-H2O':    Composition(formula='H-2O-1'),
    'y-NH3':    Composition(formula='N-1H-3'),
    'z':        Composition(formula='H-2O-1' + 'ON-1H-1'),
    'z-dot':    Composition(formula='H-2O-1' + 'ON-1'),
    'z+1':      Composition(formula='H-2O-1' + 'ON-1H1'),
    'z+2':      Composition(formula='H-2O-1' + 'ON-1H2'),
    'z+3':      Composition(formula='H-2O-1' + 'ON-1H3'),
    'z-H2O':    Composition(formula='H-2O-1' + 'ON-1H-1' + 'H-2O-1'),
    'z-NH3':    Composition(formula='H-2O-1' + 'ON-1H-1' + 'N-1H-3'),
    })


def calculate_mass(*args, **kwargs):
    """Calculates the monoisotopic mass of a polypeptide defined by a
    sequence string, parsed sequence, chemical formula or
    Composition object.

    One or none of the following keyword arguments is required:
    **formula**, **sequence**, **parsed_sequence**, **split_sequence**
    or **composition**.
    All arguments given are used to create a :py:class:`Composition` object,
    unless an existing one is passed as a keyword argument.

    Note that if a sequence string is supplied and terminal groups are not
    explicitly shown, then the mass is calculated for a polypeptide with
    standard terminal groups (NH2- and -OH).

    .. warning::

        Be careful when supplying a list with a parsed sequence. It must be
        obtained with enabled `show_unmodified_termini` option.

    Parameters
    ----------
    formula : str, optional
        A string with a chemical formula.
    sequence : str, optional
        A polypeptide sequence string in modX notation.
    parsed_sequence : list of str, optional
        A polypeptide sequence parsed into a list of amino acids.
    composition : Composition, optional
        A Composition object with the elemental composition of a substance.
    aa_comp : dict, optional
        A dict with the elemental composition of the amino acids (the
        default value is std_aa_comp).
    average : bool, optional
        If :py:const:`True` then the average mass is calculated. Note that mass
        is not averaged for elements with specified isotopes. Default is
        :py:const:`False`.
    charge : int, optional
        If not 0 then m/z is calculated: the mass is increased
        by the corresponding number of proton masses and divided
        by `charge`.
    mass_data : dict, optional
        A dict with the masses of the chemical elements (the default
        value is :py:data:`nist_mass`).
    ion_comp : dict, optional
        A dict with the relative elemental compositions of peptide ion
        fragments (default is :py:data:`std_ion_comp`).
    ion_type : str, optional
        If specified, then the polypeptide is considered to be in the form
        of the corresponding ion. Do not forget to specify the charge state!

    Returns
    -------
    mass : float
    """
    # Make a copy of `composition` keyword argument.
    composition = (Composition(kwargs['composition'])
                   if 'composition' in kwargs
                   else Composition(*args, **kwargs))
    return composition.mass(**kwargs)

def most_probable_isotopic_composition(*args, **kwargs):
    """Calculate the most probable isotopic composition of a peptide
    molecule/ion defined by a sequence string, parsed sequence,
    chemical formula or :py:class:`Composition` object.

    Note that if a sequence string without terminal groups is supplied then the
    isotopic composition is calculated for a polypeptide with standard
    terminal groups (H- and -OH).

    For each element, only two most abundant isotopes are considered.

    Parameters
    ----------
    formula : str, optional
        A string with a chemical formula.
    sequence : str, optional
        A polypeptide sequence string in modX notation.
    parsed_sequence : list of str, optional
        A polypeptide sequence parsed into a list of amino acids.
    composition : :py:class:`Composition`, optional
        A :py:class:`Composition` object with the elemental composition of a
        substance.
    elements_with_isotopes : list of str
        A list of elements to be considered in isotopic distribution
        (by default, every element has a isotopic distribution).
    aa_comp : dict, optional
        A dict with the elemental composition of the amino acids (the
        default value is :py:data:`std_aa_comp`).
    mass_data : dict, optional
        A dict with the masses of chemical elements (the default
        value is :py:data:`nist_mass`).
    ion_comp : dict, optional
        A dict with the relative elemental compositions of peptide ion
        fragments (default is :py:data:`std_ion_comp`).

    Returns
    -------
    out: tuple (Composition, float)
        A tuple with the most probable isotopic composition and its
        relative abundance.
    """

    composition = (dict(kwargs['composition']) if 'composition' in kwargs
                   else Composition(*args, **kwargs))

    # Removing isotopes from the composition.
    for isotope_string in composition:
        element_name, isotope_num = _parse_isotope_string(isotope_string)
        if isotope_num:
            composition[element_name] += composition.pop(isotope_string)

    mass_data = kwargs.get('mass_data', nist_mass)
    elements_with_isotopes = kwargs.get('elements_with_isotopes')
    isotopic_composition = Composition()

    for element_name in composition:
        if (not elements_with_isotopes
        or (element_name in elements_with_isotopes)):
            # Take the two most abundant isotopes.
            first_iso, second_iso = sorted(
                [(i[0], i[1][1])
                     for i in mass_data[element_name].items() if i[0]],
                key=lambda x: -x[1])[:2]

            # Write the number of isotopes of the most abundant type.
            first_iso_str = _make_isotope_string(element_name, first_iso[0])
            isotopic_composition[first_iso_str] = int(math.ceil(
                composition[element_name])) * first_iso[1]

            # Write the number of the second isotopes.
            second_iso_str = _make_isotope_string(element_name, second_iso[0])
            isotopic_composition[second_iso_str] = (
                composition[element_name]
                - isotopic_composition[first_iso_str])
        else:
            isotopic_composition[element_name] = composition[element_name]

    return (isotopic_composition,
            isotopic_composition_abundance(
                composition=isotopic_composition,
                mass_data=mass_data))

def isotopic_composition_abundance(*args, **kwargs):
    """Calculate the relative abundance of a given isotopic composition
    of a molecule.

    Parameters
    ----------
    formula : str, optional
        A string with a chemical formula.
    composition : Composition, optional
        A Composition object with the isotopic composition of a substance.
    mass_data : dict, optional
        A dict with the masses of chemical elements (the default
        value is :py:data:`nist_mass`).

    Returns
    -------
    relative_abundance : float
        The relative abundance of a given isotopic composition.
    """

    composition = (Composition(kwargs['composition'])
                   if 'composition' in kwargs
                   else Composition(*args, **kwargs))

    isotopic_composition = defaultdict(dict)

    # Check if there are default and non-default isotopes of the same
    # element and rearrange the elements.
    for element in composition:
        element_name, isotope_num = _parse_isotope_string(element)

        # If there is already an entry for this element and either it
        # contains a default isotope or newly added isotope is default
        # then raise an exception.
        if ((element_name in isotopic_composition)
             and (isotope_num == 0
                  or 0 in isotopic_composition[element_name])):
            raise PyteomicsError(
                'Please specify the isotopic states of all atoms of '
                '%s or do not specify them at all.' % element_name)
        else:
            isotopic_composition[element_name][isotope_num] = (
                composition[element])

    # Calculate relative abundance.
    mass_data = kwargs.get('mass_data', nist_mass)
    num1, num2, denom = 1, 1, 1
    for element_name, isotope_dict in isotopic_composition.items():
        num1 *= math.factorial(sum(isotope_dict.values()))
        for isotope_num, isotope_content in isotope_dict.items():
            denom *= math.factorial(isotope_content)
            if isotope_num:
                num2 *= (mass_data[element_name][isotope_num][1]
                        ** isotope_content)

    return num2 * (num1 / denom)

def isotopologues(*args, **kwargs):
    """Iterate over possible isotopic states of a molecule.
    The molecule can be defined by formula, sequence, parsed sequence, or composition.
    The space of possible isotopic compositions is restrained by parameters
    ``elements_with_isotopes``, ``isotope_threshold``, ``overall_threshold``.

    Parameters
    ----------
    formula : str, optional
        A string with a chemical formula.
    sequence : str, optional
        A polypeptide sequence string in modX notation.
    parsed_sequence : list of str, optional
        A polypeptide sequence parsed into a list of amino acids.
    composition : :py:class:`Composition`, optional
        A :py:class:`Composition` object with the elemental composition of a
        substance.
    report_abundance : bool, optional
        If :py:const:`True`, the output will contain 2-tuples: `(composition, abundance)`.
        Otherwise, only compositions are yielded. Default is :py:const:`False`.
    elements_with_isotopes : container of str, optional
        A set of elements to be considered in isotopic distribution
        (by default, every element has an isotopic distribution).
    isotope_threshold : float, optional
        The threshold abundance of a specific isotope to be considered.
        Default is :py:const:`5e-4`.
    overall_threshold : float, optional
        The threshold abundance of the calculateed isotopic composition.
        Default is :py:const:`0`.
    aa_comp : dict, optional
        A dict with the elemental composition of the amino acids (the
        default value is :py:data:`std_aa_comp`).
    mass_data : dict, optional
        A dict with the masses of chemical elements (the default
        value is :py:data:`nist_mass`).

    Returns
    -------
    out : iterator
        Iterator over possible isotopic compositions.
    """
    iso_threshold = kwargs.pop('isotope_threshold', 5e-4)
    overall_threshold = kwargs.pop('overall_threshold', 0.0)
    mass_data = kwargs.get('mass_data', nist_mass)
    elements_with_isotopes = kwargs.get('elements_with_isotopes')
    report_abundance = kwargs.get('report_abundance', False)
    composition = Composition(kwargs['composition']) if 'composition' in kwargs else Composition(*args, **kwargs)
    other_kw = kwargs.copy()
    for k in Composition._kw_sources:
        other_kw.pop(k, None)

    dict_elem_isotopes = {}
    for element in composition:
        if elements_with_isotopes is None or element in elements_with_isotopes:
            element_name, isotope_num = _parse_isotope_string(element)
            isotopes = {k: v for k, v in mass_data[element_name].items() if k != 0 and v[1] >= iso_threshold}
            list_isotopes = [_make_isotope_string(element_name, k) for k in isotopes]
            dict_elem_isotopes[element] = list_isotopes
        else:
            dict_elem_isotopes[element] = [element]
    all_isotoplogues = []
    for element, list_isotopes in dict_elem_isotopes.items():
        n = composition[element]
        list_comb_element_n = []
        for elementXn in combinations_with_replacement(list_isotopes, n):
            list_comb_element_n.append(elementXn)
        all_isotoplogues.append(list_comb_element_n)

    for isotopologue in product(*all_isotoplogues):
        ic = Composition(formula=''.join(atom for el in isotopologue for atom in el), **other_kw)
        if report_abundance or overall_threshold > 0.0:
            abundance = isotopic_composition_abundance(composition=ic, **other_kw)
            if abundance > overall_threshold:
                if report_abundance:
                    yield (ic, abundance)
                else:
                    yield ic
        else:
            yield ic

std_aa_mass = {
    'G': 57.02146,
    'A': 71.03711,
    'S': 87.03203,
    'P': 97.05276,
    'V': 99.06841,
    'T': 101.04768,
    'C': 103.00919,
    'L': 113.08406,
    'I': 113.08406,
    'N': 114.04293,
    'D': 115.02694,
    'Q': 128.05858,
    'K': 128.09496,
    'E': 129.04259,
    'M': 131.04049,
    'H': 137.05891,
    'F': 147.06841,
    'U': 150.95364,
    'R': 156.10111,
    'Y': 163.06333,
    'W': 186.07931,
    'O': 255.15829,
    }
"""A dictionary with monoisotopic masses of the twenty standard
amino acid residues, selenocysteine and pyrrolysine.
"""

def fast_mass(sequence, ion_type=None, charge=None, **kwargs):
    """Calculate monoisotopic mass of an ion using the fast
    algorithm. May be used only if amino acid residues are presented in
    one-letter code.

    Parameters
    ----------
    sequence : str
        A polypeptide sequence string.
    ion_type : str, optional
        If specified, then the polypeptide is considered to be
        in a form of corresponding ion. Do not forget to
        specify the charge state!
    charge : int, optional
        If not 0 then m/z is calculated: the mass is increased
        by the corresponding number of proton masses and divided
        by z.
    mass_data : dict, optional
        A dict with the masses of chemical elements (the default
        value is :py:data:`nist_mass`).
    aa_mass : dict, optional
        A dict with the monoisotopic mass of amino acid residues
        (default is std_aa_mass);
    ion_comp : dict, optional
        A dict with the relative elemental compositions of peptide ion
        fragments (default is :py:data:`std_ion_comp`).

    Returns
    -------
    mass : float
        Monoisotopic mass or m/z of a peptide molecule/ion.
    """
    aa_mass = kwargs.get('aa_mass', std_aa_mass)
    try:
        mass = sum(aa_mass[i] for i in sequence)
    except KeyError as e:
        raise PyteomicsError('No mass data for residue: ' + e.args[0])

    mass_data = kwargs.get('mass_data', nist_mass)
    mass += mass_data['H'][0][0] * 2 + mass_data['O'][0][0]

    if ion_type:
        try:
            icomp = kwargs.get('ion_comp', std_ion_comp)[ion_type]
        except KeyError:
            raise PyteomicsError('Unknown ion type: {}'.format(ion_type))

        mass += sum(mass_data[element][0][0] * num
             for element, num in icomp.items())

    if charge:
        mass = (mass + mass_data['H+'][0][0] * charge) / charge

    return mass

def fast_mass2(sequence, ion_type=None, charge=None, **kwargs):
    """Calculate monoisotopic mass of an ion using the fast
    algorithm. *modX* notation is fully supported.

    Parameters
    ----------
    sequence : str
        A polypeptide sequence string.
    ion_type : str, optional
        If specified, then the polypeptide is considered to be
        in a form of corresponding ion. Do not forget to
        specify the charge state!
    charge : int, optional
        If not 0 then m/z is calculated: the mass is increased
        by the corresponding number of proton masses and divided
        by z.
    mass_data : dict, optional
        A dict with the masses of chemical elements (the default
        value is :py:data:`nist_mass`).
    aa_mass : dict, optional
        A dict with the monoisotopic mass of amino acid residues
        (default is std_aa_mass);
    ion_comp : dict, optional
        A dict with the relative elemental compositions of peptide ion
        fragments (default is :py:data:`std_ion_comp`).

    Returns
    -------
    mass : float
        Monoisotopic mass or m/z of a peptide molecule/ion.
    """
    aa_mass = kwargs.get('aa_mass', std_aa_mass)
    mass_data = kwargs.get('mass_data', nist_mass)
    aa_mass.setdefault('H-', mass_data['H'][0][0])
    aa_mass.setdefault('-OH', mass_data['H'][0][0] + mass_data['O'][0][0])
    try:
        comp = parser.amino_acid_composition(sequence,
                show_unmodified_termini=True,
                allow_unknown_modifications=True,
                labels=aa_mass)
    except PyteomicsError:
        raise PyteomicsError('Mass not specified for label(s): {}'.format(
            ', '.join(set(parser.parse(sequence)).difference(aa_mass))))

    try:
        mass = 0
        for aa, num in comp.items():
            if aa in aa_mass:
                mass += aa_mass[aa] * num
            else:
                mod, X = parser._split_label(aa)
                mass += (aa_mass[mod] + aa_mass[X]) * num
    except KeyError as e:
        raise PyteomicsError(
                'Unspecified mass for modification: "{}"'.format(e.args[0]))

    if ion_type:
        try:
            icomp = kwargs.get('ion_comp', std_ion_comp)[ion_type]
        except KeyError:
            raise PyteomicsError('Unknown ion type: {}'.format(ion_type))

        mass += sum(mass_data[element][0][0] * num
             for element, num in icomp.items())

    if charge:
        mass = (mass + mass_data['H+'][0][0] * charge) / charge

    return mass

class Unimod():
    """A class for Unimod database of modifications.
    The list of all modifications can be retrieved via `mods` attribute.
    Methods for convenient searching are `by_title` and `by_name`.
    For more elaborate filtering, iterate manually over the list.

    .. note::
        See :py:mod:`pyteomics.mass.unimod` for a new alternative class with
        more features.
    """

    def __init__(self, source='http://www.unimod.org/xml/unimod.xml'):
        """Create a database and fill it from XML file retrieved from `source`.

        Parameters
        ----------

        source : str or file, optional
            A file-like object or a URL to read from. Don't forget the ``'file://'``
            prefix when pointing to local files.
        """
        from lxml import etree
        from ..xml import _local_name
        def process_mod(mod):
            d = mod.attrib
            new_d = {}
            for key in ('date_time_modified', 'date_time_posted'):
                new_d[key] = datetime.strptime(d.pop(key),
                    '%Y-%m-%d %H:%M:%S')
            comp = Composition()
            for delta in self._xpath('delta', mod): # executed 1 time
                for key in ('avge_mass', 'mono_mass'):
                    new_d[key] = float(delta.attrib.pop(key))
                for elem in self._xpath('element', delta):
                    e_d = elem.attrib
                    amount = int(e_d.pop('number'))
                    label = e_d.pop('symbol')
                    isotope, symbol = re.match('^(\d*)(\D+)$', label).groups()
                    if not isotope: isotope = 0
                    else: isotope = int(isotope)
                    comp += Composition(
                            formula = _make_isotope_string(symbol, isotope),
                            mass_data = self._massdata) * amount
            new_d['composition'] = comp
            new_d['record_id'] = int(d.pop('record_id'))
            new_d['approved'] = (d.pop('approved') == '1')
            new_d.update(d)
            spec = []
            for sp in self._xpath('specificity', mod):
                sp_d = sp.attrib
                sp_new_d = {}
                sp_new_d['hidden'] = (sp_d.pop('hidden') == '1')
                sp_new_d['spec_group'] = int(sp_d.pop('spec_group'))
                sp_new_d.update(sp_d)
                notes = []
                for note in self._xpath('*', sp):
                    if note.text and note.text.strip():
                        notes.append(note.text.strip())
                if notes:
                    sp_new_d['note'] = '\n'.join(notes)
                spec.append(sp_new_d)
            new_d['specificity'] = spec

            alt_names = []
            for alt_name in self._xpath('alt_name', mod):
                alt_names.append(alt_name.text)
            if alt_names:
                new_d['alt_names'] = alt_names

            refs = []
            for ref in self._xpath('xref', mod):
                ref_d = {}
                for sub in ref.iterchildren():
                    ref_d[_local_name(sub)] = sub.text
                for key in ('text', 'source', 'url'):
                    if key not in ref_d:
                        ref_d[key] = None
                refs.append(ref_d)
            new_d['refs'] = refs
            return new_d

        if isinstance(source, str):
            self._tree = etree.parse(urlopen(source))
        else:
            self._tree = etree.parse(source)
        self._massdata = self._mass_data()
        self._mods = []
        for mod in self._xpath('/unimod/modifications/mod'):
            self._mods.append(process_mod(mod))

    def _xpath(self, path, element=None):
        from ..xml import xpath
        if element is None:
            return xpath(self._tree, path, 'umod')
        return xpath(element, path, 'umod')

    def _mass_data(self):
        massdata = defaultdict(dict)
        elements = [x.attrib for x in self._xpath('/unimod/elements/elem')]
        avg = {}
        for elem in elements:
            i, label = re.match('^(\d*)(\D+)$', elem['title']).groups()
            if not i:
                iso = 0
            else:
                iso = int(i)
            massdata[label][iso] = (float(elem['mono_mass']), float(iso == 0))
            if not iso:
                avg[label] = float(elem['avge_mass'])
        for elem, isotopes in massdata.items():
            isotopes[int(round(isotopes[0][0]))] = isotopes[0]
            if len(isotopes) == 3:
                m1, m2 = (x[1][0] for x in sorted(isotopes.items())[1:])
                m_avg = avg[elem]
                a = (m2 - m_avg) / (m2 - m1)
                b = (m_avg - m1) / (m2 - m1)
                for state, abundance in zip(sorted(isotopes)[1:], (a, b)):
                    isotopes[state] = (isotopes[state][0], abundance)
        return massdata

    @property
    def mods(self):
        """Get the list of Unimod modifications"""
        return self._mods

    @property
    def mass_data(self):
        """Get element mass data extracted from the database"""
        return self._massdata

    def by_title(self, title, strict=True):
        """Search modifications by title. If a single modification is found,
        it is returned. Otherwise, a list will be returned.

        Parameters
        ----------
        title : str
            The modification title.
        strict : bool, optional
            If :py:const:`False`, the search will return all modifications
            whose title **contains** `title`, otherwise equality is required.
            :py:const:`True` by default.

        Returns
        -------
        out : dict or list
            A single modification or a list of modifications.
        """
        f = {True: operator.eq, False: operator.contains}
        func = f[strict]
        result = [m for m in self._mods if func(m['title'], title)]
        if len(result) == 1:
            return result[0]
        return result

    def by_name(self, name, strict=True):
        """Search modifications by name. If a single modification is found,
        it is returned. Otherwise, a list will be returned.

        Parameters
        ----------
        name : str
            The full name of the modification(s).
        strict : bool, optional
            If :py:const:`False`, the search will return all modifications
            whose full name **contains** `title`, otherwise equality is
            required. :py:const:`True` by default.

        Returns
        -------
        out : dict or list
            A single modification or a list of modifications.
        """
        f = {True: operator.eq, False: operator.contains}
        func = f[strict]
        result = [m for m in self._mods if func(m['full_name'], name)]
        if len(result) == 1:
            return result[0]
        return result
