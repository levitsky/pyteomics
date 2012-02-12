"""
mass - molecular masses and isotope distributions
=================================================

Summary
-------

This module defines general functions for mass and isotope abundance
calculations. For most of the functions, a user may define a studied
substance in various formats but all of them would be reduced to the
:py:func:`Composition <Composition.__init__>` object describing its
chemical composition.


Classes
-------

  :py:func:`Composition <Composition.__init__>` - a class storing chemical
  composition of a substance.

Mass calculations
-----------------

  :py:func:`calculate_mass` - a general routine for mass / m/z
  calculation. Can calculate mass for a polypeptide sequence, chemical
  formula or elemental composition. Supplied with an ion type and
  charge, the function would calculate m/z.
  
  :py:func:`fast_mass` - a less powerful but much faster function for
  polypeptide mass calculation.
  
Isotopic abundances
-------------------

  :py:func:`isotopic_composition_abundance` - calculate the relative
  abundance of a given isotopic composition.

  :py:func:`most_probable_isotopic_composition` - finds the most
  abundant isotopic composition for a molecule defined by a
  polypeptide sequence, chemical formula or elemental composition.

Data
----

  :py:data:`nist_mass` - a dict with exact masses of the most abundant
  isotopes.

  :py:data:`std_aa_comp` - a dict with the elemental compositions
  of the standard twenty amino acid residues.

  :py:data:`std_ion_comp` - a dict with the relative elemental
  compositions of the standard peptide fragment ions.

  :py:data:`std_aa_mass` - a dict with the monoisotopic masses
  of the standard twenty amino acid residues.

.. ipython::
   :suppress:

   In [1]: import pyteomics.mass; from pprint import pprint

-----------------------------------------------------------------------------
"""
# Licensed under the MIT license:
# http://www.opensource.org/licenses/mit-license.php

import math
import numpy
import parser
from auxiliary import PyteomicsError

nist_mass = {
    'H': {1: (1.0078250320710, 0.99988570),
          2: (2.01410177784, 0.00011570),
          3: (3.016049277725, 0.0),
          0: (1.0078250320710, 1.0)},
    
    'H+': {1: (1.00727646677, 1.0),
           0: (1.00727646677, 1.0)},

    'e*': {0: (0.00054857990943, 1.0)}, 

    'C': {12: (12.0000000, 0.98938),
          13: (13.0033548378, 0.01078),
          14: (14.0032419894, 0.0),
           0: (12.0000000, 1.0)},

    'N': {14: (14.00307400486, 0.9963620),
          15: (15.00010889827, 0.0036420),
           0: (14.00307400486, 1.0)},

    'O': {16: (15.9949146195616, 0.9975716),
          17: (16.9991317012, 0.000381),
          18: (17.99916107, 0.0020514),
           0: (15.9949146195616, 1.0)},
    
    'P': {31: (30.9737616320, 1.0000),
           0: (30.9737616320, 1.0000)},

    'S': {32: (31.9720710015, 0.949926),
          33: (32.9714587615, 0.00752),
          34: (33.9678669012, 0.042524),
          36: (35.9670807620, 0.00011),
           0: (31.9720710015, 1.0)},
    }
"""
A dict with the exact element masses downloaded from the NIST website:
http://www.nist.gov/pml/data/comp.cfm . There are entries for each
element containing the masses and relative abundances of several
abundant isotopes and a separate entry for undefined isotope with zero
key, mass of the most abundant isotope and 1.0 abundance.

.. ipython::
   
   In [2]: pprint(pyteomics.mass.nist_mass)
"""

def _make_isotope_string(element_name, isotope_num):
    """Form a string label for an isotope."""
    if isotope_num == 0:
        return element_name
    else:
        return '%s[%d]' % (element_name, isotope_num)

def _parse_isotope_string(label):
    """Parse an string with an isotope label and return the element name and
    the isotope number.

    >>> _parse_isotope_string('C')
    ('C', 0)
    >>> _parse_isotope_string('C[12]')
    ('C', 12)    
    """
    if label.endswith(']'):
        isotope_num = int(label[label.find('[')+1:-1])
        element_name = label[:label.find('[')]
    else:
        isotope_num = 0
        element_name = label
    return (element_name, isotope_num)

# Initialize std_aa_comp before the Composition class
# description, fill it later.
std_aa_comp = {}
"""A dictionary with elemental compositions of the twenty standard
amino acid residues and standard H- and -OH terminal groups.

.. ipython::
   
    In [2]: pprint(pyteomics.mass.std_aa_comp)
"""

class Composition(dict):
    """
    A Composition object stores a chemical composition of a
    substance. Basically, it is a dict object, with the names
    of chemical elements as keys and values equal to an integer number of
    atoms of the corresponding element in a substance.

    The main improvement over dict is that Composition objects allow
    adding and subtraction.
    """
        
    def __add__(self, other):
        result = self.copy()
        for elem, cnt in other.items():
            result[elem] = result.get(elem, 0) + cnt
        return result

    def __sub__(self, other):
        result = self.copy()
        for elem, cnt in other.items():
            result[elem] = result.get(elem, 0) - cnt
        return result
    
    def _from_parsed_sequence(self, parsed_sequence, aa_comp=std_aa_comp):
        self.clear()
        for aa in parsed_sequence:
            for elem, cnt in aa_comp[aa].items():
                self[elem] = self.get(elem, 0) + aa_comp[aa][elem]
    
    def _from_sequence(self, sequence, aa_comp=std_aa_comp):
        self.clear()
        parsed_sequence = parser.parse_sequence(
            sequence,
            labels=aa_comp.keys(),
            show_unmodified_termini=True)
        self._from_parsed_sequence(parsed_sequence, aa_comp)
        
    def _from_formula(self, formula, mass_data=nist_mass):
        # Parsing a formula backwards.
        prev_chem_symbol_start = len(formula)
        i = len(formula) - 1
        
        while i >= 0:
            # Read backwards until a non-number character is met.
            if (formula[i].isdigit() or formula[i] == '-'):
                i -= 1
                continue
            
            else:
                # If the number of atoms is omitted than it is 1.
                if i+1 == prev_chem_symbol_start:
                    num_atoms = 1
                else:
                    try:
                        num_atoms = int(formula[i+1:prev_chem_symbol_start])
                    except ValueError:
                        raise PyteomicsError(
                            'Badly-formed number of atoms: %s' % formula)

                # Read isotope number if specified else it is undefined (=0).
                if formula[i] == ']':
                    bra_pos = formula.rfind('[', 0, i)
                    if bra_pos == -1:
                        raise PyteomicsError(
                            'Badly-formed isotope number: %s' % formula)
                    try:
                        isotope_num = int(formula[bra_pos+1:i])
                    except ValueError:
                        raise PyteomicsError(
                            'Badly-formed isotope number: %s' % formula)
                    i = bra_pos - 1
                else:
                    isotope_num = 0

                # Match the element name to the mass_data.
                element_found = False
                for element_name in mass_data:
                    if formula.endswith(element_name, 0, i+1):
                        isotope_string = _make_isotope_string(
                            element_name, isotope_num)
                        self[isotope_string] = (
                            self.get(isotope_string, 0) + num_atoms)
                        i -= len(element_name)
                        prev_chem_symbol_start = i + 1
                        element_found = True
                        break
                    
                if not element_found:
                    raise PyteomicsError(
                        'Unknown chemical element in the formula: %s' %formula)

    def __init__(self, *args, **kwargs):
        """
        A Composition object stores a chemical composition of a
        substance. Basically it is a dict object, in which keys are the names
        of chemical elements and values contain integer numbers of
        corresponding atoms in a substance.
        
        The main improvement over dict is that Composition objects allow
        addition and subtraction.

        A Composition object can be initialized with one of the
        following arguments: formula, sequence or parsed_sequence.
         
        .. warning::
         
            Be cafeful when supplying a list with a parsed sequence. It must be
            obtained with enabled `show_unmodified_termini` option.
        
        Parameters
        ----------
        formula : str, optional
            A string with a chemical formula.
        sequence : str, optional
            A polypeptide sequence string in modX notation.
        parsed_sequence : list of str, optional
            A polypeptide sequence parsed into a list of amino acids.
        aa_comp : dict, optional
            A dict with the elemental composition of the amino acids (the
            default value is std_aa_comp).
        mass_data : dict, optional
            A dict with the masses of the chemical elements (the default
            value is nist_mass). It is used for formulae parsing only. 
        """
        if len(args) == 0 and len(kwargs) == 0:
            pass
        elif 'sequence' in kwargs:
            aa_comp = kwargs.get('aa_comp', std_aa_comp)
            self._from_sequence(kwargs['sequence'], aa_comp)
            
        elif 'parsed_sequence' in kwargs:
            aa_comp = kwargs.get('aa_comp', std_aa_comp)
            self._from_parsed_sequence(kwargs['parsed_sequence'], aa_comp)
            
        elif 'formula' in kwargs:
            mass_data = kwargs.get('mass_data', nist_mass)
            self._from_formula(kwargs['formula'], mass_data)
            
        elif isinstance(args[0], dict):
            mass_data = kwargs.get('mass_data', nist_mass)
            for isotope_string, num_atoms in args[0].iteritems():
                element_name, isotope_num = _parse_isotope_string(
                    isotope_string)

                if element_name in mass_data:
                    # Remove explicitly undefined isotopes (e.g. X[0]).
                    self[_make_isotope_string(element_name, isotope_num)] = (
                        num_atoms)                    
                else:
                    raise PyteomicsError('Unknown chemical element: %s' %
                                         (element_name,))
            
        else:
            raise PyteomicsError('A Composition object must be specified by '
                                 'a sequence, parsed sequence, formula or '
                                 'a dict')
    
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
    'H-':  Composition({'H': 1}),
    '-OH': Composition({'O': 1, 'H': 1}),
    })

std_ion_comp = {
    'M':        Composition(formula=''),
    'a':        Composition(formula='H-2O-1' + 'C-1O-1'),
    'a-H2O':    Composition(formula='H-2O-1' + 'C-1O-1' + 'H-2O-1'),
    'a-NH3':    Composition(formula='H-2O-1' + 'C-1O-1' + 'N-1H-3'),
    'b':        Composition(formula='H-2O-1'),
    'b-H2O':    Composition(formula='H-2O-1' + 'H-2O-1'),
    'b-NH3':    Composition(formula='H-2O-1' + 'N-1H-3'),
    'c':        Composition(formula='H-2O-1' + 'NH3'),
    'c-H2O':    Composition(formula='H-2O-1' + 'NH3' + 'H-2O-1'),
    'c-NH3':    Composition(formula='H-2O-1'),
    'x':        Composition(formula='H-2O-1' + 'CO2'),
    'x-H2O':    Composition(formula='H-2O-1' + 'CO2' + 'H-2O-1'),
    'x-NH3':    Composition(formula='H-2O-1' + 'CO2' + 'N-1H-3'),
    'y':        Composition(formula=''),
    'y-H2O':    Composition(formula='H-2O-1'),
    'y-NH3':    Composition(formula='N-1H-3'),
    'z':        Composition(formula='H-2O-1' + 'ON-1H-1'),
    'z-H2O':    Composition(formula='H-2O-1' + 'ON-1H-1' + 'H-2O-1'),
    'z-NH3':    Composition(formula='H-2O-1' + 'ON-1H-1' + 'N-1H-3'),
    }
"""A dict with relative elemental compositions of the standard peptide
fragment ions. An elemental composition of a fragment ion is calculated as a
difference between the total elemental composition of an ion
and the sum of elemental compositions of its constituting amino acid residues.

.. ipython::
   
   In [2]: pprint(pyteomics.mass.std_ion_comp)
"""

def calculate_mass(**kwargs):
    """Calculates the monoisotopic mass of a polypeptide defined by a 
    sequence string, parsed sequence, chemical formula or
    Composition object.

    One and only one of the following keyword arguments is required:
    **formula**, **sequence**, **parsed_sequence** or **composition**.
    
    Note that if a sequence string is supplied then the mass is
    calculated for a polypeptide with standard terminal groups (NH2-
    and -OH).

    .. warning::

        Be cafeful when supplying a list with a parsed sequence. It must be
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
    average : bool, optional
        If True then the average mass is calculated. Note that mass is not
        averaged for elements with specified isotopes.
    ion_type : str, optional
        If specified, then the polypeptide is considered to be in a form
        of the corresponding ion. Do not forget to specify the charge state!
    charge : int, optional
        If not 0 then m/z is calculated: the mass is increased
        by the corresponding number of proton masses and divided
        by z.
    aa_comp : dict, optional
        A dict with the elemental composition of the amino acids (the
        default value is std_aa_comp).
    mass_data : dict, optional
        A dict with the masses of the chemical elements (the default
        value is nist_mass). 
    ion_comp : dict, optional
        A dict with the relative elemental compositions of peptide ion
        fragments (default is std_ion_comp).
    
    Returns
    -------
        mass : float
    """

    mass_data = kwargs.get('mass_data', nist_mass)
    ion_comp = kwargs.get('ion_comp', std_ion_comp)
    # Make a deep copy of `composition` keyword argument.
    composition = (dict(kwargs['composition'])
                   if 'composition' in kwargs
                   else Composition(**kwargs))

    if 'ion_type' in kwargs:
        composition += ion_comp[kwargs['ion_type']]

    # Get charge.
    charge = 0
    if 'H+' in composition:
        charge = composition['H+']
    if 'charge' in kwargs:
        if charge:
            raise PyteomicsError(
                'Charge is specified both by the number of protons and '
                '\'charge\' in kwargs: %s' % str(kwargs))
        charge = kwargs['charge']
        composition['H+'] = charge

    # Calculate mass.
    mass = 0.0
    for isotope_string in composition:
        element_name, isotope_num = _parse_isotope_string(isotope_string)
        # Calculate average mass if required and the isotope number is
        # not specified.
        if (not isotope_num) and kwargs.get('average', False):
            for isotope in mass_data[element_name]:
                if isotope != 0:
                    mass += (composition[element_name]
                             * mass_data[element_name][isotope][0]
                             * mass_data[element_name][isotope][1])
        else:
            mass += (composition[isotope_string]
                     * mass_data[element_name][isotope_num][0])

    # Calculate m/z if required.
    if charge:
        mass = mass / charge
    return mass
 
def most_probable_isotopic_composition(**kwargs):
    """Calculate the most probable isotopic composition of a peptide
    molecule/ion defined by a sequence string, parsed sequence,
    chemical formula or Composition object.

    Note that if a sequence string is supplied then the isotopic
    composition is calculated for a polypeptide with standard
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
    composition : Composition, optional
        A Composition object with the elemental composition of a substance.
    elements_with_isotopes : list of str
        A list of elements to be considered in isotopic distribution
        (by default, every element has a isotopic distribution).
    aa_comp : dict, optional
        A dict with the elemental composition of the amino acids (the
        default value is std_aa_comp).
    mass_data : dict, optional
        A dict with the masses of the chemical elements (the default
        value is nist_mass). 
    ion_comp : dict, optional
        A dict with the relative elemental compositions of peptide ion
        fragments (default is std_ion_comp).

    Returns
    -------
    out: tuple (Composition, float)
        A tuple with the most probable isotopic composition and its
        relative abundance.
    """

    composition = (dict(kwargs['composition']) if 'composition' in kwargs
                   else Composition(**kwargs))

    # Removing isotopes from the composition.
    for isotope_string in composition:
        element_name, isotope_num = _parse_isotope_string(isotope_string)
        if isotope_num:
            composition[element_name] = (composition.get(element_name, 0)
                                         + composition.pop(isotope_string))

    mass_data = kwargs.get('mass_data', nist_mass)
    elements_with_isotopes = kwargs.get('elements_with_isotopes', None)
    isotopic_composition = Composition()
    
    for element_name in composition:
        if (not elements_with_isotopes
        or (element_name in elements_with_isotopes)):
            # Take the two most abundant isotopes.
            first_iso, second_iso = sorted(
                [(i[0], i[1][1])
                 for i in mass_data[element_name].iteritems() if i[0]])[:2]
            
            # Write the number of isotopes of the most abundant type.
            first_iso_str = _make_isotope_string(element_name, first_iso[0])
            isotopic_composition[first_iso_str] = math.floor(
                (composition[element_name] + 1) * first_iso[1])

            # Write the number of the second isotopes.
            second_iso_str = _make_isotope_string(element_name, second_iso[0])
            isotopic_composition[second_iso_str] = (
                composition[element_name]
                - isotopic_composition[first_iso_str])
        else:
            isotopic_composition[element_name] = composition[element_name]

    return isotopic_composition

def factorial_stirling(n):
    return math.sqrt(2 * math.pi / n) * (n / math.e) ** float(n)

def factorial_s3(n):
    # Calculate factorial using Stieltjes' third-order approximation.
    # http://www.luschny.de/math/factorial/approx/SimpleCases.html
    # The accuracy is ~1e-7 even for n=2. 
    
    N = n + 1.0
    return (factorial_stirling(N)
            * math.exp((1/12.0)/(N+(1/30.0)/(N+(53.0/210.0)/N))))

def binomial_coeff(n, k):
    return factorial_s3(n) / factorial_s3(k) / factorial_s3(n-k)

def binomial_dist(k, n, p):
    return binomial_coeff(n, k) * (p ** k) * (1.0 - p) ** (n - k)

def isotopic_composition_abundance(**kwargs):
    """Calculate the relative abundance of a given isotopic composition
    of a molecule.
    
    Parameters
    ----------
    formula : str, optional
        A string with a chemical formula.
    composition : Composition, optional
        A Composition object with the isotopic composition of a substance.
    mass_data : dict, optional
        A dict with the masses of the chemical elements (the default
        value is nist_mass). 

    Returns
    -------
    relative_abundance : float
        The relative abundance of a given isotopic composition.
    """

    composition = (dict(kwargs['composition'])
                   if 'composition' in kwargs
                   else Composition(**kwargs))

    isotopic_composition = {}

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
            if element_name not in isotopic_composition:
                isotopic_composition[element_name] = {}
            isotopic_composition[element_name][isotope_num] = (
                composition[element])

    # Calculate relative abundance.
    mass_data = kwargs.get('mass_data', nist_mass)
    relative_abundance = 1.0
    for element_name, isotope_dict in isotopic_composition.items():
        relative_abundance *= factorial_s3(sum(isotope_dict.values()))
        for isotope_num, isotope_content in isotope_dict.items():
            relative_abundance /= factorial_s3(isotope_content)            
            if isotope_num:
                relative_abundance *= (
                    mass_data[element_name][isotope_num][1]
                    ** isotope_content)
            
    return relative_abundance

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
    'R': 156.10111,
    'Y': 163.06333,
    'W': 186.07931,
    }
"""A dictionary with monoisotopic masses of the twenty standard
amino acid residues.

.. ipython::
   
    In [2]: pprint(pyteomics.mass.std_aa_mass)
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
        A dict with the masses of the chemical elements (the default
        value is nist_mass). 
    aa_mass : dict, optional
        A dict with the monoisotopic mass of amino acid residues
        (default is std_aa_mass);
    ion_comp : dict, optional
        A dict with the relative elemental compositions of peptide ion
        fragments (default is std_ion_comp).

    Returns
    -------
    mass : float
        Monoisotopic mass or m/z of a peptide molecule/ion.
    """
    aa_mass = kwargs.get('aa_mass', std_aa_mass)
    mass = sum([aa_mass[i] for i in sequence])

    mass_data = kwargs.get('mass_data', nist_mass)
    mass += mass_data['H'][0][0] * 2.0 + mass_data['O'][0][0]

    if ion_type:
        ion_comp = kwargs.get('ion_comp', std_ion_comp)
        mass += sum(
            [mass_data[element][0][0] * ion_comp[ion_type][element]
             for element in ion_comp[ion_type]])
        
    if charge:
        mass = (mass + mass_data['H+'][0][0] * charge) / charge

    return mass
 
