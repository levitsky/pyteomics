# Licensed under the MIT license:
# http://www.opensource.org/licenses/mit-license.php

import math
import numpy

import modx

# This set of data was downloaded from the NIST website:
# http://www.nist.gov/pml/data/comp.cfm
NIST_MASS = {
    'H': {1: (1.0078250320710, 0.99988570),
          2: (2.01410177784, 0.00011570),
          3: (3.016049277725, 0.0),
          'default': (1.0078250320710, 0.99988570)},
    
    'H+': {1: (1.00727646677, 1.0),
           'default': (1.00727646677, 1.0)},

    'e*': {'default': (0.00054857990943, 1.0)}, 

    'C': {12: (12.0000000, 0.98938),
          13: (13.0033548378, 0.01078),
          14: (14.0032419894, 0.0),
          'default': (12.0000000, 0.98938)},

    'N': {14: (14.00307400486, 0.9963620),
          15: (15.00010889827, 0.0036420),
          'default': (14.00307400486, 0.9963620)},

    'O': {16: (15.9949146195616, 0.9975716),
          17: (16.9991317012, 0.000381),
          18: (17.99916107, 0.0020514),
          'default': (15.9949146195616, 0.9975716)},
    
    'S': {32: (31.9720710015, 0.949926),
          33: (32.9714587615, 0.00752),
          34: (33.9678669012, 0.042524),
          36: (35.9670807620, 0.00011),
          'default': (31.9720710015, 0.949926)},
    }

def _make_isotope_string(element_name, isotope_num):
    return '%s[%d]' % element_name, isotope_num

def _parse_isotope_string(label):
    '''Parse a string with an isotope and return the element name and
    the isotope number.

    >>> _parse_isotope_string('C')
    ('C', 'default')
    >>> _parse_isotope_string('C[12]')
    ('C', 12)    
    '''
    if label.endswith(']'):
        isotope_num = int(label[label.find('[')+1:-1])
        element_name = label[:label.find('[')]
    else:
        isotope_num = 'default'
        element_name = label
    return (element_name, isotope_num)


# Initialize std_amino_acid_comp before the Composition class
# description, fill it later.
STD_AMINO_ACID_COMP = {}

class Composition(dict):
    '''A Composition object stores a chemical composition of a
    substance. Basically it is a dict object, which keys are the names
    of chemical elements and values contains integer number of
    corresponding atoms in a substance.

    The main difference from a dict it that Composition objects allows
    adding and subtraction.
    '''
        
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
    
    def _from_parsed_sequence(self, parsed_sequence,
                              amino_acid_comp=STD_AMINO_ACID_COMP):
        self.clear()
        for aa in parsed_sequence:
            for elem, cnt in amino_acid_comp[aa].items():
                self[elem] = self.get(elem, 0) + amino_acid_comp[aa][elem]
    
    def _from_sequence(self, sequence, amino_acid_comp=STD_AMINO_ACID_COMP):
        self.clear()
        parsed_sequence = modx.parse_sequence(sequence,
                                              show_standard_term_groups=True)
        self._from_parsed_sequence(parsed_sequence)

    def _from_formula(self, formula, mass_data=NIST_MASS):
        prev_chem_symbol_start = len(formula)
        i = len(formula) - 1
        while i >= 0:
            if not (formula[i].isdigit() or formula[i] == '-') :
                if i+1 == prev_chem_symbol_start:
                    num_atoms = 1
                else:
                    try:
                        num_atoms = int(formula[i+1:prev_chem_symbol_start])
                    except ValueError:
                        raise Exception(
                            'Chemical formula is not parseable: %s' % formula)
                element_found = False
                for chem_element in mass_data:
                    if formula.endswith(chem_element, 0, i+1):
                        self[chem_element] = (
                            self.get(chem_element, 0) + num_atoms)
                        i -= len(chem_element)
                        prev_chem_symbol_start = i + 1
                        element_found = True
                        break
                if not element_found:
                    raise Exception(
                        'Chemical formula is not parseable: %s' % formula)
            else:
                i -= 1
                continue

    def __init__(self, *args, **kwargs):
        """
        Keyword arguments:
        sequence -- a polypeptide sequence string;
        parsed_sequence -- a list of amino acids;
        formula -- a string with a chemical formula;
        amino_acid_comp -- a dict with amino acids elemental
                            composition (the default value is
                            std_amino_acids_comp);
        mass_data -- a dict with mass data (the default value is
                     NIST_MASS). It is used for formula parsing only. 
        """

        if 'sequence' in kwargs:
            amino_acid_comp = kwargs.get('amino_acid_comp',
                                         STD_AMINO_ACID_COMP)
            self._from_sequence(kwargs['sequence'], amino_acid_comp)
        elif 'parsed_sequence' in kwargs:
            amino_acid_comp = kwargs.get('amino_acid_comp',
                                         STD_AMINO_ACID_COMP)
            self._from_parsed_sequence(kwargs['parsed_sequence'],
                                       amino_acid_comp)
        elif 'formula' in kwargs:
            mass_data = kwargs.get('mass_data', NIST_MASS)
            self._from_formula(kwargs['formula'], mass_data)
        else:
            self.update(dict(*args, **kwargs))
    
STD_AMINO_ACID_COMP.update({
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

STD_ION_COMP = {
    'M':        Composition(formula=''),
    'a':        Composition(formula='H-2O-1' + 'C-1O-1'),
    'a-H2O':    Composition(formula='H-2O-1' + 'C-1O-1' + 'H-2O-1'),
    'a-NH3':    Composition(formula='H-2O-1' + 'C-1O-1' + 'N-1H-3'),
    'a-NH3-CO': Composition(formula='H-2O-1' + 'C-1O-1' + 'N-1H-3C-1O-1'),
    'b':        Composition(formula='H-2O-1'),
    'b-H2O':    Composition(formula='H-2O-1' + 'H-2O-1'),
    'b-NH3':    Composition(formula='H-2O-1' + 'N-1H-3'),
    'c':        Composition(formula='H-2O-1' + 'NH3'),
    'x':        Composition(formula='H-2O-1' + 'CO2'),
    'y':        Composition(formula=''),
    'y-H2O':    Composition(formula='H-2O-1'),
    'y-NH3':    Composition(formula='N-1H-3'),
    'z':        Composition(formula='H-2O-1' + 'ON-1H-1'),
    }

def calculate_mass(**kwargs):
    '''Calculates the monoisotopic mass of a given polypeptide
    sequence string, a parsed sequence, a chemical formula or
    Composition object.
    Note that if a sequence string is supplied then the mass is
    calculated for a polypeptide with standard terminal groups (NH2-
    and -OH).

    An input may be supplied with any of these keyword arguments:
    sequence -- a polypeptide sequence string;
    parsed_sequence -- a list of amino acids and terminal groups;
    formula -- a string with a chemical formula;
    composition -- a Composition object.

    Other keyword arguments:
    average -- if True then the average mass is calculated;
    ion_type -- if specified, then the polypeptide is considered to be
                in a form of corresponding ion. Do not forget to
                specify the charge state!                
    charge -- if not 0 then m/z is calculated: the mass is increased
              by the corresponding number of proton mass and divided
              by z.

    mass_data -- a dict with mass data (default NIST_MASS);
    amino_acid_comp -- a dict with the elemental composition of amino
                        acid residues (default STD_AMINO_ACID_COMP);
    ion_comp -- a dict with elemental composition of ion losses
                 (default STD_ION_COMP).
    
    Returns a float with the mass of the polypeptide.
    '''

    mass_data = kwargs.get('mass_data', NIST_MASS)
    ion_comp = kwargs.get('ion_comp', STD_ION_COMP)
    # Make a deep copy of `composition` keyword argument.
    composition = (dict(kwargs['composition']) if 'composition' in kwargs
                   else Composition(**kwargs))

    if 'ion_type' in kwargs:
        composition += ion_comp(kwargs['ion_type'])

    # Calculate mass.
    mass = 0.0
    for element in composition:
        element_name, isotope_num = _parse_isotope_string(element)
        if kwargs.get('average', False):
            for isotope_num in mass_data[element_name]:
                if isotope_num != 'default':
                    mass += (composition[element_name]
                             * mass_data[element_name][isotope_num][0]
                             * mass_data[element_name][isotope_num][1])
        else:
            mass += (composition[element]
                     * mass_data[element_name][isotope_num][0])

    # Calculate m/z if required.
    charge = kwargs.get('charge', 0)
    if charge:
        mass = (mass + mass_data['H+']['default'][0] * charge) / charge
    return mass
 
def most_probable_isotopic_composition():
    """Calculates the most probable isotopic composition of a given
    polypeptide sequence string, a parsed sequence, a chemical formula
    or Composition object. For each element, only two most abundant
    isotopes are considered.
    
    Note that if a sequence string is supplied then the isotopic
    composition is calculated for a polypeptide with standard
    terminal groups (NH2- and -OH).

    An input may be supplied with any of these keyword arguments:
    sequence -- a polypeptide sequence string;
    parsed_sequence -- a list of amino acids and terminal groups;
    formula -- a string with a chemical formula;
    composition -- a Composition object.

    Other keyword arguments:
    dynamic_range -- the ratio of the maximal and minimal abudance of
                     generated isotopic compositions (default 100);
    elements_with_isotopes -- a list of elements to be considered as
                              having isotopic distribution. (by default,
                              every element has isotopic distribution);
    mass_data -- a dict with mass data (default NIST_MASS);
    amino_acid_comp -- a dict with the elemental composition of amino
                        acid residues (default STD_AMINO_ACID_COMP);
    ion_comp -- a dict with elemental composition of ion losses
                 (default STD_ION_COMP).

    Returns a tuple (Composition object; abundance).
    """

    composition = (dict(kwargs['composition']) if 'composition' in kwargs
                   else Composition(**kwargs))

    # Removing isotopes from the composition.
    for element in composition:
        element_name, isotope_num = _parse_isotope_string(element)
        if isotope_num != 'default':
            composition[element_name] = (composition.get(element_name, 0)
                                         + composition.pop(element))

    mass_data = kwargs.get('mass_data', NIST_MASS)
    elements_with_isotopes = kwargs.get('elements_with_isotopes', None)
    isotopic_composition = Composition()
    
    for element in composition:
        if not elements_with_isotopes or element in elements_with_isotopes:
            # Take the two most abundant isotopes.
            first_iso, second_iso = sorted(
                [(i[0], i[1][1]) for i in mass_data[element].iteritems()
                 if i[0] != 'default'])[:2]
            
            # Write the number of isotopes of the most abundant type.
            first_iso_str = _make_isotope_string(element, first_iso[0])
            isotopic_composition[first_iso_str] = math.floor(
                (composition[element] + 1) * first_iso[1])

            # Write the number of the second isotopes.
            second_iso_str = _make_isotope_string(element, second_iso[0])
            isotopic_composition[second_iso_str] = (
                composition[element] - isotopic_composition[first_iso_str])

    return isotopic_composition

def factorial_stirling(n):
    return math.sqrt(2 * math.pi / n) * (n / math.e) ** float(n)

def factorial_s3(n):
    """Calculates factorial using Stieltjes' third-order approximation.
    http://www.luschny.de/math/factorial/approx/SimpleCases.html
    The accuracy is ~1e-7 even for n=2.
    """
    N = n + 1.0
    return (factorial_stirling(N)
            * math.exp((1/12.0)/(N+(1/30.0)/(N+(53.0/210.0)/N))))

def binomial_coeff(n, k):
    return factorial_s3(n) / factorial_s3(k) / factorial_s3(n-k)

def binomial_dist(k, n, p):
    return binomial_coeff(n, k) * (p ** k) * (1.0 - p) ** (n - k)

def isotopic_composition_abundance(**kwargs):
    """Calculate the relative abundance of given isotopic composition
    of a molecule.
    
    An input may be supplied with any of these keyword arguments:
    formula -- a string with a chemical formula;
    composition -- a Composition object.

    Other keyword arguments:
    mass_data -- a dict with mass data (default NIST_MASS);

    Returns the relative abundance of given isotopic composition.
    """

    composition = (dict(kwargs['composition']) if 'composition' in kwargs
                   else Composition(**kwargs))

    isotopic_composition = {}

    # Check if there are default and non-default isotopes of the same
    # element and rearrange the elements.
    for element in composition:
        element_name, isotope_num = _parse_isotope_string(element)
        if ((element_name in isotopic_composition) and
        (isotope_num == 'default' or
         'default' in isotopic_composition[element_name])):
            raise Exception('Either fully specify isotopic composition of '
                            '%s or do not specify it at all.' % element_name)
        else:
            if element_name not in isotopic_composition:
                isotopic_composition[element_name] = {}
            isotopic_composition[element_name][isotope_num] = (
                composition[element])

    # Calculate relative abundance.
    mass_data = kwargs.get('mass_data', NIST_MASS)
    relative_abundance = 1.0
    for element_name, isotope_dict in isotopic_composition.items():
        relative_abundance *= factorial_s3(sum(isotope_dict.values()))
        for isotope_num, isotope_content in isotope_dict.items():
            relative_abundance /= factorial_s3(isotope_content)            
            if isotope_num != 'default':
                relative_abundance *= (
                    mass_data[element_name][isotope_num][1]
                    ** isotope_content)
            
    return relative_abundance
    
# def isotopic_distribution(**kwargs):
#     """Generates the isotopic distribution of a given polypeptide
#     sequence string, a parsed sequence, a chemical formula or
#     Composition object.
    
#     Note that if a sequence string is supplied then the isotopic
#     distribution is calculated for a polypeptide with standard
#     terminal groups (NH2- and -OH).

#     An input may be supplied with any of these keyword arguments:
#     sequence -- a polypeptide sequence string;
#     parsed_sequence -- a list of amino acids and terminal groups;
#     formula -- a string with a chemical formula or a Composition object;
#     composition -- a Composition object.

#     Other keyword arguments:
#     dynamic_range -- the ratio of the maximal and minimal abudance of
#                      generated isotopic compositions (default 100);
#     elements_with_isotopes -- a list of elements to be considered as
#                               having isotopic distribution. (default
#                               is None, i.e. every element has isotopic
#                               distribution)
#     mass_data -- a dict with mass data (default NIST_MASS);
#     amino_acid_comp -- a dict with the elemental composition of amino
#                         acid residues (default STD_AMINO_ACID_COMP);
#     ion_comp -- a dict with elemental composition of ion losses
#                  (default STD_ION_COMP).

#     Returns a list of tuples (Composition object; abundance).
#     """

STD_AMINO_ACID_MASS = {
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

def fast_mass(sequence, ion_type=None, charge=None, **kwargs):
    """Calculate monoisotopic mass of an ion using the fast
    algorithm. May be used only if all amino acids are presented in
    one-letter code.

    Keyword arguments:
    sequence -- a polypeptide sequence string;
    ion_type -- if specified, then the polypeptide is considered to be
                in a form of corresponding ion. Do not forget to
                specify the charge state!                
    charge -- if not 0 then m/z is calculated: the mass is increased
              by the corresponding number of proton mass and divided
              by z.

    Other keyword arguments:
    mass_data -- a dict with mass data (default NIST_MASS);
    amino_acid_mass -- a dict with monoisotopic mass of amino
                        acid residues (default STD_AMINO_ACID_MASS);
    ion_comp -- a dict with elemental composition of ion losses
                 (default STD_ION_COMP).
    """
    amino_acid_mass = kwargs.get('amino_acid_mass', STD_AMINO_ACID_MASS)
    mass = sum([amino_acid_mass[i] for i in sequence])

    mass_data = kwargs.get('mass_data', NIST_MASS)
    mass += mass_data['H']['default'][0] * 2.0 + mass_data['O']['default'][0]

    if ion_type:
        ion_comp = kwargs.get('ion_comp', STD_ION_COMP)
        mass += sum(
            [mass_data[element]['default'][0] * ion_comp[ion_type][element]
             for element in ion_comp[ion_type]])
        
    if charge:
        mass = (mass + mass_data['H+']['default'][0] * charge) / charge

    return mass
 
