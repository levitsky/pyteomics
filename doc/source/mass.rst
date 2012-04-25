Mass and isotopes
=================

The functions related to mass calculations and isotopic distributions are 
organized into the :py:mod:`pyteomics.mass` module. 

Basic mass calculations
-----------------------

The most common task in mass spectrometry data analysis is to calculate the 
mass of an organic molecule or peptide or *m/z* ratio of an ion. 
The tasks of this kind can be 
performed with the :py:func:`pyteomics.mass.calculate_mass` function. It works with
chemical formulas, polypeptide sequences in modX notation, pre-parsed sequences
and dictionaries of chemical compositions:

.. code-block:: python

    >>> from pyteomics.mass import calculate_mass
    >>> mass.calculate_mass(formula='H2O')
    18.0105646837036

    >>> mass.calculate_mass(formula='C2H5OH')
    46.0418648119876

    >>> mass.calculate_mass(composition={'H':2, 'O':1})
    18.0105646837036

    >>> mass.calculate_mass(sequence='PEPTIDE')
    799.359964027207

    >>> from pyteomics import parser
    >>> ps = parser.parse_sequence('PEPTIDE', show_unmodified_termini=True)
    >>> mass.calculate_mass(parsed_sequence=ps)
    799.359964027207

.. warning::

    Always set **show_unmodified_termini=True** when parsing a
    sequence, if you want to use the result to calculate the mass. Otherwise,
    the mass of the terminal hydrogen and hydroxyl will not be taken into account.

Mass-to-charge ratio of ions
----------------------------

:py:func:`pyteomics.mass.calculate_mass` can be used to calculate the mass/charge ratio of 
peptide ions and ionized fragments. To do that, simply supply the type of the 
peptide ionized fragment and its charge:

.. code-block:: python

    >>> from pyteomics import mass
    >>> mass.calculate_mass(sequence='PEPTIDE', ion_type='M', charge=2)
    400.6872584803735

    >>> mass.calculate_mass(sequence='PEP', ion_type='b', charge=1)
    324.15539725264904

    >>> mass.calculate_mass(sequence='TIDE', ion_type='y', charge=1)
    477.219119708098
   
Mass of modified peptides
-------------------------

With :py:func:`pyteomics.mass.calculate_mass` you can calculate masses of modified peptides
as well. For the function to recognize the modified residue, you need to add the 
information about its elemental composition to the :py:data:`pyteomics.mass.std_aa_comp` 
dictionary used in the calculations by default.

.. code-block:: python

    >>> from pyteomics import mass
    >>> mass.std_aa_comp['pT'] = mass.Composition(
    ...    {'C': 4, 'H': 8, 'N': 1, 'O': 5, 'P': 1})
    >>> mass.calculate_mass(sequence='PEPpTIDE')
    879.3262945499629

To add information about modified amino acids to a user-defined `aa_comp` *dict*
one can either add the composition info for a specific modified residue or just
for a modification:

.. code-block:: python

    >>> from pyteomics import mass
    >>> aa_comp = dict(mass.std_aa_comp)
    >>> aa_comp['p'] = mass.Composition('HPO3')
    >>> mass.calculate_mass('pT', aa_comp=aa_comp)
    199.02457367493957

In this example we call :py:func:`calculate_mass` with a positional
(non-keyword) argument ('pT'). This feature was added in version
1.2.4. When you provide a non-keyword argument, it will be treated as a sequence;
if it fails, it will be treated as a formula; in case it fails as well, a
:py:class:`PyteomicsError` will be raised.
Note that 'pT' is treated as a sequence here, so default terminal groups are
implied when calculating the composition and mass:

.. code-block:: python

    >>> mass.calculate_mass('pT', aa_comp=aa_comp) == mass.calculate_mass(aa_comp['p']) + mass.calculate_mass(aa_comp['T']) + mass.calculate_mass('H2O')
    True

You can create a specific entry for a modified amino acid to override the
modification on a specific residue:

.. code-block:: python

    >>> aa_comp['pT'] = mass.Composition({'N': 2})
    >>> mass.Composition('pT', aa_comp=aa_comp)
    {'H': 2, 'O': 1, 'N': 2}
    >>> mass.Composition('pS', aa_comp=aa_comp)
    {'H': 8, 'C': 3, 'N': 1, 'O': 6, 'P': 1}

`Unimod database <http://www.unimod.org>`_ is an 
excellent resource for the information on the chemical compositions of 
known protein modifications.

Chemical compositions
---------------------

Some problems in organic mass spectrometry deal with molecules made by 
addition or subtraction of standard chemical 'building blocks'. 
In :py:mod:`pyteomics.mass` there are two ways to approach these problems.

* There is a :py:class:`pyteomics.mass.Composition` class intended to store chemical formulas.
  :py:class:`pyteomics.mass.Composition` objects are dicts that can be added or subtracted
  from one another.

  .. code-block:: python

     >>> from pyteomics import mass
     >>> p = mass.Composition(formula='HO3P') # Phosphate group 
     {'H': 1, 'O': 3, 'P': 1}
     >>> print mass.std_aa_comp['T']
     {'C': 4, 'H': 7, 'N': 1, 'O': 2}
     >>> print p + mass.std_aa_comp['T']
     {'C': 4, 'H': 8, 'N': 1, 'O': 5, 'P': 1}

  The values of :py:data:`pyteomics.mass.std_aa_comp` are :py:class:`pyteomics.mass.Composition` objects.

* All functions that accept a **formula** keyword argument sum and 
  subtract numbers following the same atom in the formula:

  .. code-block:: python

     >>> from pyteomics import mass
     >>> mass.calculate_mass(formula='C2H6') # Ethane
     30.046950192426 
     >>> mass.calculate_mass(formula='C2H6H-2') # Ethylene
     28.031300128284002

Faster mass calculations
------------------------

While :py:func:`pyteomics.mass.calculate_mass` has flexible and convenient interface, it may be 
too slow for large-scale calculations. There is an optimized and simplified 
version of this function named :py:func:`pyteomics.mass.fast_mass`. It works only with 
unmodified sequences in standard one-letter IUPAC notation. Like 
:py:func:`pyteomics.mass.calculate_mass`, :py:func:`pyteomics.mass.fast_mass` can calculate *m/z* when
provided with ion type and charge.

.. code-block:: python

    >>> from pyteomicss import mass
    >>> mass.fast_mass('PEPTIDE')
    799.3599446837036

Isotopes
--------

If not specified, :py:mod:`pyteomics.mass` assumes that the substances are in
the pure isotopic state. However, you may specify particular isotopic state in
brackets (e.g. O[18], N[15]) in a chemical formula. An element with unspecified 
isotopic state is assumed to have the mass of the most stable isotope and
abundance of 100%.

.. code-block:: python 

    >>> mass.calculate_mass(formula='H[2]2O') # Heavy water
    20.0231181752416
    >>> mass.calculate_mass(formula='H[2]HO') # Semiheavy water
    19.0168414294726

:py:func:`pyteomics.mass.isotopic_composition_abundance` function calculates the relative 
abundance of a given isotopic state of a molecule. The input can be provided
as a formula or as a Composition/dict. 

.. code-block:: python 

    >>> from pyteomics import mass
    >>> mass.isotopic_composition_abundance(formula='H2O') # Water with an unspecified isotopic state
    1.0
    >>> mass.isotopic_composition_abundance(formula='H[2]2O') # Heavy water
    1.3386489999999999e-08 
    >>> mass.isotopic_composition_abundance(formula='H[2]H[1]O') # Semiheavy water
    0.0002313727050147582
    >>> mass.isotopic_composition_abundance(composition={'H[2]’: 1, ‘H[1]’: 1, ‘O': 1}) # Semiheavy water
    0.0002313727050147582
    >>> mass.isotopic_composition_abundance(formula='H[2]2O[18]') # Heavy-hydrogen heavy-oxygen water
    2.7461045585999998e-11

.. warning::

    You cannot mix specified and unspecified states of the same element in one 
    formula in :py:func:`pyteomics.mass.isotopic_composition_abundance` due to ambiguity.

    .. code-block:: python

        >>> mass.isotopic_composition_abundance(formula='H[2]HO') 
        ...
        PyteomicsError: Pyteomics error, message: 'Please specify the isotopic states of all atoms of H or do not specify them at all.'
  
Finally, you can find the most probable isotopic composition for a substance
with :py:func:`pyteomics.mass.most_probable_isotopic_composition` function. The substance is
specified as a formula, a :py:class:`pyteomics.mass.Composition` object or a modX sequence string.

.. code-block:: python

    >>> from pyteomics import mass
    >>> mass.most_probable_isotopic_composition(formula='H2SO4')
    {'H[1]': 2.0,  'H[2]': 0.0,  'O[16]': 4.0,  'O[17]': 0.0,  'S[32]': 1.0,  'S[33]': 0.0}
    >>> mass.most_probable_isotopic_composition(formula='C300H602')
    {'C[12]': 297.0, 'C[13]': 3.0, 'H[1]': 602.0, 'H[2]': 0.0}
    >>> mass.most_probable_isotopic_composition(sequence='PEPTIDE'*100)
    {'C[12]': 3364.0,  'C[13]': 36.0,  'H[1]': 5102.0,  'H[2]': 0.0, 'N[14]': 698.0,  'N[15]': 2.0,  'O[16]':  398.0,  'O[17]': 3.0}

The information about chemical elements, their isotopes and relative abundances
is stored in the :py:data:`pyteomics.mass.nist_mass` dictionary defined in :py:mod:`pyteomics.mass.mass`.

.. code-block:: python

    >>> from pyteomics import mass
    >>> print mass.nist_mass['C']
    {0: (12.0, 1.0), 12: (12.0, 0.98938), 13: (13.0033548378, 0.01078), 14: (14.0032419894, 0.0)}

The zero key stands for the unspecified isotopic state. The data about isotopes 
are stored as tuples *(accurate mass, relative abundance)*.

At the moment, :py:data:`pyteomics.mass.nist_mass` has the data only for the atoms of organic
chemistry, the proton and electron:

.. code-block:: python

    >>> print mass.nist_mass.keys()
    ['H+', 'C', 'P', 'e*', 'H', 'S', 'O', 'N']

