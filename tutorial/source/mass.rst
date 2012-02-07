Mass and isotopes
=================

Basic functions
---------------

For mass-related problems we propose the :py:mod:`mass` module. It allows
calculation of monoisotopic masses and most probable isotopic states, as well
as estimation of the relative abundance of certain isotopic states for
peptides.

.. seealso:: For further information, please refer to :py:mod:`mass` documentation page.

Mass calculation is the most basic task for this module.

::

    >>> import pyteomics.mass as m
    >>> m.calculate_mass(sequence='PEPTIDE')
    799.359964027207
    >>> calculate_mass(parsed_sequence=parse_sequence('PEPTIDE', show_unmodified_termini=True))
    799.359964027207
    >>> calculate_mass(formula='H2O')
    18.0105646837036

.. warning::
    Don’t forget to set **show_unmodified_termini=True** when parsing a
    sequence if you want to use the result for mass calculation. Otherwise,
    the terminal groups won’t be taken into account!


With :py:func:`calculate_mass` you can calculate masses of modified peptides
as well. Just add the information about the elemental composition of modified
amino acids to the :py:data:`aa_comp` dictionary and supply it to
:py:func:`calculate_mass`. Consider using http://www.unimod.org to get the
composition info.

.. note::

    :py:func:`calculate_mass` can also calculate *m/z* when provided ion type
    and charge.

:py:func:`fast_mass` is a simpler and much faster analog of
:py:func:`calculate_mass`. It works only with unmodified sequences in standard
one-letter IUPAC notation::

    >>> from pyteomics.mass import fast_mass
    >>> fast_mass('PEPTIDE')
    799.3599446837036

.. note::

    Like :py:func:`calculate_mass`, :py:func:`fast_mass` can also calculate *m/z*
    when provided ion type and charge.


Isotopes
--------

If not specified, :py:mod:`pyteomics.mass` assumes that the substances are in
the pure isotopic state. However, you may specify particular isotopic state in
brackets (e.g. O[18], N[15]), while supplying chemical formula to
:py:func:`calculate_mass`::

    >>> calculate_mass(formula='H[2]2O')
    20.0231181752416
    >>> calculate_mass(formula='H[2]HO')
    19.0168414294726

To calculate the relative abundance of a given isotopic state one can use
:py:func:`isotopic_composition_abundance` function. The input can be provided
as a formula or as a Composition object.

::

    >>> from pyteomics.mass import isotopic_composition_abundance, Composition
    >>> isotopic_composition_abundance(formula='H2O')
    1.0
   
.. note::

    If you don’t specify the isotope, the function assumes you don’t care
    about this atom’s isotopic state.

::

    >>> isotopic_composition_abundance(formula='H[1]2O[18]')
    0.002050931076760495
    >>> isotopic_composition_abundance(formula='H[2]H[1]O')
    0.0002313727050147582
    >>> isotopic_composition_abundance(formula='H[2]2O')
    1.3386489999999999e-08
    >>> isotopic_composition_abundance(composition={'H[2]’: 1, ‘H[1]’: 1, ‘O': 1})
    0.0002313727050147582

Finally, you can find the most probable isotopic composition for a substance
with :py:func:`most_probable_isotopic_composition` function. The substance is
specified as a formula of a :py:obj:`Composition` object. Polypeptide
sequences in *modX* notation are also accepted.

::

    >>> from pyteomics.mass import most_probable_isotopic_composition
    >>> most_probable_isotopic_composition(formula='H2SO4')
    {'H[1]': 2.0,  'H[2]': 0.0,  'O[16]': 4.0,  'O[17]': 0.0,  'S[32]': 1.0,  'S[33]': 0.0}
    >>> most_probable_isotopic_composition(formula='C300H602')
    {'C[12]': 297.0, 'C[13]': 3.0, 'H[1]': 602.0, 'H[2]': 0.0}
    >>> most_probable_isotopic_composition(sequence='PEPTIDE'*100)
    {'C[12]': 3364.0,  'C[13]': 36.0,  'H[1]': 5102.0,  'H[2]': 0.0, 'N[14]': 698.0,  'N[15]': 2.0,  'O[16]':  398.0,  'O[17]': 3.0}

Keep in mind that any moment you need information about chemical elements and 
their isotopes the first place to look is the :py:data:`nist_mass` dictionary
defined in :py:mod:`mass.py`. For example::

    >>> from pyteomics.mass import nist_mass
    >>> print nist_mass['C']
    {0: (12.0, 1.0), 12: (12.0, 0.98938), 13: (13.0033548378, 0.01078), 14: (14.0032419894, 0.0)}

The zero key stands for the default (most abundant) isotope. The value for each
key is a *tuple* in the form *(accurate mass, relative abundance)*. The
"relative abundance" in the “default” entry is always 1.0.

::

    >>> print nist_mass['C'][12]
    (12.0, 0.98938)
    >>> print nist_mass['O'][0]
    (15.9949146195616, 1.0)

:py:data:`nist_mass` does not contatin all the periodic table::

    >>> 'Zn' in nist_mass
    False

