
.. module:: pyteomics.proforma

proforma - Proteoform and Peptidoform Notation
==============================================

ProForma is a notation for defining modified amino acid sequences using
a set of controlled vocabularies, as well as encoding uncertain or partial
information about localization. See `ProForma specification <https://www.psidev.info/proforma>`_
for more up-to-date information.

Strictly speaking, this implementation supports ProForma v2.

Data Access
-----------

:py:func:`parse` - The primary interface for parsing ProForma strings.

    >>> parse("EM[Oxidation]EVT[#g1(0.01)]S[#g1(0.09)]ES[Phospho#g1(0.90)]PEK")
        ([('E', None),
          ('M', [GenericModification('Oxidation', None, None)]),
          ('E', None),
          ('V', None),
          ('T', [LocalizationMarker(0.01, None, '#g1')]),
          ('S', [LocalizationMarker(0.09, None, '#g1')]),
          ('E', None),
          ('S',
          [GenericModification('Phospho', [LocalizationMarker(0.9, None, '#g1')], '#g1')]),
          ('P', None),
          ('E', None),
          ('K', None)],
         {'n_term': None,
          'c_term': None,
          'unlocalized_modifications': [],
          'labile_modifications': [],
          'fixed_modifications': [],
          'intervals': [],
          'isotopes': [],
          'group_ids': ['#g1']})

:py:func:`to_proforma` - Format a sequence and set of properties as ProForma text.


Classes
-------

:py:class:`ProForma` - An object oriented version of the parsing and formatting code,
coupled with minimal information about mass and position data.

    >>> seq = ProForma.parse("EM[Oxidation]EVT[#g1(0.01)]S[#g1(0.09)]ES[Phospho#g1(0.90)]PEK")
    >>> seq
    ProForma([('E', None), ('M', [GenericModification('Oxidation', None, None)]), ('E', None),
              ('V', None), ('T', [LocalizationMarker(0.01, None, '#g1')]), ('S', [LocalizationMarker(0.09, None, '#g1')]),
              ('E', None), ('S', [GenericModification('Phospho', [LocalizationMarker(0.9, None, '#g1')], '#g1')]),
              ('P', None), ('E', None), ('K', None)],
              {'n_term': None, 'c_term': None, 'unlocalized_modifications': [],
               'labile_modifications': [], 'fixed_modifications': [], 'intervals': [],
               'isotopes': [], 'group_ids': ['#g1'], 'charge_state': None}
            )
    >>> seq.mass
    1360.51054400136
    >>> seq.tags
    [GenericModification('Oxidation', None, None),
     LocalizationMarker(0.01, None, '#g1'),
     LocalizationMarker(0.09, None, '#g1'),
     GenericModification('Phospho', [LocalizationMarker(0.9, None, '#g1')], '#g1')]
    >>> str(seq)
    'EM[Oxidation]EVT[#g1(0.01)]S[#g1(0.09)]ES[Phospho|#g1(0.9)]PEK'

Dependencies
------------

To resolve PSI-MOD, XL-MOD, and GNO identifiers, :mod:`psims` is required. By default,
:mod:`psims` retrieves the most recent version of each controlled vocabulary from the internet, but
includes a fall-back version to use when the network is unavailable. It can also create
an application cache on disk.

CV Disk Caching
~~~~~~~~~~~~~~~
ProForma uses several different controlled vocabularies (CVs) that are each versioned separately.
Internally, the Unimod controlled vocabulary is accessed using :class:`~pyteomics.mass.mass.Unimod`
and all other controlled vocabularies are accessed using :mod:`psims`. Unless otherwise stated,
the machinery will download fresh copies of each CV when first queried.

To avoid this slow operation, you can keep a cached copy of the CV source file on disk and tell
:mod:`pyteomics` and :mod:`psims` where to find them:

.. code-block:: python

    from pyteomics import proforma

    # set the path for Unimod loading via pyteomics
    proforma.set_unimod_path("path/to/unimod.xml")

    # set the cache directory for downloading and reloading OBOs via psims
    proforma.obo_cache.cache_path = "obo/cache/dir/"
    proforma.obo_cache.enabled = True


Compliance Levels
-----------------

1. Base Level Support
Represents the lowest level of compliance, this level involves providing support for:

    - [x] Amino acid sequences
    - [x] Protein modifications using two of the supported CVs/ontologies: Unimod and PSI-MOD.
    - [x] Protein modifications using delta masses (without prefixes)
    - [x] N-terminal, C-terminal and labile modifications.
    - [x] Ambiguity in the modification position, including support for localisation scores.
    - [x] INFO tag.

2. Additional Separate Support
These features are independent from each other:

    - [x] Unusual amino acids (O and U).
    - [x] Ambiguous amino acids (e.g. X, B, Z). This would include support for sequence tags of known mass (using the character X).
    - [x] Protein modifications using delta masses (using prefixes for the different CVs/ontologies).
    - [x] Use of prefixes for Unimod (U:) and PSI-MOD (M:) names.
    - [x] Support for the joint representation of experimental data and its interpretation.

3. Top Down Extensions

    - [ ] Additional CV/ontologies for protein modifications: RESID (the prefix R MUST be used for RESID CV/ontology term names)
    - [x] Chemical formulas (this feature occurs in two places in this list).

4. Cross-Linking Extensions

    - [ ]  Cross-linked peptides (using the XL-MOD CV/ontology, the prefix X MUST be used for XL-MOD CV/ontology term names).

5. Glycan Extensions

    - [x] Additional CV/ontologies for protein modifications: GNO (the prefix G MUST be used for GNO CV/ontology term names)
    - [x] Glycan composition.
    - [x] Chemical formulas (this feature occurs in two places in this list).

6. Spectral Support

    - [x] Charge state and adducts
    - [ ] Chimeric spectra are special cases.
    - [x] Global modifications (e.g., every C is C13).


Functions
---------

.. autofunction:: parse

.. autofunction:: to_proforma

Helpers
~~~~~~~

.. autofunction:: set_unimod_path


High Level Interface
--------------------

.. autoclass:: ProForma


Tag Types
---------

.. autoclass:: TagBase

.. autoclass:: TagTypeEnum


Modification Tags
~~~~~~~~~~~~~~~~~

.. autoclass:: MassModification

.. autoclass:: ModificationBase

.. autoclass:: GenericModification

.. autoclass:: FormulaModification

.. autoclass:: UnimodModification

.. autoclass:: PSIModModification

.. autoclass:: XLMODModification

.. autoclass:: GNOmeModification

.. autoclass:: GlycanModification

.. autoclass:: ModificationToken


Label Tags
~~~~~~~~~~

.. autoclass:: InformationTag

.. autoclass:: PositionLabelTag

.. autoclass:: LocalizationMarker


Supporting Types
----------------

.. autoclass:: ModificationRule

.. autoclass:: StableIsotope

.. autoclass:: TaggedInterval

.. autoclass:: ChargeState

Modification Resolvers
~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: ModificationResolver

.. autoclass:: GenericResolver

.. autoclass:: UnimodResolver

.. autoclass:: PSIModResolver

.. autoclass:: XLMODResolver

.. autoclass:: GNOResolver
