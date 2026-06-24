The *ProForma* standard and implementation
------------------------------------------

`ProForma <https://www.psidev.info/proforma>`_ is a standard for representing proteoforms and peptidoforms,
developed by the `PSI <https://www.psidev.info/>`_ (Proteomics Standards Initiative).
It provides a structured way to represent peptide sequences a wide variety of modifications and uncertainties.

Pyteomics supports ProForma v2.0. The core functions and classes related to ProForma support are located in the :py:doc:`api/proforma`,
see there for more information.

Basic usage
~~~~~~~~~~~

The **ProForma** parser is object-oriented, with a primary class :py:class:`ProForma` representing a parsed ProForma sequence.
To instantiate a :py:class:`ProForma` object, use the class method :py:meth:`ProForma.parse`::

    .. code-block:: python

        >>> seq = ProForma.parse("EM[Oxidation]EVT[Phospho]SES[Phospho]PEK")
        >>> seq
        ProForma([('E', None), ('M', [GenericModification('Oxidation', None, None)]), ('E', None), ('V', None), ('T', [GenericModification('Phospho', None, None)]), ('S', None), ('E', None), ('S', [GenericModification('Phospho', None, None)]), ('P', None), ('E', None), ('K', None)], {'n_term': None, 'c_term': None, 'unlocalized_modifications': [], 'labile_modifications': [], 'fixed_modifications': [], 'intervals': [], 'isotopes': [], 'group_ids': [], 'charge_state': None})

        >>> seq.mass
        1440.47687500136

        >>> seq.composition()
        Composition({'H': 86, 'C': 51, 'O': 30, 'N': 12, 'S': 1, 'P': 2})

Chimeric spectra
~~~~~~~~~~~~~~~~

Top-level ``+`` in ProForma is treated as a chimeric separator only when
``chimeric=True`` is passed. The return value is then a list of parsed
components:

.. code-block:: python

    >>> forms = ProForma.parse("<[Carbamidomethyl]@C>AC+CC", chimeric=True)
    >>> len(forms)
    2
    >>> [str(form) for form in forms]
    ['<[Carbamidomethyl]@C>AC', '<[Carbamidomethyl]@C>CC']

Fixed modification rules, isotope labels, and peptidoform names are shared
across all chimeric components.

Other APIs such as mass calculation, fragment series generation, and spectrum
annotation operate on one peptidoform at a time. Use the parsed components
individually:

.. code-block:: python

    >>> from pyteomics import mass
    >>> masses = [mass.calculate_mass(proforma=str(form)) for form in forms]
    >>> fragments = [mass.fragment_series(str(form)) for form in forms]
