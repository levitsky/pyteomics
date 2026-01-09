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
