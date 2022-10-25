.. module:: pyteomics.mass.unimod

unimod - interface to the Unimod database
=========================================

This module provides an interface to the relational Unimod database.
The main class is :py:class:`Unimod`, which provides an identical interface
to that of the in-memory implementation of the same name in :mod:`pyteomics.mass`.

Dependencies
------------

This module requires :py:mod:`lxml` and :py:mod:`sqlalchemy`.


Primary Interface
-----------------

    .. autoclass:: Unimod


Relational Entities
~~~~~~~~~~~~~~~~~~~
There are many tables that are described as object-relationally mapped (ORM) types in this module. The most important two are shown
here.

    .. class:: Modification

        A single modification record from Unimod, having an :attr:`id`, :attr:`full_name`, :attr:`code_name`,
        and :attr:`ex_code_name` as identifiers, and :attr:`monoisotopic_mass`, :attr:`average_mass`, and
        :attr:`composition` as mass-describing properties.

        Additional relationships may be loaded through :attr:`specificities` (see :class:`~.Specificity`), :attr:`alternative_names`,
        :attr:`fragments`, and :attr:`notes`.


    .. class:: Specificity

        Describes the relationship between a :class:`~.Modification` and an amino acid/position rule, along with the
        chemical process type that gives rise to that modification event.

Other ORM Types
***************

The following ORM types may be useful when composing a more detailed query. Additional types may be found in the source.

    .. class:: AminoAcid

    .. class:: Position

    .. class:: Classification

    .. class:: Fragment

    .. class:: AlternativeName

    .. class:: Crossreference