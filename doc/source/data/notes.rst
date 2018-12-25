General Notes
=============

- Each module mentioned below corresponds to a file format. In each module, the
  top-level function :py:func:`read` allows iteration over entries in a file.
  It works like the built-in :py:func:`open`, allowing direct iteration and
  supporting the ``with`` syntax, which we recommend using. So you can do:

  .. code-block :: python

       >>> from pyteomics import mgf
       >>> reader = mgf.read('tests/test.mgf')
       >>> for spectrum in reader:
       >>>    ...
       >>> reader.close()

  ... but it is recommended to do:

  .. code-block :: python

       >>> from pyteomics import mgf
       >>> with mgf.read('tests/test.mgf') as reader:
       >>>     for spectrum in reader:
       >>>        ...

- Additionally, most modules provide one or several classes which implement different
  parsing modes, e.g. :py:class:`pyteomics.mgf.MGF` and :py:class:`pyteomics.mgf.IndexedMGF`.
  Indexed parsers build an index of file entries and thus allow **random access** in addition
  to **iteration**. See `Indexed Parsers`_ for a detailed description and examples.

- Apart from :py:func:`read`, which reads just one file, all modules described
  here have functions for reading multiple files: :py:func:`chain` and
  :py:func:`chain.from_iterable`.
  ``chain('f1', 'f2')`` is equivalent to ``chain.from_iterable(['f1', 'f2'])``.
  :py:func:`chain` and :py:func:`chain.from_iterable` only support the
  ``with`` syntax. If you don't want to use the ``with`` syntax, you can just
  use the :py:mod:`itertools` functions :py:func:`chain` and
  :py:func:`chain.from_iterable`.

- Throughout this section we use
  :py:func:`pyteomics.auxiliary.print_tree` to display the structure of the
  data returned by various parsers. Replace this call with the actual processsing
  that you need to perform on your files.