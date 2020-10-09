.. automodule:: pyteomics.mzid

    .. autofunction:: version_info
    .. autofunction:: fdr
    .. autofunction:: qvalues
    .. autofunction:: chain
    .. py:function :: chain.from_iterable(files, **kwargs)

        Chain :py:func:`read` for several files.
        Keyword arguments are passed to the :py:func:`read` function.

        :param files: Iterable of file names or file objects.
        :type param: iterable

    .. autofunction:: filter
    .. py:function :: filter.chain(*files, **kwargs)

        Chain :py:func:`filter` for several files.
        Positional arguments should be file names or file objects.
        Keyword arguments are passed to the :py:func:`filter` function.

    .. py:function :: filter.chain.from_iterable(*files, **kwargs)

        Chain :py:func:`filter` for several files.
        Keyword arguments are passed to the :py:func:`filter` function.

        :param files: Iterable of file names or file objects.
        :type param: iterable

