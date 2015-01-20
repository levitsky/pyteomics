.. automodule:: pyteomics.mzid
    :members:

    .. autofunction:: version_info
    .. autofunction:: iterfind
    .. autofunction:: fdr
    .. autofunction:: local_fdr
    .. autofunction:: chain
    .. py:function :: chain.from_iterable(files, **kwargs)

        Chain :py:func:`read` for several files.
        Keyword arguments are passed to the :py:func:`read` function.

        Parameters
        ----------
        files : iterable
            Iterable of file names or file objects.

    .. autofunction:: filter
    .. py:function :: filter.chain(*files, **kwargs)

        Chain :py:func:`filter` for several files.
        Positional arguments should be file names or file objects.
        Keyword arguments are passed to the :py:func:`filter` function.

    .. py:function :: filter.chain.from_iterable(*files, **kwargs)

        Chain :py:func:`filter` for several files.
        Keyword arguments are passed to the :py:func:`filter` function.

        Parameters
        ----------
        files : iterable
            Iterable of file names or file objects.

