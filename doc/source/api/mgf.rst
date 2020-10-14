.. automodule:: pyteomics.mgf

    .. autofunction:: chain
    .. py:function :: chain.from_iterable(files, **kwargs)

        Chain :py:func:`read` for several files.
        Keyword arguments are passed to the :py:func:`read` function.

        :param files: Iterable of file names or file objects.
        :type files: iterable

