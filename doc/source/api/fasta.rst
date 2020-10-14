.. automodule:: pyteomics.fasta
   :exclude-members: Protein

   .. autofunction:: chain
   .. py:function :: chain.from_iterable(files, **kwargs)

        Chain :py:func:`read` for several files.
        Keyword arguments are passed to the :py:func:`read` function.

        :param files: Iterable of file names or file objects.
        :type param: iterable

   .. autofunction:: decoy_chain
   .. py:function :: decoy_chain.from_iterable(files, **kwargs)

        Chain :py:func:`decoy_db` for several files.
        Keyword arguments are passed to the :py:func:`decoy_db` function.

        :param files: Iterable of file names or file objects.
        :type param: iterable
