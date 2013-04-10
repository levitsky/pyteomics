Example 2: Fragmentation
========================

In this example, we are going to retrieve MS/MS data from an MGF file and
compare it to identification info we read from a pepXML file. We are going to
compare the MS/MS spectrum in the file with the theoretical spectrum of a
peptide assigned to this spectrum by the search engine.

The script source can be downloaded :download:`here <example_msms.py>`. We will
also need the :download:`example MGF file <example.mgf>` and the
:download:`example pepXML file <example.pep.xml>`, but the script will download
them for you.

The MGF file has a single MS/MS spectrum in it. This spectrum is taken from the
`SwedCAD database of annotated MS/MS spectra <http://www.ncbi.nlm.nih.gov/pubmed/17711326>`_.
The pepXML file was obtained by running X!Tandem against the MGF file and
converting the results to pepXML with
`the Tandem2XML tool from TPP <http://tools.proteomecenter.org/wiki/index.php?title=Software:Tandem2XML>`_.

Let's start with importing the modules.

.. literalinclude:: example_msms.py
    :language: python
    :lines: 1-7

Then we'll download the files, if needed:

.. literalinclude:: example_msms.py
    :language: python
    :lines: 10-13

Now it's time to define the function that will give us *m/z* of theoretical
fragments for a given sequence. We will use
:py:func:`pyteomics.mass.fast_mass` to calculate the values.
All we need to do is split the sequence at every bond and iterate
over possible charges and ion types:

.. literalinclude:: example_msms.py
    :language: python
    :lines: 15-28

So, the outer loop is over "fragmentation sites", the next one is over ion
types, then over charges, and lastly over two parts of the sequence
(C- and N-terminal).

All right, now it's time to extract the info from the files.
We are going to use the `with` statement syntax, which is not required, but
recommended.

.. literalinclude:: example_msms.py
    :language: python
    :lines: 29-32

Now prepare the figure...

.. literalinclude:: example_msms.py
    :language: python
    :lines: 33-37

... plot the real spectrum:

.. literalinclude:: example_msms.py
    :language: python
    :lines: 38-39

... calculate and plot the theoretical spectrum, and show everything:

.. literalinclude:: example_msms.py
    :language: python
    :lines: 40-45

You will see something like :download:`this <example_msms.png>`.

That's it, as you can see, the most intensive peaks in the spectrum are indeed
matched by the theoretical spectrum.
