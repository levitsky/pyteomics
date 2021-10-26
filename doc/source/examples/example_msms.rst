Example 2: Fragmentation
========================

.. note::
    Pyteomics has come a long way since this example was written.
    Check `example 4 <example_annotation.html>`_ for new Pyteomics tools you should know about.

In this example, we are going to retrieve MS/MS data from an MGF file and
compare it to identification info we read from a pepXML file. We are going to
compare the MS/MS spectrum in the file with the theoretical spectrum of a
peptide assigned to this spectrum by the search engine.

The script source can be downloaded :download:`here <../_static/example_msms.py>`. We will
also need the :download:`example MGF file <../_static/example.mgf>` and the
:download:`example pepXML file <../_static/example.pep.xml>`, but the script will download
them for you.

The MGF file has a single MS/MS spectrum in it. This spectrum is taken from the
`SwedCAD database of annotated MS/MS spectra <http://www.ncbi.nlm.nih.gov/pubmed/17711326>`_.
The pepXML file was obtained by running X!Tandem against the MGF file and
converting the results to pepXML with
`the Tandem2XML tool from TPP <http://tools.proteomecenter.org/wiki/index.php?title=Software:Tandem2XML>`_.

Let's start with importing the modules.

.. literalinclude:: ../_static/example_msms.py
    :language: python
    :lines: 1-4

Then we'll download the files, if needed:

.. literalinclude:: ../_static/example_msms.py
    :language: python
    :lines: 7-15

Now it's time to define the function that will give us *m/z* of theoretical
fragments for a given sequence. We will use
:py:func:`pyteomics.mass.fast_mass` to calculate the values.
All we need to do is split the sequence at every bond and iterate
over possible charges and ion types:

.. literalinclude:: ../_static/example_msms.py
    :language: python
    :lines: 17-30

So, the outer loop is over "fragmentation sites", the next one is over ion
types, then over charges, and lastly over two parts of the sequence
(C- and N-terminal).

All right, now it's time to extract the info from the files.
We are going to use the `with` statement syntax, which is not required, but
recommended.

.. literalinclude:: ../_static/example_msms.py
    :language: python
    :lines: 32-34

Now prepare the figure...

.. literalinclude:: ../_static/example_msms.py
    :language: python
    :lines: 36-40

... plot the real spectrum:

.. literalinclude:: ../_static/example_msms.py
    :language: python
    :lines: 41-42

... calculate and plot the theoretical spectrum, and show everything:

.. literalinclude:: ../_static/example_msms.py
    :language: python
    :lines: 43-48

You will see something like this:

.. image:: ../_static/example_msms.png

That's it, as you can see, the most intensive peaks in the spectrum are indeed
matched by the theoretical spectrum.
