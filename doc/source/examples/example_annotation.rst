Example 4: Spectrum Annotation
==============================

In this example we will retrieve a spectrum and visualize it, a lot like in `example 2 <example_msms.html>`_.
However, we will not read the spectrum from a file, but retrieve it directly from an online repository.
We will also use Pyteomics to annotate fragment ions. Let's get to it!

We are going to need `spectrum_utils <https://github.com/bittremieux/spectrum_utils>`_, a tool for spectrum
processing and visualization. Pyteomics integrates with it and allows to create annotated spectrum plots easily.


.. literalinclude:: ../_static/example_annotation.py
    :language: python

