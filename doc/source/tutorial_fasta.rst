Example 1: in silico cleavage, property calculation and visualization
=====================================================================

In this tutorial, we will introduce the Pyteomics tools to predict the basic
physicochemical characteristics of peptides, such as mass, charge and 
chromatographic retention time. We will download a FASTA database with baker's
yeast proteins, digest it with trypsin and study the distributions of
various quantitative qualities that may be measured in a typical proteomic
experiment. 

The tutorial is organized as a script interrupted by comments. It is 
assumed that the reader already has experience with numpy and matplotlib 
libraries. The source code for the tutorial can be found 
:download:`here <tutorial_fasta.py>`.

Before we begin, we need to import all the modules that we may require. Besides 
pyteomics itself, we need the tools that allow to access the hard drive (os), 
download files from the Internet (urllib), process and visualize
arrays of data (numpy, matplotlib).

.. literalinclude:: tutorial_fasta.py
   :language: python
   :lines: 1-7

We also need to download a real FASTA database. For our purposes, the Uniprot 
database with Saccharomyces cerevisiae proteins will work fine.

.. literalinclude:: tutorial_fasta.py
   :language: python
   :lines: 9-14

The :py:func:`pyteomics.fasta.read_fasta` function allows to iterate over the protein 
sequences in a fasta file in a standard python cycle. In order to obtain 
the peptide sequences, we cleave each protein using the 
:py:func:`pyteomics.parser.cleave` function and combine results into a set object 
that automatically discards multiple occurrences of the same sequence.

.. literalinclude:: tutorial_fasta.py
   :language: python
   :lines: 16-21

Later in the tutorial, we will calculate different peptide properties. In order
to store them, we create a list of dicts, where each dict stores the properties
of a single peptide, including its sequence.

.. literalinclude:: tutorial_fasta.py
   :language: python
   :lines: 23

It is also more efficient to pre-parse the sequences into individual amino acids
and supply the parsed structures into the functions that calculate m/z, charge, 
etc. During parsing, we explicitly save the terminal groups of peptides so
that they are taken into the account when calculating m/z and charge of a peptide.

.. literalinclude:: tutorial_fasta.py
   :language: python
   :lines: 25-31

For our purposes, we will limit ourselves to reasonably short peptides with
the length less than 100 residues.

.. literalinclude:: tutorial_fasta.py
   :language: python
   :lines: 33

We use :py:func:`pyteomics.electrochem.charge` to calculate the charge at pH=2.0.
The neural mass and m/z of an ion is found via 
:py:func:`pyteomics.mass.calculate_mass`.

.. literalinclude:: tutorial_fasta.py
   :language: python
   :lines: 35-42

Next, we calculate the retention time in the reversed- and normal-phase 
chromatography using :py:func:`pyteomics.achrom.calculate_RT` for two different
sets of retention coefficients. 
The phase is specified by supplying corresponding sets of retention
coefficients, :py:data:`pyteomics.achrom.RCs_zubarev` and
:py:data:`pyteomics.achrom.RCs_yoshida_lc` for the reversed and normal phases,
correspondingly.

.. literalinclude:: tutorial_fasta.py
   :language: python
   :lines: 44-52

Now, as we have all the numbers we can estimate the complexity of a sample
by plotting the distributions of parameters measurable in a typical proteomic
experiment. First, we show the distribution of m/z using the standard histogram 
plotting function from matplotlib.

.. literalinclude:: tutorial_fasta.py
   :language: python
   :lines: 54-60

The same set of commands allows us to plot the distribution of charge states 
in the sample:

.. literalinclude:: tutorial_fasta.py
   :language: python
   :lines: 62-68

Next we want to visualize the statistical correlation
between m/z and retention time in reversed-phase chromatography.

The standard approach would be to use a scatter plot. 
However, with a sample of our size that would be uninformative. Instead,
we will plot a 2d-histogram. There is no standard matplotlib command for that
and we have to use a combination of numpy and matplotlib. The function 
`numpy.histogram2d()` bins a set of (x,y) points on a plane and returns 
the matrix of numbers in each individual bin and the borders of the bins.
We also use a trick of replacing zeros in this matrix with the not-a-number 
value so that on the final figure empty bins are highlighted with white color 
instead of the darkest blue. We suggest removing the fourth line in this code 
snippet to see how that affects the final plot. At the last line, we also
apply the linear regression to obtain the coefficient of correlation between
m/z and retention time.

.. literalinclude:: tutorial_fasta.py
   :language: python
   :lines: 70-74

The obtained heatmap is plotted with `matplotlib.pyplot.imshow()` function
that visualizes matrices. 

.. literalinclude:: tutorial_fasta.py
   :language: python
   :lines: 76-85

The same code can also be applied to compare the retention times obtained on 
different chromatographic phases. 
As you can see upon execution of the code, the retention times obtained on 
different chromatographic phases seem to be uncorrelated.

.. literalinclude:: tutorial_fasta.py
   :language: python
   :lines: 87-105

Finally, let us check whether the retention times remain uncorrelated when
we narrow down the sample of peptides. We select the peptides with m/z lying in
a 700-701 Th window and plot two chromatographic retention times. This time
the sample allows us to use a scatter plot.

.. literalinclude:: tutorial_fasta.py
   :language: python
   :lines: 107-118

As you can see, the retention times of peptides lying in a narrow mass window
turn out to be substantially correlated.

At this point we stop. The next example will cover the modules allowing access
to experimental proteomic datasets stored in XML-based formats.
