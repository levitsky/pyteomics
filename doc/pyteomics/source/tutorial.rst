========
Tutorial
========

Before we begin
***************

The following help is written both for libBioLCCC and pyBioLCCC. The only
difference between two these packages lies in the syntax of commands. That is
why we supply code snippets both for C++ and Python. Here is an example:

.. list-table:: Example of a code snippet
   :widths: 40 40
   :header-rows: 1

   * - C++
     - Python
   * - 

       .. literalinclude:: examples/snippet_example.cpp
          :language: cpp

     - 

       .. literalinclude:: examples/snippet_example.py
          :language: python


The Python examples are specific to Python 2.x, since our project doesn't
support Python 3.x.

Basic conceptions
*****************

There are a few simple conceptions which are widely used in libBioLCCC. Most of
them are represented by a corresponding class. Here they are:

**Polymer model** - the set of assumptions used to describe a polymer molecule.
This version of the BioLCCC model contains two models:

    **ROD** - this model represents a peptide as an absolutely rigid rod.
    Amino acids are modelled as regularly spaced beads
    threaded on this rod. This model describes peptides better comparing to
    long protein molecules.

    The equations for the ROD model are to be published in the upcoming
    paper.
    
    **CHAIN** - in this model a protein molecule is described as
    a free-joint chain of rods. The conformations of this molecule in a pore 
    can be modelled as a random walk in the field of adsorbing walls.
    This assumption should work better for long protein molecules.

    The CHAIN model was described in ''Liquid Chromatography at Critical 
    Conditions: Comprehensive Approach to Sequence-Dependent Retention Time 
    Prediction'', Alexander V. Gorshkov et al, Analytical Chemistry, 2006, 78
    (22), 7770-7777. `Link <http://dx.doi.org/10.1021/ac060913x>`_.

**Chemical group** - in libBioLCCC that is an amino acid residue OR a peptide
terminal group in a peptide chain. Examples are a histidine residue, 
phosphoserine residue and N-Terminal hydrogen that closes a peptide chain. The
properties of a chemical group are stored in the ChemicalGroup class. 

**Chemical basis** - a set of all physicochemical constants involved into the
BioLCCC equations. This set contains:

    - the list of all chemical groups, i.e. amino acids and terminal groups.
      Any peptide can be represented as a series of these, that is why it is
      a *basis* similar to the mathematical basis;
    - which terminal groups are set by default (cannon be changed);
    - the chemical properties of solvents: densities, molar mass and
      adsorption energies (adsorption energy of the first solvent always
      equals zero);
    - the model of a polymer molecule being used in calculations and
      approximations used in the equations;
    - peptide geometry: the length of amino acid and the Kuhn length;
    - the range of an interaction between an amino acid and the surface of 
      the solid phase (a.k.a. the width of the adsorbing layer).
       
The properties of a chemical basis are stored in the 
`ChemicalBasis class <./API/classBioLCCC_1_1ChemicalBasis.html>`_.

A chemical basis is specific to a type of retention chemistry, solvents
and ion paring agent being used in the experiment. In addition, it must be used
only with the same polymer model as the one used in the calibration of the
chemical basis.

**Predefined chemical basis** - a chemical basis, calculated (or, more
precisely, calibrated) for the specific retention chemistry and model of a
polymer molecule. The current version of libBioLCCC contains two predefined
chemical bases:

    **rpAcnFaRod** - a ChemicalBasis calibrated for the reversed phase,
    ACN as a second solvent, 0.1% FA in both solvents and the ROD polymer model.
    The data was obtained in the joint research of Harvard University and 
    Institute for Energy Problems for Chemical Physics, Russian Academy of
    Science.

    **rpAcnTfaChain** - a chemical basis calibrated for the reversed phase,
    ACN as a second solvent, 0.1% TFA in both solvents and the CHAIN model. 
    The initial data were taken from Guo et al, Journal of 
    Chromatography, 359 (1986) 449-517.

**Chromatographic conditions** - a description of a chromatographic equipment 
and its settings. Contains:

    - the geometry of the column.
    - the properties of the adsorbent: average size of the pores, porosity
      (i.e. percentage of volume not filled with the solid phase),
      (volume of pores)/(total volume of column) ratio, relative adsorption
      strength.
    - elution parameters: the shape of the gradient, the composition of
      components, flow rate, delay time.
    - the step of integration over volume.
    - temperature of a column (EXPERIMENTAL).

The default values were set rather arbitrarily.

Peptide sequence notation
*************************

In libBioLCCC we use the extended peptide notation. It is based on the
`one-letter IUPAC notation <http://www.chem.qmul.ac.uk/iupac/AminoAcid/>`_, 
but borrows only letters for the standard 20 aminoacid (i.e. no B, Z, X). 
We extended it in the following way:

- Modified amino acids are denoted as **xyzX**, i.e. their labels start with an 
  arbitrary number of lower-case letters and terminate with a single
  upper-case letter. The upper-case letter shows the base amino acid, while the
  lower-case letters describe the type of modification. The examples are:

    - **oxM** for oxidated methionine
    - **pS** for phosphorylated serine
    - **pT** for phosphorylated threonine
    - **camC** for carboxyamidomethylated cysteine

- The non-standard peptide terminal groups are denoted as **XxXx-** and
  **-XxXx**
  for N-terminal and C-terminal groups correspondingly. The label could contain
  an arbitrary number of mixed lower-case and upper-case letters and numbers, 
  but it should not be
  a valid peptide sequence. If a terminal group is not specified, it is
  assumed to be the standard one (i.e. an N-terminal hydrogen atom or C-terminal
  acidic group). The examples:
  
    - **Ac-** for N-Terminal acetylation
    - **H-** for N-Terminal hydrogen
    - **-NH2** for C-Terminal amidation
    - **-OH** for C-Terminal carboxyl group

- If a sequence contains two dots, then only the substring between them is
  parsed. This notation is used in several MS/MS search engines to show the
  adjacent amino acid residues for a peptide cleaved out of a protein. The
  examples are:

    -  K.APGFGDNR.K
    -  K.VGEVIVTK.D

Calculating retention time
**************************

calculateRT is the first libBioLCCC function you may need.
It requires three arguments: a peptide sequence,
a chemical basis, and and a description of chromatographic conditions. Supplied 
with these data, it calculates the retention time of the peptide.

.. list-table:: Calculating the retention time of a peptide
   :widths: 40 40
   :header-rows: 1

   * - C++
     - Python
   * - 

       .. literalinclude:: examples/rt_calculation.cpp
          :language: cpp

     - 

       .. literalinclude:: examples/rt_calculation.py
          :language: python

Please, consult with the 
`libBioLCCC API documentation <./API/namespaceBioLCCC.html>`_
for the details of calculateRT function.

Specifying chromatographic conditions
*************************************

The next thing you may need to learn is how to specify the chromatographic
conditions. In order to do that, create a new instance of ChromoConditions and
replace the default parameters with your own.

.. list-table:: Specifying chromatographic conditions
   :widths: 40 40
   :header-rows: 1

   * - C++
     - Python
   * - 

       .. literalinclude:: examples/chromoconditions.cpp
          :language: cpp

     - 

       .. literalinclude:: examples/chromoconditions.py
          :language: python


pyBioLCCC adds another way to interact with ChromoConditions. You can use its
instances as Python dictionaries:

.. list-table:: Dict-like syntax of ChromoConditions
   :widths: 40
   :header-rows: 1

   * - Python
   * - 

       .. literalinclude:: examples/chromoconditions_dict.py
          :language: python


Besides being more convenient and compact, this syntax allows ChromoConditions 
to be pickled. 

If you want to see the full list of parameters stored in a ChromoConditions
instance, please, take a look at the 
`class description <./API/classBioLCCC_1_1ChromoConditions.html>`_ 
in the libBioLCCC API documentation.

Calculating mass
****************

libBioLCCC contains functions to calculate the monoisotopic and average masses
of a peptide. Besides the sequence of a peptide, you need to specify a
ChemicalBasis instance which contains the masses of amino acids. 

.. list-table:: Calculating mass of a peptide
   :widths: 40 40
   :header-rows: 1

   * - C++
     - Python
   * - 

       .. literalinclude:: examples/mass_calculation.cpp
          :language: cpp

     - 

       .. literalinclude:: examples/mass_calculation.py
          :language: python

Getting the list of predefined chemical groups
**********************************************

Before you begin to work with libBioLCCC/pyBioLCCC, it is useful to know which
amino acids and terminal groups are predefined in this version of library.
To get this information just iterate through the chemicalGroups() map of the
predefined chemical bases.

.. list-table:: Examining a predefined chemical basis
   :widths: 40 40
   :header-rows: 1

   * - C++
     - Python
   * - 

       .. literalinclude:: examples/chemicalbasis.cpp
          :language: cpp

     - 

       .. literalinclude:: examples/chemicalbasis.py
          :language: python

..
    .. list-table:: example of a code snippet
       :widths: 40 40
       :header-rows: 1

       * - C++
         - Python
       * - 

           .. literalinclude:: ../../../src/examples/
              :language: cpp

         - 

           .. literalinclude:: ../../../src/examples/
              :language: python
