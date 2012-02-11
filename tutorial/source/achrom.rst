Retention time prediction
=========================

:py:mod:`achrom.py` is a module that implements the *additive model of
polypeptide chromatography* and thus allows prediction of chromatographic
retention times. 

.. seealso::

    Please refer to :py:mod:`achrom.py` documentation page for more
    information about the additive model.
   
One can find a set of *retention coefficients* that fit a given set of
chromatographic data and use these coefficients (along with the standard ones
from the literature) to predict peptide retention times.

Examples
--------

The following example reads the *data.txt* file where sequences and experimental
retention times are stored in the following manner::

    PEPTIDE 35.1
    PEPTID 30.2
    DEPTH 29.1
    ...

::

    data = open(‘data.txt’)
    sequences = []
    RTs = []
    for line in data:
        seq, time = line.split()
        sequences.append(seq.strip())
        RTs.append(float(time.strip()))
        
.. note::

     1. You might need to specify the full path to your datafile, unless it
        is in the working directory.
     2. Do not type sequences = RTs = [] ! This is a common Python pitfall.

Now, if you have enough data, you can use :py:data:`sequences` and :py:data:`RTs` to calibrate the additive model:

::

    >>> import pyteomics.achrom as ac
    >>> RCs = ac.get_RCs(sequences, RTs)

or::

    >>> RCs = ac.get_RCs(sequences, RTs, length_correction_factor=-0.19)

or::

    >>> RCs = ac.get_RCs_vary_lcf(sequences, RTs)

.. seealso::
    
    For more information about optional parameters, see 
    :py:mod:`achrom.py` documentation.

Having found the retention coefficients, we can use them to predict retention
times for other peptides::

    myPeptides = [‘PEPTIDE’*2, ‘PEPTIDE’*3]
    predicted_RTs = [ac.calculate_RT(peptide, RCs) for peptide in myPeptides]

.. note::

    1. :py:data:`RCs` were found in the previous snippet.
    2. You might want to read about **list comprehension in Python** to
       understand the used syntax better.
    3. There are other sets of retention coefficients coded in the module that
       we found in the literature. You can find the references in the
       documentation. 
