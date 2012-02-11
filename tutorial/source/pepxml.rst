.pep.XML
========

**.pep.XML** is a widely used XML-based format for peptide identifications.
It contains information about the MS data, as well as about the search engine
used and the assigned sequences. To access these data, use :py:mod:`pepxml.py`
module.
The main function is :py:func:`iter_psm` which generates a sequence of
*dictioinaries* which contain the peptide-spectrum match information. The
following code snippet reads a .pep.xml.gz file and calculates a histogram of
mass differences between experimental masses and those of assigned peptides.

::

    from pyteomics.pepxml import iter_psm
    import gzip
    import pylab

    mass_diffs = [psm[‘massdiff’] for psm in iter_psm(gzip.open(‘/path/to/file/output.pep.xml.gz’))]
    pylab.figure()
    pylab.hist(mass_diffs)
    pylab.show()

