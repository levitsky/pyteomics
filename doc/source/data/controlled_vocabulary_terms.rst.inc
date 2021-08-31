Controlled Vocabulary Terms In Structured Data
==============================================
When parsing files with controlled vocabulary terms, especially in positions where
those terms may be used as keys, :mod:`pyteomics` will use an instance of :class:`~.pyteomics.auxiliary.structures.cvstr`
to represent that value. A :class:`~.pyteomics.auxiliary.structures.cvstr` is a sub-class of :class:`str`
with two additional attributes storing the accession number of the term it names, and another
holding the accession number of the unit of its value, if any.

.. autoclass:: pyteomics.auxiliary.structures.cvstr
    :no-inherited-members:

Accessing the attributes of a :class:`dict` key to find out if it matches a query
is inconvenient. To handle looking up a value by accession, the :class:`~pyteomics.auxiliary.structures.CVQueryEngine`
type can help solve the problem by either looking up a single accession value, or
convert a nested :class:`dict` structure with :class:`~.cvstr` as keys into a :class:`dict`
with accession numbers as keys, mapping to the value their owners pointed to in the original
:class:`dict`, or their naming :class:`str` if the value is empty:

For example, if we had parsed an mzML file, and read out a spectrum:

.. code:: python

    >>> from pyteomics import mzml
    >>> scan = next(mzml.read("tests/test.mzml"))
    >>> scan
    {'index': 0,
    'id': 'controllerType=0 controllerNumber=1 scan=1',
    'defaultArrayLength': 19914,
    'scanList': {'count': 1,
    'scan': [{'instrumentConfigurationRef': 'IC1',
        'scanWindowList': {'count': 1,
        'scanWindow': [{'scan window lower limit': 200.0 m/z,
        'scan window upper limit': 2000.0 m/z}]},
        'scan start time': 0.004935 minute,
        'filter string': 'FTMS + p ESI Full ms [200.00-2000.00]',
        'preset scan configuration': 1.0,
        '[Thermo Trailer Extra]Monoisotopic M/Z:': 810.4152221679688}],
    'no combination': ''},
    'ms level': 1,
    'MSn spectrum': '',
    'positive scan': '',
    'profile spectrum': '',
    'base peak m/z': 810.415283203125 m/z,
    'base peak intensity': 1471973.875 number of counts,
    'total ion current': 15245068.0,
    'lowest observed m/z': 200.00018816645022 m/z,
    'highest observed m/z': 2000.0099466203771 m/z,
    'count': 2,
    'm/z array': array([ 200.00018817,  200.00043034,  200.00067252, ..., 1999.96151259,
            1999.98572931, 2000.00994662]),
    'intensity array': array([0., 0., 0., ..., 0., 0., 0.], dtype=float32)}

Then :obj:`cvquery(scan)` would yield the following look up table:

.. code:: python

    >>> from pyteomics.auxiliary import cvquery
    >>> cvquery(scan)
    {'MS:1000501': 200.0 m/z,
     'MS:1000500': 2000.0 m/z,
     'MS:1000016': 0.004935 minute,
     'MS:1000512': 'FTMS + p ESI Full ms [200.00-2000.00]',
     'MS:1000616': 1.0,
     'MS:1000795': 'no combination',
     'MS:1000511': 1,
     'MS:1000580': 'MSn spectrum',
     'MS:1000130': 'positive scan',
     'MS:1000128': 'profile spectrum',
     'MS:1000504': 810.415283203125 m/z,
     'MS:1000505': 1471973.875 number of counts,
     'MS:1000285': 15245068.0,
     'MS:1000528': 200.00018816645022 m/z,
     'MS:1000527': 2000.0099466203771 m/z,
     'MS:1000514': array([ 200.00018817,  200.00043034,  200.00067252, ..., 1999.96151259,
             1999.98572931, 2000.00994662]),
     'MS:1000515': array([0., 0., 0., ..., 0., 0., 0.], dtype=float32)}

Alternatively, if we have a particular accession in mind, e.g. :obj:`MS:1000016` for :obj:`"scan start time"`, we could query for
that specifically:

.. code:: python

    >>> cvquery(scan, "MS:1000016")
    0.004935

Because :class:`~pyteomics.auxiliary.structures.CVQueryEngine` does not have any state of its own,
we use a pre-made instance, :obj:`pyteomics.auxiliary.structures.cvquery`.

.. autofunction:: pyteomics.auxiliary.structures.cvquery

.. autoclass:: pyteomics.auxiliary.structures.CVQueryEngine
