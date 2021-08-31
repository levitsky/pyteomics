Unit Handling
=============

When parsing parsing a data file with unit information associated with a scalar value,
:mod:`pyteomics` uses annotated data types to carry units around, with the specific type
determined by parsing the value.

For instance, given the XML string:

.. code:: xml

    <cvParam cvRef="MS" accession="MS:1000927" name="ion injection time" value="68.227485656738"
             unitCvRef="UO" unitAccession="UO:0000028" unitName="millisecond"/>

This will be parsed into a :class:`~.unitfloat`:

.. code:: python

    >>> value = unitfloat(68.227485656738, "millisecond")
    # Get the unit name, perhaps to do a conversion to another unit
    >>> value.unit_info
    "millisecond"
    # Can be coerced into a plain float without issue
    >>> float(value)
    68.227485656738
    # Can be used identically to a normal float
    >>> value > 50.0 and value < 90.0
    True

To normalize the time unit, we can write a function like this:

.. code:: python

    from pyteomics.auxiliary import unitfloat

    def in_minutes(x):
        '''Convert a time quantity to minutes

        Parameters
        ----------
        x: unitfloat
            A float representing a quantity of time annotated with a time unit

        Returns
        -------
        unitfloat:
            The time after conversion to minutes
        '''
        try:
            unit = x.unit_info
        except AttributeError:
            return x
        if unit == 'minute':
            return x
        elif unit == 'second':
            y = unitfloat(x / 60., 'minute')
            return y
        elif unit == 'hour':
            y = unitfloat(x * 60, 'minute')
            return y
        else:
            warnings.warn("Time unit %r not recognized" % unit)
        return x

.. code:: python

    >>> seconds = unitfloat(93.5, "seconds")
    >>> minutes = in_minutes(seconds)
    >>> minutes
    1.55833

.. autoclass:: pyteomics.auxiliary.structures.unitfloat
    :no-inherited-members:

.. autoclass:: pyteomics.auxiliary.structures.unitint
    :no-inherited-members:

.. autoclass:: pyteomics.auxiliary.structures.unitstr
    :no-inherited-members: