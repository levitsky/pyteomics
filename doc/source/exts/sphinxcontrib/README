=====================
Autorun 
=====================

Autorun is an extension for Sphinx_ that can execute the code from a
runblock directive and attach the output of the execution to the document. 

For example::

    .. runblock:: pycon
        
        >>> for i in range(5):
        ...    print i

Produces::

    >>> for i in range(5):
    ...    print i
    1
    2
    3
    4
    5


Another example::

    .. runblock:: console

        $ date

Produces::

    $ date 
    Thu  4 Mar 2010 22:56:49 EST

Currently autorun supports ``pycon`` and ``console`` languages. It's also
possible to configure autorun (from `conf.py`) to run other languages.


Installation
-----------------

Installing from sources::

    $ hg clone http://bitbucket.org/birkenfeld/sphinx-contrib/
    $ cd sphinx-contrib/autorun
    $ python setup.py install

To enable autorun add 'sphinxcontrib.autorun' to the ``extension`` list in
`conf.py`::

    extensions.append('sphinxcontrib.autorun')

The documentation is in the doc/ folder.
