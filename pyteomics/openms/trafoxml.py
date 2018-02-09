"""
trafoxml - reader for trafoXML files
========================================

Summary
-------

**trafoXML** is a format specified in the
`OpenMS <http://open-ms.sourceforge.net/about/>`_ project.
It defines a transformation, which is a result of retention time alignment.

This module provides a minimalistic way to extract information from **trafoXML**
files. You can use the old functional interface (:py:func:`read`) or the new
object-oriented interface (:py:class:`TrafoXML`)
to iterate over entries in ``<Pair>`` elements.

Data access
-----------

  :py:class:`TrafoXML` - a class representing a single trafoXML file.
  Other data access functions use this class internally.

  :py:func:`read` - iterate through pairs in a trafoXML file. Data from a
  single trafo are converted to a human-readable dict.

  :py:func:`chain` - read multiple trafoXML files at once.

  :py:func:`chain.from_iterable` - read multiple files at once, using an
  iterable of files.

Dependencies
------------

This module requres :py:mod:`lxml`.

--------------------------------------------------------------------------------
"""

from .. import xml, auxiliary as aux, _schema_defaults

class TrafoXML(xml.XML):
    """Parser class for trafoXML files."""
    file_format = 'trafoXML'
    _root_element = 'TrafoXML'
    _default_schema = _schema_defaults._trafoxml_schema_defaults
    _default_version = '1.0'
    _default_iter_tag = 'Pair'
    _schema_location_param = 'noNamespaceSchemaLocation'

    def _get_info_smart(self, element, **kw):
        kw['recursive'] = kw.get('recursive', True)
        info = self._get_info(element, **kw)
        return info

def read(source, read_schema=True, iterative=True):
    """Parse `source` and iterate through pairs.

    Parameters
    ----------
    source : str or file
        A path to a target trafoXML file or the file object itself.

    read_schema : bool, optional
        If :py:const:`True`, attempt to extract information from the XML schema
        mentioned in the file header (default). Otherwise, use default
        parameters. Disable this to avoid waiting on slow network connections or
        if you don't like to get the related warnings.

    iterative : bool, optional
        Defines whether iterative parsing should be used. It helps reduce
        memory usage at almost the same parsing speed. Default is
        :py:const:`True`.

    Returns
    -------
    out : iterator
       An iterator over the dicts with feature properties.
    """

    return TrafoXML(source, read_schema=read_schema, iterative=iterative)

chain = aux._make_chain(read, 'read')