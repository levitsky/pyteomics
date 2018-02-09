"""
featurexml - reader for featureXML files
========================================

Summary
-------

**featureXML** is a format specified in the
`OpenMS <http://open-ms.sourceforge.net/about/>`_ project.
It defines a list of LC-MS features observed in an experiment.

This module provides a minimalistic way to extract information from **featureXML**
files. You can use the old functional interface (:py:func:`read`) or the new
object-oriented interface (:py:class:`FeatureXML`)
to iterate over entries in ``<feature>`` elements.
:py:class:`FeatureXML` also supports direct indexing with feature IDs.

Data access
-----------

  :py:class:`FeatureXML` - a class representing a single featureXML file.
  Other data access functions use this class internally.

  :py:func:`read` - iterate through features in a featureXML file. Data from a
  single feature are converted to a human-readable dict.

  :py:func:`chain` - read multiple featureXML files at once.

  :py:func:`chain.from_iterable` - read multiple files at once, using an
  iterable of files.

Dependencies
------------

This module requres :py:mod:`lxml`.

--------------------------------------------------------------------------------
"""

from .. import xml, auxiliary as aux, _schema_defaults

class FeatureXML(xml.IndexedXML):
    """Parser class for featureXML files."""
    file_format = 'featureXML'
    _root_element = 'featureMap'
    _default_schema = _schema_defaults._featurexml_schema_defaults
    _default_version = '1.6'
    _default_iter_tag = 'feature'
    _structures_to_flatten = {}
    _indexed_tags = {'feature'}
    _schema_location_param = 'noNamespaceSchemaLocation'

    def _get_info_smart(self, element, **kw):
        kw['recursive'] = kw.get('recursive', True)
        info = self._get_info(element, **kw)
        return info

def read(source, read_schema=True, iterative=True, use_index=False):
    """Parse `source` and iterate through features.

    Parameters
    ----------
    source : str or file
        A path to a target featureXML file or the file object itself.

    read_schema : bool, optional
        If :py:const:`True`, attempt to extract information from the XML schema
        mentioned in the file header (default). Otherwise, use default
        parameters. Disable this to avoid waiting on slow network connections or
        if you don't like to get the related warnings.

    iterative : bool, optional
        Defines whether iterative parsing should be used. It helps reduce
        memory usage at almost the same parsing speed. Default is
        :py:const:`True`.

    use_index : bool, optional
        Defines whether an index of byte offsets needs to be created for
        spectrum elements. Default is :py:const:`False`.

    Returns
    -------
    out : iterator
       An iterator over the dicts with feature properties.
    """

    return FeatureXML(source, read_schema=read_schema, iterative=iterative, use_index=use_index)

chain = aux._make_chain(read, 'read')