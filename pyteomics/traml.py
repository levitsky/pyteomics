"""
traml - targeted MS transition data in TraML format
===================================================

Summary
-------

TraML is a standard rich XML-format for targeted mass spectrometry method definitions.
Please refer to `psidev.info <http://www.psidev.info/traml>`_
for the detailed specification of the format and structure of TraML files.

This module provides a minimalistic way to extract information from TraML
files. You can use the object-oriented interface (:class:`TraML` instances) to
access target definitions and transitions. :class:`TraML` objects also support
indexing with entity IDs directly.

Data access
-----------

  :py:class:`TraML` - a class representing a single TraML file.
  Other data access functions use this class internally.

  :py:func:`read` - iterate through transitions in TraML format.

  :py:func:`chain` - read multiple TraML files at once.

  :py:func:`chain.from_iterable` - read multiple files at once, using an
  iterable of files.

Controlled Vocabularies
~~~~~~~~~~~~~~~~~~~~~~~
TraML relies on controlled vocabularies to describe its contents extensibly. See
`Controlled Vocabulary Terms <../data.html#controlled-vocabulary-terms-in-structured-data>`_
for more details on how they are used.

Handling Time Units and Other Qualified Quantities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TraML contains information which may be described as using a variety of different time units.
See `Unit Handling <../data.html#unit-handling>`_ for more information.

Deprecated functions
--------------------

  :py:func:`version_info` - get version information about the TraML file.
  You can just read the corresponding attribute of the :py:class:`TraML` object.

  :py:func:`iterfind` - iterate over elements in an TraML file.
  You can just call the corresponding method of the :py:class:`TraML` object.

Dependencies
------------

This module requires :py:mod:`lxml`

-------------------------------------------------------------------------------
"""

#   Copyright 2018 Joshua Klein, Lev Levitsky
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.


import warnings
from . import xml, _schema_defaults, auxiliary as aux


class TraML(xml.MultiProcessingXML, xml.IndexSavingXML):
    """Parser class for TraML files."""
    file_format = 'TraML'
    _root_element = 'TraML'
    _default_schema = _schema_defaults._traml_schema_defaults
    _default_version = '1.0.0'

    _default_iter_tag = 'Transition'
    _indexed_tags = {
        'Transition',
        'Peptide',
        'Compound',
        'Target',
        'Protein',
        'Compound',
    }

    _element_handlers = xml.XML._element_handlers.copy()
    _element_handlers.update({
        'Modification': xml.XML._promote_empty_parameter_to_name,
        'Interpretation': xml.XML._promote_empty_parameter_to_name,
        'Software': xml.XML._promote_empty_parameter_to_name,
    })

    def __init__(self, *args, **kwargs):
        kwargs.setdefault('retrieve_refs', True)
        super(TraML, self).__init__(*args, **kwargs)

    def _get_info_smart(self, element, **kw):
        kwargs = dict(kw)
        rec = kwargs.pop('recursive', None)
        info = self._get_info(
            element,
            recursive=(rec if rec is not None else True),
            **kwargs)
        return info

    def _retrieve_refs(self, info, **kwargs):
        """Retrieves and embeds the data for each attribute in `info` that
        ends in `Ref`. Removes the id attribute from `info`"""
        for k, v in dict(info).items():
            if k[-3:] in {'Ref', 'ref'}:
                if isinstance(v, str):
                    key = v
                elif isinstance(v, dict):
                    key = v['ref']
                else:
                    if k != 'ref':
                        info[k[:-3]] = info.pop(k)
                    continue
                try:
                    by_id = self.get_by_id(key, retrieve_refs=True)
                except KeyError:
                    warnings.warn('Ignoring unresolved reference: ' + key)
                else:
                    if k == 'ref':
                        info.update(by_id)
                    else:
                        # by_id.pop('id', None)
                        info[k[:-3]] = by_id
                        del info[k]



def read(source, retrieve_refs=True, read_schema=False, iterative=True, use_index=False, huge_tree=False):
    """Parse `source` and iterate through transitions.

    Parameters
    ----------
    source : str or file
        A path to a target TraML file or the file object itself.

    retrieve_refs : bool, optional
        If :py:const:`True`, additional information from references will be
        automatically added to the results. The file processing time will
        increase. Default is :py:const:`True`.

    read_schema : bool, optional
        If :py:const:`True`, attempt to extract information from the XML schema
        mentioned in the TraML header. Otherwise, use default parameters.
        Not recommended without Internet connection or
        if you don't like to get the related warnings.

    iterative : bool, optional
        Defines whether iterative parsing should be used. It helps reduce
        memory usage at almost the same parsing speed. Default is
        :py:const:`True`.

    use_index : bool, optional
        Defines whether an index of byte offsets needs to be created for
        spectrum elements. Default is :py:const:`False`.

    huge_tree : bool, optional
        This option is passed to the `lxml` parser and defines whether
        security checks for XML tree depth and node size should be disabled.
        Default is :py:const:`False`.
        Enable this option for trusted files to avoid XMLSyntaxError exceptions
        (e.g. `XMLSyntaxError: xmlSAX2Characters: huge text node`).

    Returns
    -------
    out : TraML
       A :py:class:`TraML` object, suitable for iteration and possibly random access.
    """

    return TraML(source, retrieve_refs=retrieve_refs, read_schema=read_schema, iterative=iterative,
                 use_index=use_index, huge_tree=huge_tree)


def iterfind(source, path, **kwargs):
    """Parse `source` and yield info on elements with specified local
    name or by specified "XPath".

    .. note:: This function is provided for backward compatibility only.
        If you do multiple :py:func:`iterfind` calls on one file, you should
        create an :py:class:`TraML` object and use its
        :py:meth:`!iterfind` method.

    Parameters
    ----------
    source : str or file
        File name or file-like object.

    path : str
        Element name or XPath-like expression. Only local names separated
        with slashes are accepted. An asterisk (`*`) means any element.
        You can specify a single condition in the end, such as:
        ``"/path/to/element[some_value>1.5]"``
        Note: you can do much more powerful filtering using plain Python.
        The path can be absolute or "free". Please don't specify
        namespaces.

    recursive : bool, optional
        If :py:const:`False`, subelements will not be processed when
        extracting info from elements. Default is :py:const:`True`.

    iterative : bool, optional
        Specifies whether iterative XML parsing should be used. Iterative
        parsing significantly reduces memory usage and may be just a little
        slower. When `retrieve_refs` is :py:const:`True`, however, it is
        highly recommended to disable iterative parsing if possible.
        Default value is :py:const:`True`.

    read_schema : bool, optional
        If :py:const:`True`, attempt to extract information from the XML schema
        mentioned in the mzIdentML header. Otherwise, use default
        parameters. Not recommended without Internet connection or
        if you don't like to get the related warnings.

    Returns
    -------
    out : iterator
    """
    return TraML(source, **kwargs).iterfind(path, **kwargs)


version_info = xml._make_version_info(TraML)

chain = aux.ChainBase._make_chain(TraML)
