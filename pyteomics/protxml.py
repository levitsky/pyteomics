from . import xml, auxiliary as aux, _schema_defaults

class ProtXML(xml.XML):
    """Parser class for pepXML files."""
    file_format = 'protXML'
    _root_element = 'protein_summary'
    _default_schema = _schema_defaults._protxml_schema_defaults
    # _default_version = None
    _default_iter_tag = 'protein_group'
    _structures_to_flatten = {'annotation'}
    # attributes which contain unconverted values
    _convert_items = {'float':  {'pct_spectrum_ids'},
        'int': {'group_number', 'prot_length'},
        'bool': {'is_contributing_evidence', 'is_nondegenerate_evidence'}
        }.items()

    def _get_info_smart(self, element, **kwargs):
        """Extract the info in a smart way depending on the element type"""
        try:
            name = kwargs.pop('ename')
        except KeyError:
            name = xml._local_name(element)
        rec = kwargs.pop('recursive', None)
        if name == self._root_element:
            info = self._get_info(element, ename=name,
                    recursive=(rec if rec is not None else False),
                    **kwargs)
        else:
            info = self._get_info(element, ename=name,
                    recursive=(rec if rec is not None else True),
                    **kwargs)

        converters = {'float': float, 'int': int,
                'bool': lambda x: x.lower() in {'1', 'true', 'y'}}
        for k, v in dict(info).items():
            for t, s in self._convert_items:
                if k in s:
                    del info[k]
                    info[k] = converters[t](v)
        p = info.get('parameter')
        if isinstance(p, list) and len(p) == 1 and isinstance(p[0], dict):
            info.update(info.pop('parameter')[0])

        if 'modification_info' in info:
            # this is a list with one element
            info.update(info.pop('modification_info')[0])
        return info

def read(source, read_schema=False, iterative=True, **kwargs):
    """Parse `source` and iterate through protein groups.

    Parameters
    ----------
    source : str or file
        A path to a target protXML file or the file object itself.

    read_schema : bool, optional
        If :py:const:`True`, attempt to extract information from the XML schema
        mentioned in the protXML header. Otherwise, use default parameters.
        Not recommended without Internet connection or
        if you don't like to get the related warnings.

    iterative : bool, optional
        Defines whether iterative parsing should be used. It helps reduce
        memory usage at almost the same parsing speed. Default is
        :py:const:`True`.

    Returns
    -------
    out : ProtXML
       An iterator over dicts with protein group properties.
    """

    return ProtXML(source, read_schema=read_schema, iterative=iterative)

chain = aux._make_chain(read, 'read')

def DataFrame(*args, **kwargs):
    """Read protXML output files into a :py:class:`pandas.DataFrame`.

    Requires :py:mod:`pandas`.

    Parameters
    ----------
    *args, **kwargs : passed to :py:func:`chain`

    sep : str or None, optional
        Some values related to protein groups are variable-length lists.
        If `sep` is a :py:class:`str`, they will be packed into single string using
        this delimiter. If `sep` is :py:const:`None`, they are kept as lists. Default is
        :py:const:`None`.

    pd_kwargs : dict, optional
        Keyword arguments passed to the :py:class:`pandas.DataFrame` constructor.

    Returns
    -------
    out : pandas.DataFrame
    """
    import pandas as pd
    kwargs = kwargs.copy()
    sep = kwargs.pop('sep', None)
    pd_kwargs = kwargs.pop('pd_kwargs', {})
    def gen_items():
        with chain(*args, **kwargs) as f:
            for item in f:
                info = {}
                for k, v in item.items():
                    if isinstance(v, (str, int, float)):
                        info[k] = v
                if 'protein' in item:
                    for prot in item['protein']:
                        out = dict(info)
                        out.update(prot)
                        if 'unique_stripped_peptides' in out:
                            if sep is None:
                                out['unique_stripped_peptides'] = out['unique_stripped_peptides'].split('+')
                            else:
                                out['unique_stripped_peptides'] = sep.join(out['unique_stripped_peptides'].split('+'))
                        if 'indistinguishable_protein' in out:
                            if sep is None:
                                out['indistinguishable_protein'] = [p['protein_name'] for p in out['indistinguishable_protein']]
                            else:
                                out['indistinguishable_protein'] = sep.join(p['protein_name'] for p in out['indistinguishable_protein'])
                        yield out
    return pd.DataFrame(gen_items(), **pd_kwargs)
