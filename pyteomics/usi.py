"""
usi - Universal Spectrum Identifier (USI) parser and minimal PROXI client
=========================================================================

Summary
-------
`USI <http://www.psidev.info/usi>`_ is a standardized method of referencing a specific
spectrum in a dataset, possibly attached to an interpretation. This module includes a
:class:`USI` type which can represent these constructs, :meth:`~USI.parse` them and
reconstruct them.

One use-case for USI is to request spectrum information from a `PROXI <http://www.psidev.info/proxi>`_
service host. PROXI services are available from several of the major national proteomics data hosts,
including MassIVE, PeptideAtlas, PRIDE, and jPOST.


Data access
-----------

  :py:class:`USI` for representing Universal Spectrum Identifiers. Call :meth:`USI.parse` to parse a USI
  string.

  :py:func:`proxi` to request a USI from a remote service. Provides access to the PeptideAtlas, MassIVE,
  PRIDE and jPOST hosts.

"""
import json
from collections import namedtuple

try:
    from urllib2 import Request, urlopen
except ImportError:
    from urllib.request import Request, urlopen

try:
    import numpy as np

    def coerce_array(array_data):
        return np.array([float(v) for v in array_data])

except ImportError:

    def coerce_array(array_data):
        return [float(v) for v in array_data]


class USI(namedtuple("USI", ['protocol', 'dataset', 'datafile', 'scan_identifier_type', 'scan_identifier', 'interpretation'])):
    '''Represent a Universal Spectrum Identifier (USI).

    .. note::
        This implementation will capture the interpretation component but will not interpret it at this time.

    Attributes
    ----------
    protocol: str
        The protocol to use to access the data (usually mzspec)
    dataset: str
        The name or accession number for the dataset the spectrum residues in
    datafile: str
        The basename of the data file from :attr:`dataset` to retrieve the spectrum from
    scan_identifier_type: str
        The format of the scan identifier, one of (scan, index, nativeId, trace)
    scan_identifier: str
        A usually numerical but potentially comma separated value encoded as a string to uniquely
        identify the spectrum to be recovered from :attr:`datafile` in :attr:`dataset`.
    interpretation: str
        The trailing material of the USI, such as the ProForma peptide sequence and charge
    '''
    def __str__(self):
        return ':'.join(filter(lambda x: x is not None, self))

    @classmethod
    def parse(cls, usi):
        '''Parse a USI string into a :class:`USI` object.

        Parameters
        ----------
        usi: str
            The USI string to parse

        Returns
        -------
        USI
        '''
        return cls(*_usi_parser(str(usi)))


def _usi_parser(usi):
    tokens = usi.split(":", 5)
    protocol = tokens[0]
    dataset = tokens[1]
    datafile = tokens[2]
    scan_identifier_type = tokens[3]
    scan_identifier = tokens[4]
    try:
        interpretation = tokens[5]
    except IndexError:
        interpretation = None
    return (protocol, dataset, datafile, scan_identifier_type, scan_identifier, interpretation)



class _PROXIBackend(object):
    '''A base class for all PROXI backends to implement the gory details of HTTP requests
    and protocol parsing.

    If special processing needs to be done to interpret the spectrum returned from the service
    provider, override the :meth:`_coerce` method.

    If extra information needs to be provided to the service provider for them to fulfill the
    request not passed through the URL, override the :meth:`_request` method.

    Attributes
    ----------
    name: str
        The name of the backend service
    url_template: str
        The URL with {} fields to populate with the USI and any other relevant options, like protocol version
        or the like.
    options: dict
        Additional options to be used when preparing the request URL
    '''
    def __init__(self, name, url_template, **kwargs):
        kwargs.setdefault('version', '0.1')
        self.name = name
        self.url_template = url_template
        self.options = kwargs

    def __repr__(self):
        return "{self.__class__.__name__}({self.options})".format(self=self)

    def _request(self, usi):
        url = self.url_template.format(usi=usi, **self.options)
        req = Request(url)
        response = urlopen(req)
        if response.getcode() != 200:
            raise ValueError("PROXI Service Response Code %r" % (response.getcode()))
        data = response.read().decode("utf-8")
        data = json.loads(data)
        return data

    def get(self, usi):
        '''Retrieve a ``USI`` from the host PROXI service over the network.

        Parameters
        ----------
        usi : str or :class:`USI`
            The universal spectrum identifier to retrieve.

        Returns
        -------
        dict:
            The spectrum as represented by the requested PROXI host.
        '''
        data = self._request(usi)
        result = self._coerce(data)
        return result

    def _coerce(self, data):
        '''Override and extend this method to change how the spectrum information is refined.

        This implementation just deals with properly formatting the peak arrays and doing minor
        cosmetic name normalization.

        Parameters
        ----------
        data: dict
            The raw mzSpecML representation parsed from JSON

        Returns
        -------
        dict:
            The coerced spectrum data of appropriate types
        '''
        if isinstance(data, list):
            data_collection = data
            data = data_collection[0]
        result = {}
        result['attributes'] = data.pop('attributes', {})
        result['m/z array'] = coerce_array(data.pop('mzs', []))
        result['intensity array'] = coerce_array(data.pop('intensities', []))
        for key, value in data.items():
            if key in result:
                raise ValueError(
                    "Attempting to set explicit value for {key!r}".format(key=key))
            result[key] = value
        return result

    def __call__(self, usi):
        return self.get(usi)


class PeptideAtlasBackend(_PROXIBackend):
    _url_template = "http://www.peptideatlas.org/api/proxi/v{version}/spectra?resultType=full&usi={usi!s}"

    def __init__(self, **kwargs):

        super(PeptideAtlasBackend, self).__init__(
            'PeptideAtlas', self._url_template, **kwargs)


class MassIVEBackend(_PROXIBackend):

    _url_template = "http://massive.ucsd.edu/ProteoSAFe/proxi/v{version}/spectra?resultType=full&usi={usi}"

    def __init__(self, **kwargs):
        super(MassIVEBackend, self).__init__(
            'MassIVE', self._url_template, **kwargs)


class PRIDEBackend(_PROXIBackend):
    _url_template = "http://wwwdev.ebi.ac.uk/pride/proxi/archive/v{version}/spectra?resultType=full&usi={usi}"

    def __init__(self, **kwargs):
        super(PRIDEBackend, self).__init__(
            'PRIDE', self._url_template, **kwargs)


class JPOSTBackend(_PROXIBackend):
    _url_template = 'https://repository.jpostdb.org/spectrum/?USI={usi}'

    def __init__(self, **kwargs):
        super(JPOSTBackend, self).__init__('jPOST', self._url_template, **kwargs)
        kwargs.pop("version", None)



_proxies = {
    "peptide_atlas": PeptideAtlasBackend,
    "massive": MassIVEBackend,
    "pride": PRIDEBackend,
    "jpost": JPOSTBackend,
}

def proxi(usi, backend='peptide_atlas', **kwargs):
    '''Retrieve a ``USI`` from a `PROXI <http://www.psidev.info/proxi>`.

    Parameters
    ----------
    usi : str or :class:`USI`
        The universal spectrum identifier to request.
    backend : str or :class:`Callable`
        Either the name of a PROXI host (peptide_atlas, massive, pride, or jpost), or a
        callable object (which :class:`_PROXIBackend` instances are) which will be used
        to resolve the USI.
    **kwargs:
        extra arguments passed when constructing the backend by name.

    Returns
    -------
    dict :
        The spectrum as represented by the requested PROXI host.
    '''
    if isinstance(backend, str):
        backend = _proxies[backend](**kwargs)
    elif issubclass(backend, _PROXIBackend):
        backend = backend(**kwargs)
    elif callable(backend):
        backend = backend
    else:
        raise TypeError("Unrecognized backend type")
    return backend(usi)
