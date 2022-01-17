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

.. seealso::
   LeDuc, Richard D., Eric W. Deutsch, Pierre-Alain Binz, Ryan T. Fellers, Anthony J. Cesnik,
   Joshua A. Klein, Tim Van Den Bossche, et al.
   "Proteomics Standards Initiative's ProForma 2.0: Unifying the Encoding of Proteoforms and Peptidoforms."
   ArXiv:2109.11352 [q-Bio], September 23, 2021. http://arxiv.org/abs/2109.11352.



Data access
-----------

  :py:class:`USI` for representing Universal Spectrum Identifiers. Call :meth:`USI.parse` to parse a USI
  string.

  :py:func:`proxi` to request a USI from a remote service. Provides access to the PeptideAtlas, MassIVE,
  PRIDE and jPOST hosts.

"""
import json
import warnings
import threading
import multiprocessing

from collections import namedtuple, defaultdict

try:
    from multiprocessing.dummy import Pool as ThreadPool
except ImportError:
    ThreadPool = None

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

from .auxiliary import PyteomicsError


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


def cast_numeric(value):
    try:
        return int(value)
    except ValueError:
        pass
    try:
        return float(value)
    except ValueError:
        return value


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
        Additional options to be used when preparing the request URL.
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
        result['attributes'] = data.pop('attributes', [])
        for attrib in result['attributes']:
            if 'value' in attrib and isinstance(attrib['value'], str) and attrib['value'][0].isdigit():
                try:
                    attrib['value'] = cast_numeric(attrib['value'])
                except TypeError:
                    continue
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
    _url_template = 'https://repository.jpostdb.org/proxi/spectra?resultType=full&usi={usi}'

    def __init__(self, **kwargs):
        super(JPOSTBackend, self).__init__('jPOST', self._url_template, **kwargs)
        kwargs.pop("version", None)


class ProteomeExchangeBackend(_PROXIBackend):
    _url_template = 'http://proteomecentral.proteomexchange.org/api/proxi/v{version}/spectra?resultType=full&usi={usi!s}'

    def __init__(self, **kwargs):

        super(ProteomeExchangeBackend, self).__init__(
            'ProteomeExchange', self._url_template, **kwargs)


class PROXIAggregator(object):
    '''Aggregate across requests across multiple PROXI servers.

    Will attempt to coalesce responses from responding servers into a single spectrum
    representation.

    Attributes
    ----------
    backends : :class:`dict` mapping :class:`str` to :class:`_PROXIBackend`
        The backend servers to query. Defaults to the set of all available backends.
    n_threads : int
        The number of threads to run concurrently to while making requests. Defaults
        to the number of servers to query.
    timeout : float
        The number of seconds to wait for a response.
    ephemeral_pool : bool
        Whether or not to tear down the thread pool between requests.
    '''

    _coalesce_resolution_methods = ("first", )

    def __init__(self, backends=None, n_threads=None, timeout=15, merge=True, ephemeral_pool=True, **kwargs):
        if backends is None:
            backends = {k: v() for k, v in _proxies.items()}
        if n_threads is None:
            n_threads = len(backends)

        self.lock = threading.RLock()

        self.timeout = timeout
        self.backends = backends
        self.n_threads = n_threads
        self.ephemeral_pool = ephemeral_pool
        self.pool = None
        self.merge = merge

    def _init_pool(self):
        if ThreadPool is None:
            return False
        if self.pool is not None:
            return True
        with self.lock:
            if self.pool is None:
                self.pool = ThreadPool(self.n_threads)
        return True

    def _clean_up_pool(self):
        if self.pool:
            self.pool.close()
            self.pool.terminate()
            self.pool = None

    def _fetch_usi(self, usi):
        use_pool = self._init_pool()
        agg = []
        if use_pool:
            with self.lock:
                for backend in self.backends.values():
                    result = self.pool.apply_async(backend.get, (usi, ))
                    agg.append((backend, result))
                tmp = []
                for backend, res in agg:
                    try:
                        res = res.get(self.timeout)
                        tmp.append((backend, res))
                    except (multiprocessing.TimeoutError, Exception) as err:
                        tmp.append((backend, err))
                agg = tmp
                if self.ephemeral_pool:
                    self._clean_up_pool()
        else:
            for backend in self.backends.values():
                try:
                    agg.append(backend, backend.get(usi))
                except Exception as err:
                    agg.append((backend, err))
                    continue
        return agg

    def coalesce(self, responses, method='first'):
        '''Merge responses from disparate servers into a single spectrum representation.

        The merging process will use the first of every array encountered, and all unique
        attributes.

        Parameters
        ----------
        responses : list
            A list of response values, pairs (:class:`_PROXIBackend` and either
            :class:`dict` or :class:`Exception`).
        method : str
            The name of the coalescence technique to use. Currently only "first" is
            supported.

        Returns
        -------
        result : :class:`dict`
            The coalesced spectrum
        '''
        if method not in self._coalesce_resolution_methods:
            raise ValueError("Coalescence method %r not recognized" % (method, ))

        def collapse_attribute(values):
            try:
                acc = list(set(v['value'] for v in values))
            except TypeError:
                acc = []
                for v in values:
                    if v['value'] not in acc:
                        acc.append(v['value'])

            result = []
            template = values[0].copy()
            for v in acc:
                t = template.copy()
                t['value'] = v
                result.append(t)
            return result

        arrays = {}
        attributes = defaultdict(list)

        found = []
        error = []

        for backend, response in responses:
            if isinstance(response, Exception):
                error.append((backend.name, (response)))
                continue
            else:
                found.append(backend.name)
            for array_name in ('m/z array', 'intensity array'):
                if array_name not in arrays:
                    arrays[array_name] = response[array_name]
                else:
                    array = response[array_name]
                    if len(array) != len(arrays[array_name]):
                        warnings.warn("Length mismatch from %s for %s" %
                            (backend.name, array_name))
                        arrays[array_name] = max((array, arrays[array_name]), key=len)
                    elif not np.allclose(array, arrays[array_name]):
                        warnings.warn("Value mismatch from %s for %s" %
                            (backend.name, array_name))
            for attr in response['attributes']:
                attributes[attr.get('accession', attr.get('name'))].append(attr)

        finalized_attributes = []
        for k, v in attributes.items():
            finalized_attributes.extend(collapse_attribute(v))

        result = {"responders": found, 'errors': error, 'attributes': finalized_attributes}
        result.update(arrays)
        if 'm/z array' not in result:
            raise ValueError("No valid responses found")
        return result

    def tag_with_source(self, responses):
        '''Mark each response with it's source.

        Parameters
        ----------
        responses : list
            A list of response values, pairs (:class:`_PROXIBackend` and either
            :class:`dict` or :class:`Exception`).

        Returns
        -------
        result : list[dict]
            The tagged :class:`dict` for each response.
        '''
        output = []
        for backend, response in responses:
            if isinstance(response, dict):
                response['source'] = backend
            else:
                response = {
                    "source": backend,
                    "error": response
                }
            output.append(response)
        return output

    def get(self, usi):
        '''Retrieve a ``USI`` from each PROXI service over the network.

        Parameters
        ----------
        usi : str or :class:`USI`
            The universal spectrum identifier to retrieve.

        Returns
        -------
        result : dict or list[dict]
            The spectrum coalesced from all responding PROXI hosts if :attr:`merge` is :const:`True`,
            or a list of responses marked by host.
        '''
        agg = self._fetch_usi(usi)
        if self.merge:
            return self.coalesce(agg)
        else:
            return self.tag_with_source(agg)

    def __call__(self, usi):
        return self.get(usi)

    def __del__(self):
        self._clean_up_pool()

_proxies = {
    "peptide_atlas": PeptideAtlasBackend,
    "massive": MassIVEBackend,
    "pride": PRIDEBackend,
    "jpost": JPOSTBackend,
    'proteome_exchange': ProteomeExchangeBackend,
}

default_backend = 'peptide_atlas'

AGGREGATOR_KEY = "aggregator"
AGGREGATOR = PROXIAggregator()


def proxi(usi, backend=default_backend, **kwargs):
    '''Retrieve a ``USI`` from a `PROXI <http://www.psidev.info/proxi>`.

    Parameters
    ----------
    usi : str or :class:`USI`
        The universal spectrum identifier to request.
    backend : str or :class:`Callable`
        Either the name of a PROXI host (peptide_atlas, massive, pride, jpost, or aggregator),
        or a callable object (which :class:`_PROXIBackend` instances are) which will be used
        to resolve the USI. The "aggregator" backend will use a :class:`PROXIAggregator` instance
        which will request the same USI from all the registered servers and attempt to merge their
        responses into a single whole. See :meth:`PROXIAggregator.coalesce` for more details on the
        merging process.
    **kwargs:
        extra arguments passed when constructing the backend by name.

    Returns
    -------
    dict :
        The spectrum as represented by the requested PROXI host.
    '''
    if isinstance(backend, str):
        if backend == AGGREGATOR_KEY:
            backend = AGGREGATOR
        elif backend in _proxies:
            backend = _proxies[backend](**kwargs)
        else:
            raise PyteomicsError("Unknown PROXI backend name: {}.".format(backend))
    elif isinstance(backend, type) and issubclass(backend, (_PROXIBackend, PROXIAggregator)):
        backend = backend(**kwargs)
    elif callable(backend):
        backend = backend
    else:
        raise TypeError("Unrecognized backend type: {0.__name__}".format(type(backend)))
    return backend(usi)
