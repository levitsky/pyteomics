from __future__ import absolute_import
import re
import operator as op
import math

try:
    basestring
except NameError:
    basestring = (str, bytes)

try:
    from collections import Container, Sized
except ImportError:
    from collections.abc import Container, Sized
from bisect import bisect_right
from contextlib import contextmanager


from .structures import PyteomicsError
from .file_helpers import _keepstate, IteratorContextManager, _make_chain
from .patch import pd


def _fix_docstring(f, **defaults):
    for argname, v in defaults.items():
        if v is not None:
            f.__doc__ = re.sub('{} : .*'.format(argname),
                               lambda m: m.group() + ', optional', f.__doc__)


def _calculate_qvalues(scores, isdecoy, peps=False, **kwargs):
    """Actual q-value calculation.

    Parameters
    ----------
    scores : numpy.ndarray
        Sorted array of PSMs.
    isdecoy : numpy.ndarray
        Sorted array of bools (decoy/target) or floats (PEPs).

    Returns
    -------
    out : numpy.ndarray
        Calculated q-values.
    """
    correction = kwargs.pop('correction', 0)
    ratio = kwargs.pop('ratio', 1)
    remove_decoy = kwargs.get('remove_decoy', False)
    formula = kwargs.pop('formula', (2, 1)[bool(remove_decoy)])
    if formula not in {1, 2}:
        raise PyteomicsError('`formula` must be either 1 or 2')

    # score_label = kwargs['score_label']
    cumsum = isdecoy.cumsum(dtype=np.float64)
    tfalse = cumsum.copy()
    ind = np.arange(1., scores.shape[0] + 1., dtype=np.float64)

    if peps:
        q = cumsum / ind
    else:
        if isinstance(correction, int):
            if correction == 1:
                tfalse += 1
            elif correction == 2:
                p = 1. / (1. + ratio)
                targ = ind - cumsum
                for i in range(tfalse.size):
                    tfalse[i] = _expectation(cumsum[i], targ[i], p)
        elif 0 < correction < 1:
            p = 1. / (1. + ratio)
            targ = ind - cumsum
            for i in range(tfalse.size):
                tfalse[i] = _confidence_value(
                    correction, cumsum[i], targ[i], p)
        elif correction:
            raise PyteomicsError('Invalid value for `correction`.')

        if formula == 1:
            q = tfalse / (ind - cumsum) / ratio
        else:
            q = (cumsum + tfalse / ratio) / ind

    # Make sure that q-values are equal for equal scores (conservatively)
    # and that q-values are monotonic
    for i in range(scores.size - 1, 0, -1):
        if (scores[i] == scores[i - 1] or
                q[i - 1] > q[i]):
            q[i - 1] = q[i]

    return q


def _qvalues_df(psms, keyf, isdecoy, **kwargs):
    full = kwargs.get('full_output', False)
    remove_decoy = kwargs.get('remove_decoy', False)
    peps = kwargs.get('pep')
    decoy_or_pep_label = _decoy_or_pep_label(**kwargs)
    q_label = kwargs.setdefault('q_label', 'q')
    score_label = kwargs.setdefault('score_label', 'score')
    if callable(keyf):
        keyf = psms.apply(keyf, axis=1)
    if callable(isdecoy):
        isdecoy = psms.apply(isdecoy, axis=1)
    if not isinstance(keyf, basestring):
        if psms.shape[0]:
            psms[score_label] = keyf
        else:
            psms[score_label] = []
        keyf = kwargs['score_label']
    if not isinstance(isdecoy, basestring):
        if psms.shape[0]:
            psms[decoy_or_pep_label] = isdecoy
        else:
            psms[decoy_or_pep_label] = []
        isdecoy = decoy_or_pep_label
    reverse = kwargs.get('reverse', False)

    if not full:  # create fields early
        if peps is None:
            fields = [(keyf, np.float64), (isdecoy, np.bool_),
                      (q_label, np.float64)]
        else:
            fields = [(isdecoy, np.float64), (q_label, np.float64)]
        dtype = np.dtype(fields)

    psms.sort_values([keyf, isdecoy], ascending=[
                     not reverse, True], inplace=True)

    if not psms.shape[0]:
        if full:
            psms[q_label] = []
            return psms
        else:
            return np.array([], dtype=dtype)

    q = _calculate_qvalues(psms[keyf].values, psms[
                           isdecoy].values, peps is not None, **kwargs)
    if remove_decoy:
        q = q[~psms[isdecoy].values]
        psms = psms[~psms[isdecoy]].copy()
    if not full:
        psms_ = np.empty_like(q, dtype=dtype)
        if peps is None:
            psms_[keyf] = psms[keyf]
        psms_[isdecoy] = psms[isdecoy]
        psms_[q_label] = q
        psms = psms_
    else:
        q_label = kwargs['q_label']
        psms[q_label] = q
    return psms


def _decoy_or_pep_label(**kwargs):
    peps = kwargs.get('pep')
    return kwargs.get('decoy_label', 'is decoy') if peps is None else kwargs.get(
        'pep_label', peps if isinstance(peps, basestring) else 'PEP')


def _construct_dtype(*args, **kwargs):
    full = kwargs.pop('full_output', False)
    peps = kwargs.get('pep')
    q_label = kwargs.setdefault('q_label', 'q')
    score_label = kwargs.setdefault('score_label', 'score')

    fields = [(score_label, np.float64),
              (_decoy_or_pep_label(**kwargs),
               np.bool_ if peps is None else np.float64),
              (q_label, np.float64)]
    # if all args are NumPy arrays with common dtype, use it in the output
    if full:
        dtypes = {getattr(arg, 'dtype', None) for arg in args}
        if len(dtypes) == 1 and None not in dtypes:
            psm_dtype = dtypes.pop()
        else:
            psm_dtype = np.object_
        dtype = np.dtype(fields + [('psm', psm_dtype)])
    else:
        dtype = np.dtype(fields)
    return dtype


def _make_qvalues(read, is_decoy_prefix, is_decoy_suffix, key):
    """Create a function that reads PSMs from a file and calculates q-values
    for each value of `key`."""

    def qvalues(*args, **kwargs):
        """Read `args` and return a NumPy array with scores and q-values.
        q-values are calculated either using TDA or based on provided values of PEP.

        Requires :py:mod:`numpy` (and optionally :py:mod:`pandas`).

        Parameters
        ----------

        positional args : file or str
            Files to read PSMs from. All positional arguments are treated as
            files. The rest of the arguments must be named.

        key : callable / array-like / iterable / str, keyword only
            If callable, a function used for sorting of PSMs. Should accept
            exactly one argument (PSM) and return a number (the smaller the better).
            If array-like, should contain scores for all given PSMs.
            If string, it is used as a field name (PSMs must be in a record array
            or a :py:class:`DataFrame`).

            .. warning::
                The default function may not work
                with your files, because format flavours are diverse.

        reverse : bool, keyword only, optional
            If :py:const:`True`, then PSMs are sorted in descending order,
            i.e. the value of the key function is higher for better PSMs.
            Default is :py:const:`False`.

        is_decoy : callable / array-like / iterable / str, keyword only
            If callable, a function used to determine if the PSM is decoy or not.
            Should accept exactly one argument (PSM) and return a truthy value if the
            PSM should be considered decoy.
            If array-like, should contain boolean values for all given PSMs.
            If string, it is used as a field name (PSMs must be in a record array
            or a :py:class:`DataFrame`).

            .. warning::
                The default function may not work
                with your files, because format flavours are diverse.

        decoy_prefix : str, optional
            If the default `is_decoy` function works for you, this parameter specifies which
            protein name prefix to use to detect decoy matches. If you provide your own
            `is_decoy`, or if you specify `decoy_suffix`, this parameter has no effect.
            Default is `"DECOY_"`.

        decoy_suffix : str, optional
            If the default `is_decoy` function works for you, this parameter specifies which
            protein name suffix to use to detect decoy matches. If you provide your own
            `is_decoy`, this parameter has no effect. Mutually exclusive with `decoy_prefix`.

        pep : callable / array-like / iterable / str, keyword only, optional
            If callable, a function used to determine the posterior error probability (PEP).
            Should accept exactly one argument (PSM) and return a float.
            If array-like, should contain float values for all given PSMs.
            If string, it is used as a field name (PSMs must be in a record array
            or a :py:class:`DataFrame`).

            .. note:: If this parameter is given, then PEP values will be used to calculate
               q-values. Otherwise, decoy PSMs will be used instead. This option conflicts with:
               `is_decoy`, `remove_decoy`, `formula`, `ratio`, `correction`.
               `key` can still be provided. Without `key`, PSMs will be sorted by PEP.

        remove_decoy : bool, keyword only, optional
            Defines whether decoy matches should be removed from the output.
            Default is :py:const:`False`.

            .. note:: If set to :py:const:`False`, then by default the decoy
               PSMs will be taken into account when estimating FDR. Refer to the
               documentation of :py:func:`fdr` for math; basically, if
               `remove_decoy` is :py:const:`True`, then formula 1 is used
               to control output FDR, otherwise it's formula 2. This can be
               changed by overriding the `formula` argument.

        formula : int, keyword only, optional
            Can be either 1 or 2, defines which formula should be used for FDR
            estimation. Default is 1 if `remove_decoy` is :py:const:`True`,
            else 2 (see :py:func:`fdr` for definitions).

        ratio : float, keyword only, optional
            The size ratio between the decoy and target databases. Default is
            1. In theory, the "size" of the database is the number of
            theoretical peptides eligible for assignment to spectra that are
            produced by *in silico* cleavage of that database.

        correction : int or float, keyword only, optional
            Possible values are 0, 1 and 2, or floating point numbers between 0 and 1.

            0 (default): no correction;

            1: enable "+1" correction. This accounts for the probability that a false
            positive scores better than the first excluded decoy PSM;

            2: this also corrects that probability for finite size of the sample,
            so the correction will be slightly less than "+1".

            If a floating point number
            is given, then instead of the expectation value for the number of false PSMs,
            the confidence value is used. The value of `correction` is then interpreted as
            desired confidence level. E.g., if correction=0.95, then the calculated q-values
            do not exceed the "real" q-values with 95% probability.

            See `this paper <http://dx.doi.org/10.1021/acs.jproteome.6b00144>`_ for further explanation.

        q_label : str, optional
            Field name for q-value in the output. Default is ``'q'``.

        score_label : str, optional
            Field name for score in the output. Default is ``'score'``.

        decoy_label : str, optional
            Field name for the decoy flag in the output. Default is ``'is decoy'``.

        pep_label : str, optional
            Field name for PEP in the output. Default is ``'PEP'``.

        full_output : bool, keyword only, optional
            If :py:const:`True`, then the returned array has PSM objects along
            with scores and q-values. Default is :py:const:`False`.

        **kwargs : passed to the :py:func:`chain` function.

        Returns
        -------
        out : numpy.ndarray
            A sorted array of records with the following fields:

            - 'score': :py:class:`np.float64`
            - 'is decoy': :py:class:`np.bool_`
            - 'q': :py:class:`np.float64`
            - 'psm': :py:class:`np.object_` (if `full_output` is :py:const:`True`)
        """
        import numpy as np

        @_keepstate
        def get_scores(*args, **kwargs):
            scores = []
            with read(*args, **kwargs) as f:
                for i, psm in enumerate(f):
                    row = []
                    for func in (keyf, isdecoy):
                        if callable(func):
                            row.append(func(psm))
                        elif isinstance(func, basestring):
                            row.append(psm[func])
                        else:
                            row.append(func[i])
                    row.append(None)
                    if full:
                        row.append(psm)
                    scores.append(tuple(row))
            return scores

        peps = kwargs.get('pep', None)
        if peps is not None:
            x = {'is_decoy', 'remove_decoy', 'formula',
                 'ratio', 'correction'}.intersection(kwargs)
            if x:
                raise PyteomicsError(
                    "Can't use these parameters with `pep`: " + ', '.join(x))
        keyf = kwargs.pop('key', key)
        reverse = kwargs.get('reverse', False)
        if keyf is None:
            keyf = peps
            if reverse:
                raise PyteomicsError(
                    'reverse = True when using PEPs for sorting')

        if not callable(keyf) and not isinstance(keyf, (Sized, Container)):
            keyf = np.array(list(keyf))

        if peps is None:
            if 'is_decoy' not in kwargs:
                if 'decoy_suffix' in kwargs:
                    isdecoy = lambda x: is_decoy_suffix(x, kwargs['decoy_suffix'])
                elif 'decoy_prefix' in kwargs:
                    isdecoy = lambda x: is_decoy_prefix(x, kwargs['decoy_prefix'])
                else:
                    isdecoy = is_decoy_prefix
            else:
                isdecoy = kwargs['is_decoy']
        else:
            isdecoy = peps

        if not callable(isdecoy) and not isinstance(isdecoy, (Sized, Container)):
            isdecoy = np.array(list(isdecoy))

        remove_decoy = kwargs.get('remove_decoy', False)
        decoy_or_pep_label = _decoy_or_pep_label(**kwargs)
        score_label = kwargs.setdefault('score_label', 'score')
        q_label = kwargs.setdefault('q_label', 'q')
        dtype = _construct_dtype(*args, **kwargs)
        full = kwargs.get('full_output', False)
        arr_flag = False
        psms = None

        # time to check arg type
        if pd is not None and all(isinstance(arg, pd.DataFrame) for arg in args):
            psms = pd.concat(args)
            return _qvalues_df(psms, keyf, isdecoy, **kwargs)

        if not all(isinstance(arg, np.ndarray) for arg in args):
            if isinstance(keyf, basestring):
                keyf = op.itemgetter(keyf)
            if isinstance(isdecoy, basestring):
                isdecoy = op.itemgetter(isdecoy)
            if isinstance(peps, basestring):
                peps = op.itemgetter(peps)

        if callable(keyf) or callable(isdecoy):
            kwargs.pop('full_output', None)
            scores = np.array(get_scores(*args, **kwargs), dtype=dtype)
        else:
            if all(isinstance(arg, np.ndarray) for arg in args):
                psms = np.concatenate(args)

            if not isinstance(keyf, basestring):
                keyf = np.array(keyf)
                arr_flag = True
            if not isinstance(isdecoy, basestring):
                isdecoy = np.array(isdecoy)
                arr_flag = True

            if arr_flag:
                scores = np.empty(keyf.size if hasattr(
                    keyf, 'size') else isdecoy.size, dtype=dtype)
                for func, label in zip((keyf, isdecoy), (score_label, decoy_or_pep_label)):
                    if not isinstance(func, basestring):
                        scores[label] = func
                    else:
                        scores[label] = psms[func]
            else:
                scores = np.empty(psms.shape[0], dtype=dtype)
                scores[score_label] = psms[keyf]
                scores[decoy_or_pep_label] = psms[isdecoy]

        if not scores.size:
            if full and psms is not None:
                return psms
            return scores

        if not reverse:
            keys = scores[decoy_or_pep_label], scores[score_label]
        else:
            keys = scores[decoy_or_pep_label], -scores[score_label]
        lexsort = np.lexsort(keys)
        scores = scores[lexsort]
        if psms is not None:
            psms = psms[lexsort]

        scores[q_label] = _calculate_qvalues(scores[score_label], scores[
                                             decoy_or_pep_label], peps is not None, **kwargs)
        if remove_decoy:
            if psms is not None:
                psms = psms[~scores[decoy_or_pep_label]]
            scores = scores[~scores[decoy_or_pep_label]]

        if full and psms is not None:
            if isinstance(psms, np.ndarray):
                fields = sorted(psms.dtype.fields,
                                key=lambda x: psms.dtype.fields[x][1])
                extra = []
                for func, label in zip((keyf, isdecoy), ('score', decoy_or_pep_label)):
                    if not (isinstance(func, basestring) or label in psms.dtype.fields):
                        extra.append(label)
                    elif label in psms.dtype.fields:
                        psms[label] = scores[label]
                newdt = [(name, psms.dtype.fields[name][0]) for name in fields] + [
                    (name, np.float64) for name in extra] + [(q_label, np.float64)]
                psms_ = psms
                psms = np.empty_like(psms_, dtype=newdt)
                for f in fields:
                    psms[f] = psms_[f]
                for f in extra:
                    psms[f] = scores[f]
            else:
                for func, label in zip((keyf, isdecoy), ('score', decoy_or_pep_label)):
                    if not isinstance(label, basestring):
                        psms[label] = scores[label]
            psms[q_label] = scores[q_label]
            return psms
        return scores

    _fix_docstring(qvalues, is_decoy=is_decoy_prefix, key=key)
    if read is _iter:
        qvalues.__doc__ = qvalues.__doc__.replace("""positional args : file or str
            Files to read PSMs from. All positional arguments are treated as
            files.""", """positional args : iterables
            Iterables to read PSMs from. All positional arguments are chained."""
                ).replace("""\n            .. warning::
                The default function may not work
                with your files, because format flavours are diverse.

        decoy_prefix : str, optional
            If the default `is_decoy` function works for you, this parameter specifies which
            protein name prefix to use to detect decoy matches. If you provide your own
            `is_decoy`, or if you specify `decoy_suffix`, this parameter has no effect.
            Default is `"DECOY_"`.

        decoy_suffix : str, optional
            If the default `is_decoy` function works for you, this parameter specifies which
            protein name suffix to use to detect decoy matches. If you provide your own
            `is_decoy`, this parameter has no effect. Mutually exclusive with `decoy_prefix`.\n""", "")

    return qvalues


def _make_filter(read, is_decoy_prefix, is_decoy_suffix, key, qvalues):
    """Create a function that reads PSMs from a file and filters them to
    the desired FDR level (estimated by TDA), returning the top PSMs
    sorted by `key`.
    """
    def filter(*args, **kwargs):
        try:
            fdr = kwargs.pop('fdr')
        except KeyError:
            raise PyteomicsError('Keyword argument required: fdr')

        args = [list(arg) if not isinstance(
            arg, (Container, Sized)) else arg for arg in args]
        peps = kwargs.get('pep')
        if peps is None:
            remove_decoy = kwargs.pop('remove_decoy', True)
            scores = qvalues(*args, remove_decoy=remove_decoy, **kwargs)
        else:
            scores = qvalues(*args, **kwargs)
        keyf = kwargs.pop('key', key)
        if keyf is None:
            keyf = peps
        reverse = kwargs.pop('reverse', False)
        better = [op.lt, op.gt][bool(reverse)]
        if 'is_decoy' not in kwargs:
            if 'decoy_suffix' in kwargs:
                isdecoy = lambda x: is_decoy_suffix(x, kwargs['decoy_suffix'])
            elif 'decoy_prefix' in kwargs:
                isdecoy = lambda x: is_decoy_prefix(x, kwargs['decoy_prefix'])
            else:
                isdecoy = is_decoy_prefix
        else:
            isdecoy = kwargs['is_decoy']
        kwargs.pop('formula', None)
        decoy_or_pep_label = _decoy_or_pep_label(**kwargs)
        score_label = kwargs.setdefault('score_label', 'score')
        q_label = kwargs.get('q_label', 'q')

        try:
            i = scores[q_label].searchsorted(fdr, side='right')
            if isinstance(i, Sized):
                i = i[0]
        except AttributeError:
            i = bisect_right(scores['q'], fdr)
        if kwargs.pop('full_output', False):
            if pd is not None and isinstance(scores, pd.DataFrame):
                return scores.iloc[:i]
            elif callable(keyf) or callable(isdecoy):
                return scores['psm'][:i]
            else:
                return scores[:i]
        elif not scores.size:
            return (_ for _ in ())
        if peps is None:
            label = score_label
        else:
            label = decoy_or_pep_label
        cutoff = scores[label][i] if i < scores.size else (
            scores[label][-1] + (1, -1)[bool(reverse)])

        def out():
            with read(*args, **kwargs) as f:
                for p, s in zip(f, scores):
                    if peps is not None or not remove_decoy or not s[decoy_or_pep_label]:
                        if better(s[label], cutoff):
                            yield p
        return out()

    def _filter(*args, **kwargs):
        """Read `args` and yield only the PSMs that form a set with
        estimated false discovery rate (FDR) not exceeding `fdr`.

        Requires :py:mod:`numpy` and, optionally, :py:mod:`pandas`.

        Parameters
        ----------
        positional args : file or str
            Files to read PSMs from. All positional arguments are treated as
            files. The rest of the arguments must be named.

        fdr : float, keyword only, 0 <= fdr <= 1
            Desired FDR level.

        key : callable / array-like / iterable / str, keyword only
            A function used for sorting of PSMs. Should accept exactly one
            argument (PSM) and return a number (the smaller the better). The
            default is a function that tries to extract e-value from the PSM.

            .. warning::
                The default function may not work
                with your files, because format flavours are diverse.

        reverse : bool, keyword only, optional
            If :py:const:`True`, then PSMs are sorted in descending order,
            i.e. the value of the key function is higher for better PSMs.
            Default is :py:const:`False`.

        is_decoy : callable / array-like / iterable / str, keyword only
            A function used to determine if the PSM is decoy or not. Should
            accept exactly one argument (PSM) and return a truthy value if the
            PSM should be considered decoy.

            .. warning::
                The default function may not work
                with your files, because format flavours are diverse.

        decoy_prefix : str, optional
            If the default `is_decoy` function works for you, this parameter specifies which
            protein name prefix to use to detect decoy matches. If you provide your own
            `is_decoy`, or if you specify `decoy_suffix`, this parameter has no effect.
            Default is `"DECOY_"`.

        decoy_suffix : str, optional
            If the default `is_decoy` function works for you, this parameter specifies which
            protein name suffix to use to detect decoy matches. If you provide your own
            `is_decoy`, this parameter has no effect. Mutually exclusive with `decoy_prefix`.

        remove_decoy : bool, keyword only, optional
            Defines whether decoy matches should be removed from the output.
            Default is :py:const:`True`.

            .. note:: If set to :py:const:`False`, then by default the decoy
               PSMs will be taken into account when estimating FDR. Refer to the
               documentation of :py:func:`fdr` for math; basically, if
               `remove_decoy` is :py:const:`True`, then formula 1 is used
               to control output FDR, otherwise it's formula 2. This can be
               changed by overriding the `formula` argument.

        formula : int, keyword only, optional
            Can be either 1 or 2, defines which formula should be used for FDR
            estimation. Default is 1 if `remove_decoy` is :py:const:`True`,
            else 2 (see :py:func:`fdr` for definitions).

        ratio : float, keyword only, optional
            The size ratio between the decoy and target databases. Default is
            1. In theory, the "size" of the database is the number of
            theoretical peptides eligible for assignment to spectra that are
            produced by *in silico* cleavage of that database.

        correction : int or float, keyword only, optional
            Possible values are 0, 1 and 2, or floating point numbers between 0 and 1.

            0 (default): no correction;

            1: enable "+1" correction. This accounts for the probability that a false
            positive scores better than the first excluded decoy PSM;

            2: this also corrects that probability for finite size of the sample,
            so the correction will be slightly less than "+1".

            If a floating point number
            is given, then instead of the expectation value for the number of false PSMs,
            the confidence value is used. The value of `correction` is then interpreted as
            desired confidence level. E.g., if correction=0.95, then the calculated q-values
            do not exceed the "real" q-values with 95% probability.

            See `this paper <http://dx.doi.org/10.1021/acs.jproteome.6b00144>`_ for further explanation.

        pep : callable / array-like / iterable / str, keyword only, optional
            If callable, a function used to determine the posterior error probability (PEP).
            Should accept exactly one argument (PSM) and return a float.
            If array-like, should contain float values for all given PSMs.
            If string, it is used as a field name (PSMs must be in a record array
            or a :py:class:`DataFrame`).

            .. note:: If this parameter is given, then PEP values will be used to calculate
               q-values. Otherwise, decoy PSMs will be used instead. This option conflicts with:
               `is_decoy`, `remove_decoy`, `formula`, `ratio`, `correction`.
               `key` can still be provided. Without `key`, PSMs will be sorted by PEP.

        full_output : bool, keyword only, optional
            If :py:const:`True`, then an array of PSM objects is returned.
            Otherwise, an iterator / context manager object is returned, and the
            files are parsed twice. This saves some RAM, but is ~2x slower.
            Default is :py:const:`True`.

            .. note:: The name for the parameter comes from the fact that it is
                      internally passed to :py:func:`qvalues`.

        q_label : str, optional
            Field name for q-value in the output. Default is ``'q'``.

        score_label : str, optional
            Field name for score in the output. Default is ``'score'``.

        decoy_label : str, optional
            Field name for the decoy flag in the output. Default is ``'is decoy'``.

        pep_label : str, optional
            Field name for PEP in the output. Default is ``'PEP'``.

        **kwargs : passed to the :py:func:`chain` function.

        Returns
        -------
        out : iterator or :py:class:`numpy.ndarray` or :py:class:`pandas.DataFrame`
        """
        if kwargs.pop('full_output', True):
            return filter(*args, full_output=True, **kwargs)
        return IteratorContextManager(filter, *args, **kwargs)

    _fix_docstring(_filter, is_decoy=is_decoy_prefix, key=key)
    if read is _iter:
        _filter.__doc__ = _filter.__doc__.replace("""positional args : file or str
            Files to read PSMs from. All positional arguments are treated as
            files.""", """positional args : iterables
            Iterables to read PSMs from. All positional arguments are chained.""").replace(
                """\n            .. warning::
                The default function may not work
                with your files, because format flavours are diverse.

        decoy_prefix : str, optional
            If the default `is_decoy` function works for you, this parameter specifies which
            protein name prefix to use to detect decoy matches. If you provide your own
            `is_decoy`, or if you specify `decoy_suffix`, this parameter has no effect.
            Default is `"DECOY_"`.

        decoy_suffix : str, optional
            If the default `is_decoy` function works for you, this parameter specifies which
            protein name suffix to use to detect decoy matches. If you provide your own
            `is_decoy`, this parameter has no effect. Mutually exclusive with `decoy_prefix`.\n""", "")
    return _filter


@contextmanager
def _itercontext(x, **kw):
    try:
        yield (row for i, row in x.iterrows())
    except AttributeError:
        yield x


_iter = _make_chain(_itercontext, 'iter')
qvalues = _make_qvalues(_iter, None, None, None)

filter = _make_filter(_iter, None, None, None, qvalues)
filter.chain = _make_chain(filter, 'filter', True)

try:
    import numpy as np
    _precalc_fact = np.log([math.factorial(n) for n in range(20)])

    def log_factorial(x):
        x = np.array(x)
        pf = _precalc_fact
        m = (x >= pf.size)
        out = np.empty(x.shape)
        out[~m] = pf[x[~m].astype(int)]
        x = x[m]
        out[m] = x * np.log(x) - x + 0.5 * np.log(2 * np.pi * x)
        return out

    def _expectation(d, T, p=0.5):
        if T is None:
            return d + 1
        T = np.array(T, dtype=int)
        m = np.arange(T.max() + 1, dtype=int)
        pi = np.exp(_log_pi(d, m, p))
        return ((m * pi).cumsum() / pi.cumsum())[T]

    def _confidence_value(conf, d, T, p=0.5):
        if T is not None:
            T = np.array(T, dtype=int)
            m = np.arange(T.max() + 1, dtype=int)
        else:
            m = np.arange(max(50 * d, 10000))
        log_pi = _log_pi(d, m, p)
        pics = np.exp(log_pi).cumsum()
        return np.searchsorted(pics, conf * (pics[T] if T is not None else 1))

except ImportError:
    def log_factorial(n):
        if n > 10:
            return n * math.log(n) - n + 0.5 * math.log(2 * math.pi * n)
        else:
            return math.log(math.factorial(n))

    def _expectation(*a, **k):
        raise NotImplementedError('NumPy required')

    def _confidence_value(*a, **k):
        raise NotImplementedError('NumPy required')


def _log_pi_r(d, k, p=0.5):
    return k * math.log(p) + log_factorial(k + d) - log_factorial(k) - log_factorial(d)


def _log_pi(d, k, p=0.5):
    return _log_pi_r(d, k, p) + (d + 1) * math.log(1 - p)


def _make_fdr(is_decoy_prefix, is_decoy_suffix):
    def fdr(psms=None, formula=1, is_decoy=None, ratio=1, correction=0, pep=None, decoy_prefix='DECOY_', decoy_suffix=None):
        """Estimate FDR of a data set using TDA or given PEP values.
        Two formulas can be used. The first one (default) is:

        .. math::

                FDR = \\frac{N_{decoy}}{N_{target} * ratio}

        The second formula is:

        .. math::

                FDR = \\frac{N_{decoy} * (1 + \\frac{1}{ratio})}{N_{total}}

        .. note::
            This function is less versatile than :py:func:`qvalues`. To obtain FDR,
            you can call :py:func:`qvalues` and take the last q-value. This function
            can be used (with `correction = 0` or `1`) when :py:mod:`numpy` is not available.

        Parameters
        ----------
        psms : iterable, optional
            An iterable of PSMs, e.g. as returned by :py:func:`read`.
            Not needed if `is_decoy` is an iterable.

        formula : int, optional
            Can be either 1 or 2, defines which formula should be used for FDR
            estimation. Default is 1.

        is_decoy : callable, iterable, or str
            If callable, should accept exactly one argument (PSM) and return a truthy value
            if the PSM is considered decoy. Default is :py:func:`is_decoy`.
            If array-like, should contain float values for all given PSMs.
            If string, it is used as a field name (PSMs must be in a record array
            or a :py:class:`pandas.DataFrame`).

            .. warning::
                The default function may not work
                with your files, because format flavours are diverse.

        decoy_prefix : str, optional
            If the default `is_decoy` function works for you, this parameter specifies which
            protein name prefix to use to detect decoy matches. If you provide your own
            `is_decoy`, or if you specify `decoy_suffix`, this parameter has no effect.
            Default is `"DECOY_"`.

        decoy_suffix : str, optional
            If the default `is_decoy` function works for you, this parameter specifies which
            protein name suffix to use to detect decoy matches. If you provide your own
            `is_decoy`, this parameter has no effect. Mutually exclusive with `decoy_prefix`.

        pep : callable, iterable, or str, optional
            If callable, a function used to determine the posterior error probability (PEP).
            Should accept exactly one argument (PSM) and return a float.
            If array-like, should contain float values for all given PSMs.
            If string, it is used as a field name (PSMs must be in a record array
            or a :py:class:`pandas.DataFrame`).

            .. note:: If this parameter is given, then PEP values will be used to calculate FDR.
               Otherwise, decoy PSMs will be used instead. This option conflicts with:
               `is_decoy`, `formula`, `ratio`, `correction`.

        ratio : float, optional
            The size ratio between the decoy and target databases. Default is 1.
            In theory, the "size" of the database is the number of
            theoretical peptides eligible for assignment to spectra that are
            produced by *in silico* cleavage of that database.

        correction : int or float, optional
            Possible values are 0, 1 and 2, or floating point numbers between 0 and 1.

            0 (default): no correction;

            1: enable "+1" correction. This accounts for the probability that a false
            positive scores better than the first excluded decoy PSM;

            2: this also corrects that probability for finite size of the sample,
            so the correction will be slightly less than "+1".

            If a floating point number
            is given, then instead of the expectation value for the number of false PSMs,
            the confidence value is used. The value of `correction` is then interpreted as
            desired confidence level. E.g., if correction=0.95, then the calculated q-values
            do not exceed the "real" q-values with 95% probability.

            See `this paper <http://dx.doi.org/10.1021/acs.jproteome.6b00144>`_ for further explanation.

            .. note::
                Requires :py:mod:`numpy`, if `correction` is a float or 2.

            .. note::
                Correction is only needed if the PSM set at hand was obtained using TDA
                filtering based on decoy counting (as done by using :py:func:`!filter` without
                `correction`).

        Returns
        -------
        out : float
            The estimation of FDR, (roughly) between 0 and 1.
        """
        if formula not in {1, 2}:
            raise PyteomicsError('`formula` must be either 1 or 2.')
        total, decoy = 0, 0
        if pep is not None:
            is_decoy = pep
        elif is_decoy is None:
            if decoy_suffix is not None:
                is_decoy = lambda x: is_decoy_suffix(x, decoy_suffix)
            else:
                is_decoy = lambda x: is_decoy_prefix(x, decoy_prefix)
        if isinstance(is_decoy, basestring):
            decoy = psms[is_decoy].sum()
            total = psms.shape[0]
        elif callable(is_decoy):
            for psm in psms:
                total += 1
                d = is_decoy(psm)
                decoy += d if pep is not None else bool(d)
        else:
            if not isinstance(is_decoy, (Sized, Container)):
                is_decoy = list(is_decoy)
            if pep is not None:
                decoy = sum(is_decoy)
            else:
                decoy = sum(map(bool, is_decoy))
            total = len(is_decoy)
        if pep is not None:
            return float(decoy) / total
        tfalse = decoy
        if correction == 1 or (correction == 2 and total / decoy > 10):
            tfalse += 1
        elif correction == 2:
            p = 1. / (1. + ratio)
            tfalse = _expectation(decoy, total - decoy, p)
        elif 0 < correction < 1:
            p = 1. / (1. + ratio)
            tfalse = _confidence_value(correction, decoy, total - decoy, p)
        if formula == 1:
            return float(tfalse) / (total - decoy) / ratio
        return (decoy + tfalse / ratio) / total

    _fix_docstring(fdr, is_decoy=is_decoy_prefix)
    if is_decoy_prefix is None:
        fdr.__doc__ = fdr.__doc__.replace(
            """\n            .. warning::
                The default function may not work
                with your files, because format flavours are diverse.

        decoy_prefix : str, optional
            If the default `is_decoy` function works for you, this parameter specifies which
            protein name prefix to use to detect decoy matches. If you provide your own
            `is_decoy`, or if you specify `decoy_suffix`, this parameter has no effect.
            Default is `"DECOY_"`.

        decoy_suffix : str, optional
            If the default `is_decoy` function works for you, this parameter specifies which
            protein name suffix to use to detect decoy matches. If you provide your own
            `is_decoy`, this parameter has no effect. Mutually exclusive with `decoy_prefix`.\n""", "")
    return fdr


fdr = _make_fdr(None, None)
