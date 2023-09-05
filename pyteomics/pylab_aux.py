"""
pylab_aux - auxiliary functions for plotting with pylab
=======================================================

This module serves as a collection of useful routines for data plotting with
matplotlib.

Generic plotting
----------------

  :py:func:`plot_line` - plot a line.

  :py:func:`scatter_trend` - plot a scatter plot with a regression line.

  :py:func:`plot_function_3d` - plot a 3D graph of a function of two variables.

  :py:func:`plot_function_contour` - plot a contour graph of a function of
  two variables.

Spectrum visualization
----------------------

  :py:func:`plot_spectrum` - plot a single spectrum (m/z vs intensity).

  :py:func:`annotate_spectrum` - plot and annotate peaks in MS/MS spectrum.

  :py:func:`mirror` - create a mirror plot of two spectra (using :py:mod:`spectrum_utils`).

FDR control
-----------

  :py:func:`plot_qvalue_curve` - plot the dependence of q-value on the amount of PSMs
  (similar to a ROC curve).

See also
--------

  - `Matplotlib cookbook <http://www.scipy.org/Cookbook/Matplotlib/>`_
  - `Matplotlib tutorial
    <http://matplotlib.sourceforge.net/mpl_toolkits/mplot3d/tutorial.html>`_

Dependencies
------------

This module requires :py:mod:`matplotlib`. Optional dependencies: :py:mod:`adjustText`, :py:mod:`spectrum_utils`.

-------------------------------------------------------------------------------

"""

#   Copyright 2012 Anton Goloborodko, Lev Levitsky
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

import pylab
import numpy as np
from .auxiliary import linear_regression, PyteomicsError
from .version import VersionInfo
from . import parser, mass, mgf, proforma

try:
    import spectrum_utils
    if VersionInfo(spectrum_utils.__version__) < VersionInfo('0.4'):
        raise ImportError("Supported spectrum_utils version is 0.4.0 or newer.")
    import spectrum_utils.spectrum as sus
    import spectrum_utils.plot as sup
except ImportError:
    sus = sup = None


def plot_line(a, b, xlim=None, *args, **kwargs):
    """Plot a line y = a * x + b.

    Parameters
    ----------
    a : float
        The slope of the line.
    b : float
        The intercept of the line.
    xlim : tuple, optional
        Minimal and maximal values of `x`. If not given, :py:func:`pylab.xlim` will be called.
    *args
        Passed to :py:func:`pylab.plot` after `x` and `y` values.
    **kwargs
        Passed to :py:func:`pylab.plot`.

    Returns
    -------
    out : matplotlib.lines.Line2D
        The line object.
    """
    if xlim is None:
        xlim = pylab.xlim()
    return pylab.plot([xlim[0], xlim[1]], [a * xlim[0] + b, a * xlim[1] + b], *args, **kwargs)


def scatter_trend(x, y=None, **kwargs):
    """Make a scatter plot with a linear regression.

    Parameters
    ----------
    x : array_like of float
        1-D array of floats. If `y` is omitted, `x` must be a 2-D array of shape (N, 2).
    y : array_like of float, optional
        1-D arrays of floats. If `y` is omitted or :py:const:`None`, `x` must be a 2-D array of shape (N, 2).
    plot_trend : bool, optional
        If :py:const:`True` then plot a trendline (default).
    plot_sigmas : bool, optional
        If :py:const:`True` then plot confidence intervals of the linear fit.
        :py:const:`False` by default.
    show_legend : bool, optional
        If :py:const:`True`, a legend will be shown with linear fit equation,
        correlation coefficient, and standard deviation from the fit. Default is
        :py:const:`True`.
    title : str, optional
        The title. Empty by default.
    xlabel, ylabel : str, optional
        The axes labels. Empty by default.
    alpha_legend : float, optional
        Legend box transparency. 1.0 by default
    scatter_kwargs : dict, optional
        Keyword arguments for :py:func:`pylab.scatter`.
        Empty by default.
    plot_kwargs : dict, optional
        Keyword arguments for :py:func:`plot_line`.
        By default, sets `xlim` and `label`.
    legend_kwargs : dict, optional
        Keyword arguments for :py:func:`pylab.legend`.
        Default is :py:const:`{'loc': 'upper left'}`.
    sigma_kwargs : dict, optional
        Keyword arguments for :py:func:`pylab.plot` used for sigma lines.
        Default is :py:const:`{'color': 'red', 'linestyle': 'dashed'}`.
    sigma_values : iterable, optional
        Each value will be multiplied with standard error of the fit, and the line
        shifted by the resulting value will be plotted. Default is :py:const:`range(-3, 4)`.
    regression : callable, optional
        Function to perform linear regression. Will be given ``x`` and ``y`` as arguments.
        Must return a 4-tuple: (a, b, r, stderr).
        Default is :py:func:`pyteomics.auxiliary.linear_regression`.

    Returns
    -------
    out : tuple
        A (scatter_plot, trend_line, sigma_lines, legend) tuple.
    """
    regression = kwargs.get('regression', linear_regression)
    a, b, r, stderr = regression(x, y)
    pylab.title(kwargs.get('title', ''))
    pylab.xlabel(kwargs.get('xlabel', ''))
    pylab.ylabel(kwargs.get('ylabel', ''))

    equation = (
        '$y\,=\,{:.3f}x\,{}\,{:.3f}$, '
        '$R^2=\,{:.3f}$ \n$\sigma\,=\,{:.3f}$'.format(
            a, '-' if b < 0 else '+', abs(b), r*r, stderr))

    if y is None:
        x = np.array(x, copy=False)
        y = x[:, 1]
        x = x[:, 0]
    else:
        x = np.array(x)
        y = np.array(y)
    sc = pylab.scatter(x, y, **kwargs.get('scatter_kwargs', {}))
    xlim = (x.min(), x.max())
    plkw = kwargs.get('plot_kwargs', {}).copy()
    plkw.setdefault('xlim', xlim)
    plkw.setdefault('label', equation)
    if kwargs.get('plot_trend', True):
        line = plot_line(a, b, **plkw)
    else:
        line = None

    if kwargs.get('plot_sigmas', False):
        s_lines = []
        sigma_kwargs = kwargs.get('sigma_kwargs', {'color': 'red', 'linestyle': 'dashed'})
        for i in kwargs.get('sigma_values', range(-3, 4)):
            s_lines.append(plot_line(a, b + i * stderr, xlim, **sigma_kwargs))
    else:
        s_lines = None

    if kwargs.get('show_legend', True):
        legend = pylab.legend(**kwargs.get('legend_kwargs', {'loc': 'upper left'}))
        legend_frame = legend.get_frame()
        legend_frame.set_alpha(kwargs.get('alpha_legend', 1.0))
    else:
        legend = None
    return sc, line, s_lines, legend


def plot_function_3d(x, y, function, **kwargs):
    """Plot values of a function of two variables in 3D.

    More on 3D plotting in pylab:

    http://www.scipy.org/Cookbook/Matplotlib/mplot3D

    Parameters
    ----------
    x : array_like of float
        The plotting range on X axis.
    y : array_like of float
        The plotting range on Y axis.
    function : function
        The function to plot.
    plot_type : {'surface', 'wireframe', 'scatter', 'contour', 'contourf'}, keyword only, optional
        The type of a plot, see
        `scipy cookbook <http://www.scipy.org/Cookbook/Matplotlib/mplot3D>`_
        for examples. The default value is 'surface'.
    num_contours : int
        The number of contours to plot, 50 by default.
    xlabel : str, keyword only, optional
        The X axis label. Empty by default.
    ylabel : str, keyword only, optional
        The Y axis label. Empty by default.
    zlabel : str, keyword only, optional
        The Z axis label. Empty by default.
    title : str, keyword only, optional
        The title. Empty by default.
    **kwargs
        Passed to the respective plotting function.
    """
    import mpl_toolkits.mplot3d.axes3d as pylab3d
    ax = pylab3d.Axes3D(pylab.gcf())
    ax.set_xlabel(kwargs.pop('xlabel', ''))
    ax.set_ylabel(kwargs.pop('ylabel', ''))
    ax.set_zlabel(kwargs.pop('zlabel', ''))
    ax.set_title(kwargs.pop('title', ''))
    X, Y = np.meshgrid(x, y)
    Z = []
    for y_value in y:
        Z.append([])
        for x_value in x:
            Z[-1].append(function(x_value, y_value))
    Z = np.array(Z)
    plot_type = kwargs.pop('plot_type', 'surface')
    if plot_type == 'surface':
        ax.plot_surface(X, Y, Z,
                rstride=kwargs.pop('rstride', 1),
                cstride=kwargs.pop('cstride', 1),
                cmap=kwargs.pop('cmap', pylab.cm.jet),
                **kwargs)
    elif plot_type == 'wireframe':
        ax.plot_wireframe(X, Y, Z,
                cmap=kwargs.pop('cmap', pylab.cm.jet), **kwargs)
    elif plot_type == 'scatter':
        ax.scatter3D(np.ravel(X), np.ravel(Y), np.ravel(Z), **kwargs)
    elif plot_type == 'contour':
        num_contours = kwargs.pop('num_contours', 50)
        ax.contour3D(X, Y, Z, num_contours,
                cmap=kwargs.pop('cmap', pylab.cm.jet), **kwargs)
    elif plot_type == 'contourf':
        num_contours = kwargs.pop('num_contours', 50)
        ax.contourf3D(X, Y, Z, num_contours,
                cmap=kwargs.pop('cmap', pylab.cm.jet), **kwargs)
    else:
        raise PyteomicsError('Unknown plot type: {}'.format(plot_type))


def plot_function_contour(x, y, function, **kwargs):
    """Make a contour plot of a function of two variables.

    Parameters
    ----------
    x, y : array_like of float
        The positions of the nodes of a plotting grid.
    function : function
        The function to plot.
    filling : bool
        Fill contours if True (default).
    num_contours : int
        The number of contours to plot, 50 by default.
    xlabel, ylabel : str, optional
        The axes labels. Empty by default.
    title : str, optional
        The title. Empty by default.
    **kwargs
        Passed to :py:func:`pylab.contour` or :py:func:`pylab.contourf`.
    """
    pylab.xlabel(kwargs.pop('xlabel', ''))
    pylab.ylabel(kwargs.pop('ylabel', ''))
    pylab.title(kwargs.pop('title', ''))
    X, Y = np.meshgrid(x, y)
    Z = []
    for y_value in y:
        Z.append([])
        for x_value in x:
            Z[-1].append(function(x_value, y_value))
    Z = np.array(Z)
    num_contours = kwargs.pop('num_contours', 50)
    if kwargs.pop('filling', True):
        pylab.contourf(X, Y, Z, num_contours,
                cmap=kwargs.pop('cmap', pylab.cm.jet), **kwargs)
    else:
        pylab.contour(X, Y, Z, num_contours,
                cmap=kwargs.pop('cmap', pylab.cm.jet), **kwargs)


def plot_qvalue_curve(qvalues, *args, **kwargs):
    """
    Plot a curve with q-values on the X axis and corresponding PSM number
    (starting with ``1``) on the Y axis.

    Parameters
    ----------
    qvalues : array-like
        An array of q-values for sorted PSMs.
    xlabel : str, keyword only, optional
        Label for the X axis. Default is "q-value".
    ylabel : str, keyword only, optional
        Label for the Y axis. Default is "# of PSMs".
    title : str, keyword only, optional
        The title. Empty by default.
    *args
        Given to :py:func:`pylab.plot` after `x` and `y`.
    **kwargs
        Given to :py:func:`pylab.plot`.

    Returns
    -------
    out : matplotlib.lines.Line2D
    """
    pylab.xlabel(kwargs.pop('xlabel', 'q-value'))
    pylab.ylabel(kwargs.pop('ylabel', '# of PSMs'))
    pylab.title(kwargs.pop('title', ''))
    return pylab.plot(qvalues, 1 + np.arange(qvalues.size), *args, **kwargs)


def _default_plot_spectrum(spectrum, *args, **kwargs):
    ax = kwargs.pop('ax', None) or pylab.gca()
    if kwargs.pop('centroided', True):
        kwargs.setdefault('align', 'center')
        kwargs.setdefault('width', 0)
        kwargs.setdefault('linewidth', 1)
        kwargs.setdefault('edgecolor', 'k')
        ax.bar(spectrum['m/z array'], spectrum['intensity array'], *args, **kwargs)
    else:
        ax.plot(spectrum['m/z array'], spectrum['intensity array'], *args, **kwargs)
    return ax


def _spectrum_utils_plot(spectrum, *args, **kwargs):

    with SpectrumUtilsColorScheme(kwargs.pop('colors', None)):
        spectrum = _spectrum_utils_create_spectrum(spectrum, None, *args, **kwargs)
        return sup.spectrum(spectrum)


def _spectrum_utils_iplot(spectrum, *args, **kwargs):
    import spectrum_utils.iplot as supi
    with SpectrumUtilsColorScheme(kwargs.pop('colors', None)):
        spectrum = _spectrum_utils_create_spectrum(spectrum, None, *args, **kwargs)
        return supi.spectrum(spectrum)


_plot_backends = {
    'default': _default_plot_spectrum,
    'spectrum_utils': _spectrum_utils_plot,
    'spectrum_utils.iplot': _spectrum_utils_iplot,
}


def plot_spectrum(spectrum, *args, **kwargs):
    """
    Plot a spectrum, assuming it is a dictionary containing "m/z array" and "intensity array".

    Parameters
    ----------
    spectrum : dict
        A dictionary, as returned by pyteomics MS data parsers.
        Must contain "m/z array" and "intensity array" keys with decoded arrays.
    backend : str, keyword only, optional
        One of `{'default', 'spectrum_utils', 'spectrum_utils.iplot'}`.
        The `spectrum_utils` backend requires installing :py:mod:`spectrum_utils`.
        The `spectrum_utils.iplot` backend requires installing :py:mod:`spectrum_utils[iplot]`.
    xlabel : str, keyword only, optional
        Label for the X axis. Default is "m/z".
    ylabel : str, keyword only, optional
        Label for the Y axis. Default is "intensity".
    title : str, keyword only, optional
        The title. Empty by default.

    centroided : bool, keyword only, optional
        Works only for the `default` backend.
        If :py:const:`True` (default), peaks of the spectrum are plotted using :py:func:`pylab.bar`.
        If :py:const:`False`, the arrays are simply plotted using :py:func:`pylab.plot`.
    *args
        When using `default` backend: given to :py:func:`pylab.plot` or :py:func:`pylab.bar` (depending on `centroided`).
    **kwargs
        When using `default` backend: given to :py:func:`pylab.plot` or :py:func:`pylab.bar` (depending on `centroided`).

    min_intensity : float, keyword only, optional
        Remove low-intensity peaks; this is a factor of maximum peak intensity. Default is 0 (no filtering).
        Only works with `spectrum_utils` and `spectrum_utils.iplot` backends.
    max_num_peaks : int or None, keyword only, optional
        Remove low-intensity peaks; this is the number of peaks to keep. Default is :py:const:`None` (no filtering).
        Only works with `spectrum_utils` and `spectrum_utils.iplot` backends.
    scaling : one of `{'root', 'log', 'rank'}` or None, keyword only, optional
        Scaling to apply to peak intensities. Only works with `spectrum_utils` and `spectrum_utils.iplot` backends.
    max_intensity : float or None, keyword only, optional
        Intensity of the most intense peak relative to which the peaks will be scaled
        (the default is :py:const:`None`, which means that no scaling
        relative to the most intense peak will be performed).
        Only works with `spectrum_utils` and `spectrum_utils.iplot` backends.

    Returns
    -------
    out : matplotlib.pyplot.Axes
    """
    bname = kwargs.pop('backend', 'default')
    backend = _plot_backends.get(bname)
    if backend is None:
        raise PyteomicsError('Unknown backend name: {}. Should be one of: {}.'.format(
            bname, '; '.join(_plot_backends)))

    pylab.xlabel(kwargs.pop('xlabel', 'm/z'))
    pylab.ylabel(kwargs.pop('ylabel', 'intensity'))
    if 'title' in kwargs:
        pylab.title(kwargs.pop('title'))
    return backend(spectrum, *args, **kwargs)


def _default_annotate_spectrum(spectrum, peptide, *args, **kwargs):

    # common kwargs
    types = kwargs.pop('ion_types', ('b', 'y'))
    aa_mass = kwargs.pop('aa_mass', mass.std_aa_mass)
    mass_data = kwargs.pop('mass_data', mass.nist_mass)
    ion_comp = kwargs.pop('ion_comp', mass.std_ion_comp)
    colors = {
        'a': '#388E3C',
        'b': '#1976D2',
        'c': '#00796B',
        'x': '#7B1FA2',
        'y': '#D32F2F',
        'z': '#F57C00',
    }
    colors.update(kwargs.pop('colors', {}))
    ftol = kwargs.pop('ftol', None)
    if ftol is None:
        rtol = kwargs.pop('rtol', 1e-5)
    text_kw = kwargs.pop('text_kw', dict(ha='center', clip_on=True, backgroundcolor='#ffffff99'))
    precursor_charge = kwargs.pop('precursor_charge', None)
    if precursor_charge is None:
        precursor_charge = _get_precursor_charge(spectrum)
    if precursor_charge is None:
        raise PyteomicsError('Could not extract precursor charge from spectrum. Please specify `precursor_charge` kwarg.')
    maxcharge = kwargs.pop('maxcharge', max(1, precursor_charge - 1))
    ax = kwargs.get('ax', None)
    # end of common kwargs

    # backend-specific kwargs
    centroided = kwargs.pop('centroided', True)
    adjust = kwargs.pop('adjust_text', None)
    if adjust or adjust is None:
        try:
            from adjustText import adjust_text
            adjust_kw = kwargs.pop('adjust_kw', dict(
                only_move={'text': 'y', 'points': 'y', 'objects': 'y'}, autoalign=False, force_text=(1, 1)))
        except ImportError:
            if adjust:
                raise PyteomicsError('Install adjustText for text adjustment')
            adjust = False
        else:
            if adjust is None:
                adjust = True
    # end of backend-specific kwargs

    parsed = parser.parse(peptide, True, labels=list(aa_mass) + [parser.std_cterm, parser.std_nterm])
    n = len(parsed)
    maxpeak = spectrum['intensity array'].max()
    mz, names = {}, {}
    for ion in types:
        for charge in range(1, maxcharge + 1):
            if ion[0] in 'abc':
                for i in range(2, n):
                    mz.setdefault(ion, []).append(mass.fast_mass2(parsed[:i] + [parser.std_cterm],
                        aa_mass=aa_mass, charge=charge, ion_type=ion, mass_data=mass_data, ion_comp=ion_comp))
                    names.setdefault(ion, []).append(ion[0] + str(i - 1) + ion[1:])
            else:
                for i in range(1, n - 1):
                    mz.setdefault(ion, []).append(mass.fast_mass2([parser.std_nterm] + parsed[n - (i + 1):],
                        aa_mass=aa_mass, charge=charge, ion_type=ion, mass_data=mass_data, ion_comp=ion_comp))
                    names.setdefault(ion, []).append(ion[0] + str(i) + ion[1:])
    texts = []
    for ion in types:
        c = colors.get(ion, colors.get(ion[0], 'blue'))
        matrix = np.abs(spectrum['m/z array'] - np.array(mz[ion]).reshape(-1, 1))
        if ftol is not None:
            match = np.where(matrix < ftol)
        else:
            match = np.where(matrix / spectrum['m/z array'] < rtol)
        pseudo_spec = {'m/z array': spectrum['m/z array'][match[1]], 'intensity array': spectrum['intensity array'][match[1]]}
        plot_spectrum(pseudo_spec, centroided=True, edgecolor=c, ax=ax)
        for j, i in zip(*match):
            x = spectrum['m/z array'][i]
            y = spectrum['intensity array'][i] + maxpeak * 0.02
            name = names[ion][j]
            texts.append(pylab.text(x, y, name, color=c, **text_kw))
    if adjust:
        adjust_text(texts, **adjust_kw)
    kwargs.setdefault('zorder', -1)
    return plot_spectrum(spectrum, *args, centroided=centroided, **kwargs)


def _get_precursor_charge(spectrum):
    try:
        return mgf.MGFBase.parse_precursor_charge(spectrum['params']['charge'], list_only=True)[0]
    except (PyteomicsError, KeyError):
        pass
    try:
        return int(spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['charge state'])
    except KeyError:
        pass
    return None


def _get_precursor_mz(spectrum):
    try:
        return spectrum['params']['pepmass'][0]
    except KeyError:
        pass
    try:
        return spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z']
    except KeyError:
        pass
    if 'attributes' in spectrum:
        for attr in spectrum['attributes']:
            if attr in {"MS:1000827", "MS:1000744", "MS:1002234"}:
                return spectrum['attributes'][attr]
    return None


def _spectrum_utils_create_spectrum(spectrum, *args, **kwargs):
    if sus is None:
        raise PyteomicsError('This backend requires `spectrum_utils>=0.4`.')

    # backend-specific parameters
    mz_range = kwargs.pop('mz_range', None)

    min_intensity = kwargs.pop('min_intensity', 0.0)
    max_num_peaks = kwargs.pop('max_num_peaks', None)
    scaling = kwargs.pop('scaling', None)
    max_intensity = kwargs.pop('max_intensity', None)
    spectrum = sus.MsmsSpectrum(
        'None', kwargs.pop('precursor_mz', None), kwargs.pop('precursor_charge', None),
        spectrum['m/z array'], spectrum['intensity array'])
    if mz_range:
        spectrum = spectrum.set_mz_range(*mz_range)

    spectrum = spectrum.filter_intensity(min_intensity=min_intensity, max_num_peaks=max_num_peaks
        ).scale_intensity(scaling, max_intensity)
    return spectrum


def _spectrum_utils_annotate_spectrum(spectrum, peptide, *args, **kwargs):

    # common kwargs
    aa_mass = kwargs.pop('aa_mass', mass.std_aa_mass)
    types = kwargs.pop('ion_types', ('b', 'y'))
    tol = kwargs.pop('ftol', None)
    if tol is None:
        tol = kwargs.pop('rtol', 1e-5) * 1e6
        tol_mode = 'ppm'
    else:
        tol_mode = 'Da'

    # kwargs.pop('text_kw', None)  # not used

    precursor_charge = kwargs.pop('precursor_charge', None)
    if precursor_charge is None:
        precursor_charge = _get_precursor_charge(spectrum)
    if precursor_charge is None:
        raise PyteomicsError('Could not extract precursor charge from spectrum. '
            'Please specify `precursor_charge` keyword argument.')

    maxcharge = kwargs.pop('maxcharge', max(1, precursor_charge - 1))
    # end of common kwargs

    # backend-specific parameters
    remove_precursor_peak = kwargs.pop('remove_precursor_peak', False)

    # peptide can be modX or proforma. spectrum_utils supports proforma only
    aa_comp = kwargs.get('aa_comp')
    mod_names = kwargs.get('mod_names')
    prefix = kwargs.get('prefix')

    try:
        parsed_proforma = proforma.ProForma.parse(peptide)
        peptide_pro = peptide
    except Exception:
        parsed_proforma = None
        try:
            peptide_pro = parser.to_proforma(peptide, aa_mass=aa_mass, aa_comp=aa_comp, mod_names=mod_names, prefix=prefix)
        except Exception:
            raise PyteomicsError("Cannot parse {} as ProForma or convert from modX".format(peptide))

    precursor_mz = kwargs.pop('precursor_mz', None)
    if precursor_mz is None:
        precursor_mz = _get_precursor_mz(spectrum)
    if precursor_mz is None:
        try:
            if aa_comp:
                precursor_mz = mass.calculate_mass(peptide, aa_comp=aa_comp, charge=precursor_charge)
            elif not parsed_proforma:
                precursor_mz = mass.fast_mass2(peptide, aa_mass=aa_mass, charge=precursor_charge)
            else:
                precursor_mz = mass.mass_charge_ratio(parsed_proforma.mass, precursor_charge)
        except PyteomicsError:
            raise PyteomicsError('Cannot obtain precursor m/z, please specify `precursor_mz` argument.')

    spectrum = _spectrum_utils_create_spectrum(spectrum, *args,
        precursor_mz=precursor_mz, precursor_charge=precursor_charge, **kwargs)
    if remove_precursor_peak:
        spectrum = spectrum.remove_precursor_peak(tol, tol_mode)
    spectrum = spectrum.annotate_proforma(peptide_pro, tol, tol_mode, types, maxcharge)

    return spectrum


class SpectrumUtilsColorScheme:
    """Context manager that temporarily changes `spectrum_utils.plot.colors`."""
    def __init__(self, colors):
        self.colors = colors
        self.previous_colors = sup.colors.copy()

    def __enter__(self):
        if self.colors:
            sup.colors.update(self.colors)

    def __exit__(self, *args, **kwargs):
        sup.colors = self.previous_colors


def _spectrum_utils_annotate_plot(spectrum, peptide, *args, **kwargs):

    with SpectrumUtilsColorScheme(kwargs.pop('colors', None)):
        spectrum = _spectrum_utils_annotate_spectrum(spectrum, peptide, *args, **kwargs)
        return sup.spectrum(spectrum, annot_kws=kwargs.pop('text_kw', None), ax=kwargs.pop('ax', None))


def _spectrum_utils_annotate_iplot(spectrum, peptide, *args, **kwargs):
    import spectrum_utils.iplot as supi
    with SpectrumUtilsColorScheme(kwargs.pop('colors', None)):
        spectrum = _spectrum_utils_annotate_spectrum(spectrum, peptide, *args, **kwargs)
        return supi.spectrum(spectrum, annot_kws=kwargs.pop('text_kw', None))


_annotation_backends = {
    'default': _default_annotate_spectrum,
    'spectrum_utils': _spectrum_utils_annotate_plot,
    'spectrum_utils.iplot': _spectrum_utils_annotate_iplot,
}


def annotate_spectrum(spectrum, peptide, *args, **kwargs):
    """Plot a spectrum and annotate matching fragment peaks.

    Parameters
    ----------
    spectrum : dict
        A spectrum as returned by Pyteomics parsers. Needs to have 'm/z array' and 'intensity array' keys.
    peptide : str
        A modX sequence.
    backend : str, keyword only, optional
        One of `{'default', 'spectrum_utils', 'spectrum_utils.iplot'}`.
        The `spectrum_utils` backend requires installing :py:mod:`spectrum_utils`.
        The `spectrum_utils.iplot` backend requires installing :py:mod:`spectrum_utils[iplot]`.
    ion_types : Container, keyword only, optional
        Ion types to be considered for annotation. Default is `('b', 'y')`.
    precursor_charge : int, keyword only, optional
        If not specified, an attempt is made to extract it from `spectrum`.
    maxcharge : int, keyword only, optional
        Maximum charge state for fragment ions to be considered. Default is `precursor_charge - 1`.
    colors : dict, keyword only, optional
        Keys are ion types, values are colors to plot the annotated peaks with. Default depends on backend.
    ftol : float, keyword only, optional
        A fixed m/z tolerance value for peak matching. Alternative to `rtol`.
    rtol : float, keyword only, optional
        A relative m/z error for peak matching. Default is 10 ppm.
    aa_mass : dict, keyword only, optional
        A dictionary of amino acid residue masses.
    text_kw : dict, keyword only, optional
        Keyword arguments for :py:func:`pylab.text`.
    xlabel : str, keyword only, optional
        Label for the X axis. Default is "m/z". Does not work with `spectrum_utils.iplot` backend.
    ylabel : str, keyword only, optional
        Label for the Y axis. Default is "intensity". Does not work with `spectrum_utils.iplot` backend.
    title : str, keyword only, optional
        The title. Empty by default. Does not work with `spectrum_utils.iplot` backend.
    ax : matplotlib.pyplot.Axes, keyword only, optional
        Axes to draw the spectrum. Does not work with `spectrum_utils.iplot` backend.

    *args
        Passed to the plotting backend.
    **kwargs
        Passed to the plotting backend.

    centroided : bool, keyword only, optional
        Passed to :py:func:`plot_spectrum`. Only works with `default` backend.
    ion_comp : dict, keyword only, optional
        A dictionary defining ion compositions to override :py:const:`pyteomics.mass.std_ion_comp`.
        Only works with `default` backend.
    mass_data : dict, keyword only, optional
        A dictionary of element masses to override :py:const:`pyteomics.mass.nist_mass`.
        Only works with `default` backend.

    adjust_text : bool, keyword only, optional
        Adjust the overlapping text annotations using :py:mod:`adjustText`. Only works with `default` backend.
    adjust_kw : dict, keyword only, optional
        Keyword arguments for :py:func:`adjust_text`. Only works with `default` backend.

    remove_precursor_peak : bool, keyword only, optional
        Remove precursor peak from spectrum before annotation. Default is :py:const:`False`.
        Only works with `spectrum_utils` backend.
    min_intensity : float, keyword only, optional
        Remove low-intensity peaks; this is a factor of maximum peak intensity. Default is 0 (no filtering).
        Only works with `spectrum_utils` and `spectrum_utils.iplot` backends.
    max_num_peaks : int or None, keyword only, optional
        Remove low-intensity peaks; this is the number of peaks to keep. Default is :py:const:`None` (no filtering).
        Only works with `spectrum_utils` and `spectrum_utils.iplot` backends.
    scaling : one of `{'root', 'log', 'rank'}` or None, keyword only, optional
        Scaling to apply to peak intensities. Only works with `spectrum_utils` and `spectrum_utils.iplot` backends.
    max_intensity : float or None, keyword only, optional
        Intensity of the most intense peak relative to which the peaks will be scaled
        (the default is :py:const:`None`, which means that no scaling
        relative to the most intense peak will be performed).
        Only works with `spectrum_utils` and `spectrum_utils.iplot` backends.
    aa_comp : dict, keyword only, optional
        Amino acid compositions, including modified ones. If given, will be used for conversion from *modX* to ProForma.
    mod_names : dict or callable, keyword only, optional
        If given, will be used for conversion from *modX* to ProForma.
    prefix : str, keyword only, optional
        If given, will be used for conversion from *modX* to ProForma.

    Returns
    -------
    out : matplotlib.pyplot.Axes
    """
    bname = kwargs.pop('backend', 'default')
    backend = _annotation_backends.get(bname)
    if backend is None:
        raise PyteomicsError('Unknown backend name: {}. Should be one of: {}.'.format(
            bname, '; '.join(_annotation_backends)))

    pylab.xlabel(kwargs.pop('xlabel', 'm/z'))
    pylab.ylabel(kwargs.pop('ylabel', 'intensity'))
    pylab.title(kwargs.pop('title', ''))
    return backend(spectrum, peptide, *args, **kwargs)


def _spectrum_utils_mirror(spec_top, spec_bottom, spectrum_kws=None, ax=None, **kwargs):
    with SpectrumUtilsColorScheme(kwargs.pop('colors', None)):
        ax = sup.mirror(spec_top, spec_bottom, spectrum_kws=spectrum_kws, ax=ax)
        ax.set_xlabel(kwargs.pop('xlabel', 'm/z'))
        ax.set_ylabel(kwargs.pop('ylabel', 'intensity'))
        ax.set_title(kwargs.pop('title', ''))
        return ax


def _spectrum_utils_iplot_mirror(spec_top, spec_bottom, spectrum_kws=None, **kwargs):
    import spectrum_utils.iplot as supi
    with SpectrumUtilsColorScheme(kwargs.pop('colors', None)):
        return supi.mirror(spec_top, spec_bottom, spectrum_kws=spectrum_kws)


_mirror_backends = {
    'spectrum_utils': _spectrum_utils_mirror,
    'spectrum_utils.iplot': _spectrum_utils_iplot_mirror,
}


def mirror(spec_top, spec_bottom, peptide=None, spectrum_kws=None, ax=None, **kwargs):
    """Create a mirror plot of two (possible annotated) spectra using `spectrum_utils`.

    Parameters
    ----------
    spec_top : dict
        A spectrum as returned by Pyteomics parsers. Needs to have 'm/z array' and 'intensity array' keys.
    spec_bottom : dict
        A spectrum as returned by Pyteomics parsers. Needs to have 'm/z array' and 'intensity array' keys.
    peptide : str or None, optional
        A modX sequence or ProForma. If provided, the peaks will be annotated as peptide fragments.
    spectrum_kws : dict or None, optional
        Passed to :py:func:`spectrum_utils.plot.mirror`.
    backend : str, keyword only, optional
        One of {'spectrum_utils', 'spectrum_utils.iplot'}. Default is 'spectrum_utils'.

        .. note ::
            Requires :py:mod:`spectrum_utils` or :py:mod:`spectrun_utils[iplot]`, respectively.

    ax : matplotlib.pyplot.Axes or None, optional
        Passed to :py:func:`spectrum_utils.plot.mirror`. Works only for the 'spectrum_utils' backend.
    xlabel : str, keyword only, optional
        Label for the X axis. Default is "m/z". Works only for the 'spectrum_utils' backend.
    ylabel : str, keyword only, optional
        Label for the Y axis. Default is "intensity". Works only for the 'spectrum_utils' backend.
    title : str, keyword only, optional
        The title. Empty by default. Works only for the 'spectrum_utils' backend.

    **kwargs : same as for :py:func:`annotate_spectrum` for `spectrum_utils` backends.

    Returns
    -------
    out : matplotlib.pyplot.Axes
    """

    spec_gen = _spectrum_utils_create_spectrum if peptide is None else _spectrum_utils_annotate_spectrum
    spec_top = spec_gen(spec_top, peptide, **kwargs)
    spec_bottom = spec_gen(spec_bottom, peptide, **kwargs)

    bname = kwargs.pop('backend', 'spectrum_utils')
    backend = _mirror_backends.get(bname)
    if backend is None:
        raise PyteomicsError('Unknown backend name: {}. Should be one of: {}.'.format(
            bname, '; '.join(_mirror_backends)))
    backend_kw = {'spectrum_kws': spectrum_kws}
    if bname == 'spectrum_utils':
        backend_kw['ax'] = ax
    backend_kw.update(kwargs)
    return backend(spec_top, spec_bottom, **backend_kw)
