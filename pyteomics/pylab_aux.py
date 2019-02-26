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

This module requires :py:mod:`matplotlib`.

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
from . import parser, mass

def plot_line(a, b, xlim=None, *args, **kwargs):
    """Plot a line y = a * x + b.

    Parameters
    ----------
    a, b : float
        The slope and intercept of the line.
    xlim : tuple, optional
        Minimal and maximal values of `x`. If not given, :py:func:`pylab.xlim`
        will be called.

    *args, **kwargs : passed to :py:func:`pylab.plot` after `x` and `y` values.

    Returns
    -------
    out : matplotlib.lines.Line2D
        The line object.
    """
    if xlim is None: xlim = pylab.xlim()
    return pylab.plot([xlim[0], xlim[1]],
               [a * xlim[0] + b, a * xlim[1] + b],
               *args, **kwargs)

def scatter_trend(x, y=None, **kwargs):
    """Make a scatter plot with a linear regression.

    Parameters
    ----------
    x, y : array_like of float
        1-D arrays of floats. If `y` is omitted, `x` must be a 2-D array of shape (N, 2).
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
            s_lines.append(plot_line(a, b + i*stderr, xlim, **sigma_kwargs))
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
    x, y : array_like of float
        The plotting range.
    function : function
        The function to plot.
    plot_type : {'surface', 'wireframe', 'scatter', 'contour', 'contourf'}
        The type of a plot, see
        `scipy cookbook <http://www.scipy.org/Cookbook/Matplotlib/mplot3D>`_
        for examples. The default value is 'surface'.
    num_contours : int
        The number of contours to plot, 50 by default.
    xlabel, ylabel, zlabel : str, optional
        The axes labels. Empty by default.
    title : str, optional
        The title. Empty by default.
    **kwargs : passed to the respective plotting function.
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
    **kwargs : passed to :py:func:`pylab.contour` or :py:func:`pylab.contourf`.
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
    xlabel : str, optional
        Label for the X axis. Default is "q-value".
    ylabel : str, optional
        Label for the Y axis. Default is "# of PSMs".
    title : str, optional
        The title. Empty by default.
    *args, **kwargs : will be given to :py:func:`pylab.plot` after `x` and `y`.

    Returns
    -------
    out : matplotlib.lines.Line2D
    """
    pylab.xlabel(kwargs.pop('xlabel', 'q-value'))
    pylab.ylabel(kwargs.pop('ylabel', '# of PSMs'))
    pylab.title(kwargs.pop('title', ''))
    return pylab.plot(qvalues, 1+np.arange(qvalues.size), *args, **kwargs)

def plot_spectrum(spectrum, centroided=True, *args, **kwargs):
    """
    Plot a spectrum, assuming it is a dictionary containing "m/z array" and "intensity array".

    Parameters
    ----------
    spectrum : dict
        A dictionary, as returned by MGF, mzML or mzXML parsers.
        Must contain "m/z array" and "intensity array" keys with decoded arrays.
    centroided : bool, optional
        If :py:const:`True` (default), peaks of the spectrum are plotted using :py:func:`pylab.bar`.
        If :py:const:`False`, the arrays are simply plotted using :py:func:`pylab.plot`.
    xlabel : str, optional
        Label for the X axis. Default is "m/z".
    ylabel : str, optional
        Label for the Y axis. Default is "intensity".
    title : str, optional
        The title. Empty by default.
    *args, **kwargs : will be given to :py:func:`pylab.plot` or :py:func:`pylab.bar` (depending on `centroided`).
    """
    pylab.xlabel(kwargs.pop('xlabel', 'm/z'))
    pylab.ylabel(kwargs.pop('ylabel', 'intensity'))
    pylab.title(kwargs.pop('title', ''))
    if centroided:
        kwargs.setdefault('align', 'center')
        kwargs.setdefault('width', 0)
        kwargs.setdefault('linewidth', 1)
        kwargs.setdefault('edgecolor', 'k')
        return pylab.bar(spectrum['m/z array'], spectrum['intensity array'], *args, **kwargs)
    return pylab.plot(spectrum['m/z array'], spectrum['intensity array'], *args, **kwargs)


def annotate_spectrum(spectrum, peptide, centroided=True, *args, **kwargs):
    """Plot a spectrum and annotate matching fragment peaks.

    Parameters
    ----------
    spectrum : dict
        A spectrum as returned by Pyteomics parsers. Needs to have 'm/z array' and 'intensity array' keys.
    peptide : str
        A modX sequence.
    centroided : bool, optional
        Passed to :py:func:`plot_spectrum`.
    types : Container, optional
        Ion types to be considered for annotation. Default is `('b', 'y')`.
    colors : dict, optional
        Keys are ion types, values are colors to plot the annotated peaks with. Defaults to a red-blue scheme.
    ftol : float, optional
        A fixed m/z tolerance value for peak matching. Alternative to `rtol`.
    rtol : float, optional
        A relative m/z error for peak matching. Default is 10 ppm.
    adjust_text : bool, optional
        Adjust the overlapping text annotations using :py:mod:`adjustText`.
    text_kw : dict, optional
        Keyword arguments for :py:func:`pylab.text`.
    adjust_kw : dict, optional
        Keyword argyuments for `:py:func:`adjust_text`.
    ion_comp : dict, optional
        A dictionary defining definitions of ion compositions to override :py:const:`pyteomics.mass.std_ion_comp`.
    mass_data : dict, optional
        A dictionary of element masses to override :py:const:`pyteomics.mass.nist_mass`.
    aa_mass : dict, optional
        A dictionary of amino acid residue masses.
    *args, **kwargs : passed to :py:func:`plot_spectrum`.
    """
    types = kwargs.pop('types', ('b', 'y'))
    maxcharge = kwargs.pop('maxcharge', 1)
    aa_mass = kwargs.pop('aa_mass', mass.std_aa_mass)
    mass_data = kwargs.pop('mass_data', mass.nist_mass)
    ion_comp = kwargs.pop('ion_comp', mass.std_ion_comp)
    std_colors = {i: 'red' for i in 'xyz'}
    std_colors.update({i: 'blue' for i in 'abc'})
    colors = kwargs.pop('colors', std_colors)
    ftol = kwargs.pop('ftol', None)
    if ftol is None:
        rtol = kwargs.pop('rtol', 1e-5)
    text_kw = kwargs.pop('text_kw', dict(ha='center', clip_on=True, backgroundcolor='#ffffff99'))
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
    parsed = parser.parse(peptide, True, labels=list(aa_mass) + [parser.std_cterm, parser.std_nterm])
    n = len(parsed)
    maxpeak = spectrum['intensity array'].max()
    mz, names = {}, {}
    for ion in types:
        for charge in range(1, maxcharge+1):
            if ion[0] in 'abc':
                for i in range(2, n):
                    mz.setdefault(ion, []).append(mass.fast_mass2(parsed[:i] + [parser.std_cterm],
                        aa_mass=aa_mass, charge=charge, ion_type=ion, mass_data=mass_data, ion_comp=ion_comp))
                    names.setdefault(ion, []).append(ion[0] + str(i-1) + ion[1:])
            else:
                for i in range(1, n-1):
                    mz.setdefault(ion, []).append(mass.fast_mass2([parser.std_nterm] + parsed[n-(i+1):],
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
        plot_spectrum(pseudo_spec, centroided=True, edgecolor=c)
        for j, i in zip(*match):
            x = spectrum['m/z array'][i]
            y = spectrum['intensity array'][i] + maxpeak * 0.02
            name = names[ion][j]
            texts.append(pylab.text(x, y, name, color=c, **text_kw))
    if adjust:
        adjust_text(texts, **adjust_kw)
    kwargs.setdefault('zorder', -1)
    plot_spectrum(spectrum, centroided, *args, **kwargs)
