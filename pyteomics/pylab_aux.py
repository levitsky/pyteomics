"""
pylab_aux - auxiliary functions for plotting with pylab
=======================================================

This module serves as a collection of useful routines for data plotting with
matplotlib.

Data plotting
-------------

  :py:func:`plot_line` - plot a line.

  :py:func:`scatter_trend` - plot a scatter plot with a regression line.

  :py:func:`plot_function_3d` - plot a 3D graph of a function of two variables.

  :py:func:`plot_function_contour` - plot a contour graph of a function of
  two variables.

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

def plot_line(a, b, **kwargs):
    """Plot a line y = a * x + b.

    Parameters
    ----------
    a, b : float
        The slope and intercept of the line.

    **kwargs : passed to :py:func:`pylab.plot`.
    """
    xlim = pylab.xlim()
    ylim = pylab.ylim()
    pylab.plot([xlim[0], xlim[1]],
               [a * xlim[0] + b, a * xlim[1] + b],
               **kwargs)

def scatter_trend(x, y, **kwargs):
    """Make a scatter plot with a linear regression.

    Parameters
    ----------
    x, y : array_like of float
    plot_trend : bool, optional
        If :py:const:`True` then plot a trendline (default).
    plot_sigmas : bool, optional
        If :py:const:`True` then plot confidence intervals of the linear fit.
        :py:const:`False` by default.
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
        Keyword arguments for :py:func:`pylab.plot`.
        Empty by default.
    legend_kwargs : dict, optional
        Keyword arguments for :py:func:`pylab.legend`.
        Default is :py:const:`{'loc': 'upper left'}`.
    """
    a, b, r, stderr = linear_regression(x, y)
    pylab.title(kwargs.get('title', ''))
    pylab.xlabel(kwargs.get('xlabel', ''))
    pylab.ylabel(kwargs.get('ylabel', ''))
    scat_plot = pylab.scatter(x, y, **kwargs.get('scatter_kwargs', {}))
    scat_plot.set_label(
        '$y\,=\,{:.3f}x\,{}\,{:.3f}$, '
        '$R^2=\,{:.3f}$ \n$\sigma\,=\,{:.3f}$'.format(
            a, '-' if b < 0 else '+', abs(b), r*r, stderr))
    legend = pylab.legend(**kwargs.get('legend_kwargs', {'loc': 'upper left'}))
    legend_frame = legend.get_frame()
    legend_frame.set_alpha(kwargs.get('alpha_legend', 1.0))
    if kwargs.get('plot_trend', True):
        pylab.plot([min(x), max(x)],
                   [a*min(x)+b, a*max(x)+b],
                   **kwargs.get('plot_kwargs', {}))
    if kwargs.get('plot_sigmas', False):
        for i in [-3.0,-2.0,-1.0,1.0,2.0,3.0]:
            pylab.plot([min(x), max(x)],
                       [a*min(x)+b+i*stderr, a*max(x)+b+i*stderr],
                       'r--')

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

def plot_qvalue_curve(qvalues, **kwargs):
    """
    Plot a curve with q-values on the X axis and corresponding PSM number
    (starting with ``1``) on the Y axis.

    Parameters
    ----------
    qvalues : structured NumPy array
        An array of q-values and scores, as returned by
        :py:func:`pyteomics.tandem.qvalues`,
        :py:func:`pyteomics.pepxml.qvalues`,
        :py:func:`pyteomics.mzid.qvalues`,
        :py:func:`pyteomics.auxiliary.qvalues` or your own function.
    xlabel : str, optional
        Label for the X axis. Default is "q-value".
    ylabel : str, optional
        Label for the Y axis. Default is "# of PSMs".
    title : str, optional
        The title. Empty by default.
    **kwargs : will be given to :py:func:`pylab.plot`.
    """
    pylab.xlabel(kwargs.pop('xlabel', 'q-value'))
    pylab.ylabel(kwargs.pop('ylabel', '# of PSMs'))
    pylab.title(kwargs.pop('title', ''))
    pylab.plot(qvalues, 1+np.arange(qvalues.size), **kwargs)
