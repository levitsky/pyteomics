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

-------------------------------------------------------------------------------

"""
# Licensed under the MIT license:
# http://www.opensource.org/licenses/mit-license.php 

import pylab
import numpy
from auxiliary import linear_regression, PyteomicsError

def plot_line(a, b, **kwargs):
    """Plot a line y = a * x + b.
    
    Parameters
    ----------
    a, b : float
        The slope and intercept of the line.
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
        If True then plot a trendline. True by default.
    plot_sigmas : bool, optional
        If True then plot confidence intervals of the linear fit.
        False by default.
    title : str, optional
        The title. Empty by default.
    xlabel, ylabel : str, optional
        The axes labels. Empty by default.
    alpha : float, optional
        Transparency of points. 1.0 by default
    alpha_legend : float, optional
        Legend box transparency. 1.0 by default
    """
    a, b, r, stderr = linear_regression(x, y)
    pylab.title(kwargs.get('title', ''))
    pylab.xlabel(kwargs.get('xlabel', ''))
    pylab.ylabel(kwargs.get('ylabel', ''))
    scat_plot = pylab.scatter(x, y,
                              c=kwargs.get('c', 'b'),
                              alpha=kwargs.get('alpha', 1.0))
    scat_plot.set_label(
        '$y\,=\,%.3fx\,+\,%.3f$, $R^2=\,%.3f$ \n$\sigma\,=\,%.3f$' % (
            a, b, r*r, stderr))
    legend = pylab.legend(loc='upper left')
    legend_frame = legend.get_frame()
    legend_frame.set_alpha(kwargs.get('alpha_legend', 1.0))
    if kwargs.get('plot_trend', True):
        pylab.plot([min(x), max(x)],
                   [a*min(x)+b, a*max(x)+b])
    if kwargs.get('plot_sigmas', False):
        for i in [-3.0,-2.0,-1.0,1.0,2.0,3.0]:
            pylab.plot([min(x), max(x)],
                       [a*min(x)+b+i*stderr, a*max(x)+b+i*stderr],
                       'r--')

def plot_function_3d(x, y, function, **kwargs):
    """Plot values of a function of two variables in 3D.

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

    See also
    --------
    More on 3D plotting in pylab:
    http://www.scipy.org/Cookbook/Matplotlib/mplot3D

    """
    import mpl_toolkits.mplot3d.axes3d as pylab3d
    ax = pylab3d.Axes3D(pylab.gcf())
    X, Y = numpy.meshgrid(x, y)
    Z = []
    for y_value in y:
        Z.append([])
        for x_value in x:
            Z[-1].append(function(x_value, y_value))
    Z = numpy.array(Z)
    plot_type = kwargs.get('plot_type', 'surface')
    if plot_type == 'surface':
        ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=pylab.cm.jet)
    elif plot_type == 'wireframe':
        ax.plot_wireframe(X, Y, Z, cmap=pylab.cm.jet)
    elif plot_type == 'scatter':
        ax.scatter3D(numpy.ravel(X), numpy.ravel(Y), numpy.ravel(Z))
    elif plot_type == 'contour':
        num_contours = kwargs.get('num_contours', 50)
        ax.contour3D(X, Y, Z, num_contours, cmap=pylab.cm.jet)
    elif plot_type == 'contourf':
        num_contours = kwargs.get('num_contours', 50)
        ax.contourf3D(X, Y, Z, num_contours, cmap=pylab.cm.jet)
    else:
        raise PyteomicsError('Unknown plot type: %s' % (plot_type,))
    ax.set_xlabel(kwargs.get('xlabel', ''))
    ax.set_ylabel(kwargs.get('ylabel', ''))
    ax.set_zlabel(kwargs.get('zlabel', ''))
    ax.set_title(kwargs.get('title', ''))

def plot_function_contour(x, y, function, **kwargs):
    """Plot values of a function of two variables in 3D.

    Parameters
    ----------
    x, y : array_like of float
        The plotting range.
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

    """
    X, Y = numpy.meshgrid(x, y)
    Z = []
    for y_value in y:
        Z.append([])
        for x_value in x:
            Z[-1].append(function(x_value, y_value))
    Z = numpy.array(Z)
    num_contours = kwargs.get('num_contours', 50)
    if kwargs.get('filling', True):
        pylab.contourf(X, Y, Z, num_contours, cmap=pylab.cm.jet)
    else:
        pylab.contour(X, Y, Z, num_contours, cmap=pylab.cm.jet)
    pylab.xlabel(kwargs.get('xlabel', ''))
    pylab.ylabel(kwargs.get('ylabel', ''))
    pylab.title(kwargs.get('title', ''))

