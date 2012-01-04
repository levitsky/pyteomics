import pylab
import numpy

def linear_regression(x, y, a=None, b=None):
    if (a!=None and b==None):
        b = numpy.mean([y[i] - a * x[i] for i in range(len(x))])
    elif (a!=None and b!= None):
        pass
    else:
        a, b = numpy.polyfit(x, y, 1)

    r = numpy.corrcoef(x, y)[0, 1]
    stderr = numpy.std([y[i] - a * x[i] - b for i in range(len(x))])

    return (a, b, r, stderr)

def plot_line(a, b, **kwargs):
    xlim = pylab.xlim()
    ylim = pylab.ylim()
    pylab.plot([xlim[0], xlim[1]], 
               [a * xlim[0] + b, a * xlim[1] + b],
               **kwargs) 

def scatter_trend(x, y, **kwargs):
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

def plot_3d(x, y, function, **kwargs):
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
    elif plot_type == 'contourf3D':
        ax.contourf3D(X, Y, Z, 20, cmap=pylab.cm.jet)
    ax.set_xlabel(kwargs.get('xlabel', ''))
    ax.set_ylabel(kwargs.get('ylabel', ''))
    ax.set_zlabel(kwargs.get('zlabel', ''))
    ax.set_title(kwargs.get('title', ''))

def contourf(x, y, function, **kwargs):
    X, Y = numpy.meshgrid(x, y)
    Z = []
    for y_value in y:
        Z.append([])
        for x_value in x:
            Z[-1].append(function(x_value, y_value))
    Z = numpy.array(Z)
    pylab.contourf(X, Y, Z, 50, cmap=pylab.cm.jet)
    pylab.xlabel(kwargs.get('xlabel', ''))
    pylab.ylabel(kwargs.get('ylabel', ''))
    pylab.title(kwargs.get('title', ''))

def scatter_3d(x, y, z, **kwargs):
    import mpl_toolkits.mplot3d.axes3d as pylab3d
    ax = pylab3d.Axes3D(pylab.gcf())
    X = numpy.array(x)
    Y = numpy.array(y)
    Z = numpy.array(z)
    ax.scatter3D(X, Y, Z, cmap=pylab.cm.jet)
    ax.set_xlabel(kwargs.get('xlabel', ''))
    ax.set_ylabel(kwargs.get('ylabel', ''))
    ax.set_zlabel(kwargs.get('zlabel', ''))
    ax.set_title(kwargs.get('title', ''))
