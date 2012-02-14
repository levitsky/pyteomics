"""
auxiliary - common functions and objects 
========================================

Math
----

  :py:func:`linear_regression` - a wrapper for numpy linear regression

Project infrastructure
----------------------

  :py:class:`PyteomicsError` - a pyteomics-specific exception

-------------------------------------------------------------------------------

"""
# Licensed under the MIT license:
# http://www.opensource.org/licenses/mit-license.php 
import numpy

class PyteomicsError(Exception):
    """Exception raised for errors in Pyteomics library.

    Attributes
    ----------
    msg : str
        Error message.
    """

    def __init__(self, msg):
        self.msg = msg
        
    def __str__(self):
        return "Pyteomics error, message: %s" % (repr(self.msg),)

def linear_regression(x, y, a=None, b=None):
    """Calculate coefficients of a linear regression y = a * x + b.

    Parameters
    ----------
    x, y : array_like of float
    a : float, optional
        If specified then the slope coefficient is fixed and equals a.
    b : float, optional        
        If specified then the free term is fixed and equals b.
    
    Returns
    -------
    out : 4-tuple of float
        The structure is (a, b, r, stderr), where
        a -- slope coefficient,
        b -- free term,
        r -- Peason correlation coefficient,
        stderr -- standard deviation.
    """

    if (a!=None and b==None):
        b = numpy.mean([y[i] - a * x[i] for i in range(len(x))])
    elif (a!=None and b!= None):
        pass
    else:
        a, b = numpy.polyfit(x, y, 1)

    r = numpy.corrcoef(x, y)[0, 1]
    stderr = numpy.std([y[i] - a * x[i] - b for i in range(len(x))])

    return (a, b, r, stderr)
