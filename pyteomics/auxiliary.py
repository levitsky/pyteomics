class PyteomicsError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

def linear_regression(x, y, a=None, b=None):
    """Calculate coefficients of a linear regression y = a * x + b.

    Parameters
    ----------
    x, y : list or 1-D numpy array of float
    a : float
        If specified then the slope coefficient is fixed and equals a.
    b : float        
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
