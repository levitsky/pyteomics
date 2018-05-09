from .structures import PyteomicsError


def linear_regression_vertical(x, y=None, a=None, b=None):
    """Calculate coefficients of a linear regression y = a * x + b.
    The fit minimizes *vertical* distances between the points and the line.

    Requires :py:mod:`numpy`.

    Parameters
    ----------
    x, y : array_like of float
        1-D arrays of floats. If `y` is omitted, `x` must be a 2-D array of shape (N, 2).
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

    import numpy as np
    x = np.array(x, copy=False)
    if y is not None:
        y = np.array(y, copy=False)
    else:
        if len(x.shape) != 2 or x.shape[-1] != 2:
            raise PyteomicsError(
                'If `y` is not given, x.shape should be (N, 2), given: {}'.format(x.shape))
        y = x[:, 1]
        x = x[:, 0]
    if (a is not None and b is None):
        b = (y - a * x).mean()
    elif (a is not None and b is not None):
        pass
    else:
        a, b = np.polyfit(x, y, 1)

    r = np.corrcoef(x, y)[0, 1]
    stderr = (y - a * x - b).std()

    return a, b, r, stderr


def linear_regression(x, y=None, a=None, b=None):
    """Alias of :py:func:`linear_regression_vertical`."""
    return linear_regression_vertical(x, y, a, b)


def linear_regression_perpendicular(x, y=None):
    """Calculate coefficients of a linear regression y = a * x + b.
    The fit minimizes *perpendicular* distances between the points and the line.

    Requires :py:mod:`numpy`.

    Parameters
    ----------
    x, y : array_like of float
        1-D arrays of floats. If `y` is omitted, `x` must be a 2-D array of shape (N, 2).

    Returns
    -------
    out : 4-tuple of float
        The structure is (a, b, r, stderr), where
        a -- slope coefficient,
        b -- free term,
        r -- Peason correlation coefficient,
        stderr -- standard deviation.
    """

    import numpy as np
    x = np.array(x, copy=False)
    if y is not None:
        y = np.array(y, copy=False)
        data = np.hstack((x.reshape((-1, 1)), y.reshape((-1, 1))))
    else:
        if len(x.shape) != 2 or x.shape[-1] != 2:
            raise PyteomicsError(
                'If `y` is not given, x.shape should be (N, 2), given: {}'.format(x.shape))
        data = x
    mu = data.mean(axis=0)
    eigenvectors, eigenvalues, V = np.linalg.svd((data - mu).T, full_matrices=False)
    a = eigenvectors[0][1] / eigenvectors[0][0]
    xm, ym = data.mean(axis=0)
    b = ym - a * xm

    r = np.corrcoef(data[:, 0], data[:, 1])[0, 1]
    stderr = ((data[:, 1] - a * data[:, 0] - b) / np.sqrt(a**2 + 1)).std()

    return a, b, r, stderr
