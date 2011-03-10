#ifndef AUXILIARY_H
#define AUXILIARY_H

namespace BioLCCC
{

//! Calculates the second derivatives of a function.
/*!
    The function should be described by its \a n values \a y at consecutive
    points \a x. The function stores its output in the array \a y2 of size \a n.
 */
void fitSpline(const double *x, const double *y, const int n, double * y2);

//! Calculates the value of a function using the cubic spline interpolation.
/*!
    The value of the function at arbitrary point \a x_in is calculated using 
    the values of the function \a y and its second derivative \a y2 
    at \a n consecutive points \a x.
 */
double calculateSpline(const double *x, const double *y, const double * y2,
    const int n, const double x_in);

//! Calculates the value of a function using a piecewise linear interpolation.
/*!
    Calculates the value of a function at the point \a x_in. The values of the
    function at \a n points \a x are known and equal \a y.
*/
double linInterpolate(const double * x, const double * y, const int n, 
                      const double x_in);

//! Calculates the value of a function using a polynomial interpolation.
/*!
    Finds the value of a function at the point \a x_in. The values of the
    function at \a n points \a x are known and equal \a y.
 */
double polInterpolate(const double * x, const double * y, const int n, 
                      const double x_in);

//! Calculates the value of a function using a partial polynomial interpolation.
/*!
    This version calculates the value of the function at the point \a x_in
    using only the values of the function at \a n_part * 2 nearest points.
 */
double partPolInterpolate(const double * x, const double * y, 
    const int n, const int n_part, const double x_in);

//! Solves a linear matrix equation m * x = rhs for square matrix m of size nxn
/*!
    The equation is solved using the Gauss-Jordan elimination method. The
    resulting solution vector x is written to \a rhs.
 */
void solveMatrixEquation(double * m, double * rhs, const int n);

//! Constructs a polynomial whose values at n points x equals y.
/*!
    The computed coefficients before terms x^0, x^1, ... , x^(n-1) are 
    stored in \a y.
 */
void fitPolynomial(double * x, double * y, const int n);

//! Calculates the value of a polynomial of (n-1)th power at point x.
/*!
    The computed coefficients before terms x^0, x^1, ... , x^(n-1) should
    stored in \a coeffs.
 */
double calculatePolynomial(const double * coeffs, const int n, const double x);

}

#endif
