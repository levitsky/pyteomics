#include <iostream>
#include <cmath>
#include "auxiliary.h"

namespace BioLCCC
{

void fitSpline(const double *x, const double *y, const int n, double * y2)
{
    double a,b,c,d;
    double * c1 = new double[n-1];
    double * d1 = new double[n-1];
    c1[0] = 0.0;
    d1[0] = 0.0;
    
    for (int i = 1; i < n - 1; i++)
    {
        a = x[i] - x[i-1];
        b = 2.0 * (x[i+1] - x[i-1]);
        c = x[i+1] - x[i];
        d = 6.0 * ((y[i+1] - y[i]) / (x[i+1] - x[i]) 
                    - (y[i] - y[i-1]) / (x[i] - x[i-1]));
        c1[i] = c / (b - c1[i-1] * a);
        d1[i] = (d - d1[i-1] * a) / (b - c1[i-1] * a);
    }

    y2[n-1] = 0.0;
    for (int i = n - 2; i >= 0; i--)
    {
        y2[i] = d1[i] - c1[i] * y2[i+1];
    }

    delete[] c1, d1;
}

double calculateSpline(const double *x, const double *y, const double * y2,
    const int n, const double x_in)
{
    int j = 0;
    int j_up = n - 1;
    while (j_up > j + 1) 
    {
        if ((x[j] <= x_in) && (x_in <= x[(j + j_up) / 2]))
        {
            j_up = (j + j_up) / 2;
        }
        else
        {
            j = (j + j_up) / 2;
        }
    }

    double dx = x[j+1] - x[j];
    double a = (x[j+1] - x_in) / dx;
    double b = (x_in - x[j]) / dx;
    return a * y[j] + b * y[j+1] 
        + ((a * a * a - a) * y2[j] + (b * b * b - b) * y2[j + 1])
          * (dx * dx) / 6.0;
}

double linInterpolate(const double * x, const double * y, const int n,
                      const double x_in)
{
    for (int i=0; i<n-1; ++i)
    {
        if ((x[i] <= x_in) && (x_in <= x[i+1]))
        {
            return y[i] + (y[i+1] - y[i]) * (x_in - x[i]) / (x[i+1] - x[i]);
        }
    }
    return y[n];
}

double polInterpolate(const double * x, const double * y, const int n, 
                      const double x_in)
{
    double * p = new double[n];
    for ( int i=0; i<n; i++)
    {
        p[i] = y[i];
    }
    for (int i=1; i<n; i++)
    {
        for (int j=0; j<n-i; j++)
        {
            p[j] = (p[j] * (x_in - x[j+i]) + p[j+1] * (x[j] - x_in))
                   / (x[j] - x[j+i]);
        }
    }
    double output = p[0];
    delete[] p;
    return output;
}

double partPolInterpolate(const double * x, const double * y, 
    const int n, const int n_part, const double x_in)
{
    int k = 0;
    int k_up = n - 1;
    while (k_up > k + 1) 
    {
        if ((x[k] <= x_in) && (x_in <= x[(k + k_up) / 2]))
        {
            k_up = (k + k_up) / 2;
        }
        else
        {
            k = (k + k_up) / 2;
        }
    }

    k = ((k - n_part + 1) > 0) ? (k - n_part + 1) : 0;
    k = (k + 2 * n_part < n ) ? k : n - 2 * n_part;

    double * p = new double[n_part*2];
    for (int i=0; i<n_part*2; i++)
    {
        p[i] = y[k+i];
    }
    for (int i=1; i<n_part*2; i++)
    {
        for (int j=0; j<n_part*2-i; j++)
        {
            p[j] = (p[j] * (x_in - x[k+j+i]) + p[j+1] * (x[k+j] - x_in))
                   / (x[k+j] - x[k+j+i]);
        }
    }
    double output = p[0];
    delete[] p;
    return output;
}

void solveMatrixEquation(double * m, double * rhs, const int n)
{
    double temp;
    bool * reduced = new bool[n];
    for (int i=0; i<n; i++)
    {
        reduced[i] = false;
    }

    for (int k=0; k<n; k++) 
    {
        int pivot_i=0;
        int pivot_j=0;
        temp = 0.0;

        for (int i=0; i<n; i++)
        {
            if (!reduced[i])
            {
                for (int j=0; j<n; j++)
                {
                    if ((!reduced[j]) && (fabs(m[i * n + j]) >= temp))
                    {
                        pivot_i = i;
                        pivot_j = j;
                        temp = fabs(m[i * n + j]);
                    }
                }
            }
        }

        reduced[pivot_j] = true;

        if (pivot_i != pivot_j)
        {
            for (int j=0; j<n; j++)
            {
                temp = m[pivot_i * n + j];
                m[pivot_i * n + j] = m[pivot_j * n + j];
                m[pivot_j * n + j] = temp;
            }

            temp = rhs[pivot_i];
            rhs[pivot_i] = rhs[pivot_j];
            rhs[pivot_j] = temp;
        }

        if (m[pivot_j * n + pivot_j] == 0.0) 
        {
            throw("The matrix is singular.");
        }

        temp = 1.0/ m[pivot_j * n + pivot_j];
        for (int j=0; j<n; j++)
        {
            m[pivot_j * n + j] *= temp;
        }
        rhs[pivot_j] *= temp;

        for (int i=0; i<n; i++)
        {
            if (i != pivot_j)
            {
                temp = m[i * n + pivot_j];
                for (int j=0; j<n; j++) 
                {
                    m[i * n + j] -= m[pivot_j * n + j] * temp;
                }
                rhs[i] -= rhs[pivot_j] * temp;
            }
        }
    }

    delete[] reduced;
}

void fitPolynomial(double * x, double * y, const int n) 
{
    double * matrix = new double[n*n];
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<n; j++)
        {
            matrix[i*n + j] = pow(x[i], j);
        }
    }

    solveMatrixEquation(matrix, y, n);
}    

double calculatePolynomial(const double * coeffs, const int n, const double x)
{
    double output = 0;
    for (int i=0; i<n; i++)
    {
        output += coeffs[i] * pow(x, i);
    }
    return output;
}

}
