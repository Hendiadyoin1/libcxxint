/*
 * Copyright (C) 2013  Qiming Sun <osirpt.sun@gmail.com>
 *
 * basic functions
 */

#include <math.h>
#include <complex>
using namespace std::literals;
#include "misc.hpp"
#include "cint_const.hpp"

void CINTdcmplx_re(const int n, std::complex<double> *z, const double *re)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = re[i] + 0i;
        }
}

void CINTdcmplx_im(const int n, std::complex<double> *z, const double *im)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = im[i]*1i;
        }
}

void CINTdcmplx_pp(const int n, std::complex<double> *z,
                   const double *re, const double *im)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = re[i] + im[i] * 1i;
        }
}
void CINTdcmplx_pn(const int n, std::complex<double> *z,
                   const double *re, const double *im)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = re[i] - im[i] * 1i;
        }
}
void CINTdcmplx_np(const int n, std::complex<double> *z,
                   const double *re, const double *im)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = -re[i] + im[i] * 1i;
        }
}
void CINTdcmplx_nn(const int n, std::complex<double> *z,
                   const double *re, const double *im)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = -re[i] - im[i] * 1i;
        }
}


double CINTsquare_dist(const double *r1, const double *r2)
{
        double r12[3];

        r12[0] = r1[0] - r2[0];
        r12[1] = r1[1] - r2[1];
        r12[2] = r1[2] - r2[2];

        return r12[0] * r12[0] + r12[1] * r12[1] + r12[2] * r12[2];
}

static int factorial(int n)
{
        int i, fact = 1;
        for (i = 1; i <= n; i++) {
                fact *= i;
        }
        return fact;
}
double CINTgto_norm(int n, double a)
{
        double nn = pow(2, (2*n+3)) * factorial(n+1) * pow((2*a), (n+1.5)) \
                / (factorial(2*n+2) * sqrt(M_PI));
        return sqrt(nn);
}
double CINTgto_norm_(int *n, double *a)
{
        return CINTgto_norm(*n, *a);
}
