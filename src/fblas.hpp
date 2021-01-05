/*
 * Copyright (C) 2013  Qiming Sun <osirpt.sun@gmail.com>
 *
 * blas interface and blas-like functions
 */
#pragma once
#include <complex>

using namespace std::literals;

#include "config.hpp"
double dasum_(const int *n, const double *dx, const int *incx);
void dscal_(const int *n, const double *da, double *dx, const int *incx);
void daxpy_(const int *n, const double *da, const double *dx,
           const int *incx, double *dy, const int *incy);
double ddot_(const int *n, const double *dx, const int *incx,
             const double *dy, const int *incy);
void dcopy_(const int *n, const double *dx, const int *incx,
            const double *dy, const int *incy);
void dgemm_(const char*, const char*,
            const int*, const int*, const int*,
            const double*, const double*, const int*,
            const double*, const int*,
            const double*, double*, const int*);
void dgemv_(const char*, const int*, const int*,
            const double*, const double*, const int*,
            const double*, const int*,
            const double*, double*, const int*);
void dger_(const int *m, const int *n,
           const double *alpha, const double *x,
           const int *incx, const double *y, const int *incy,
           double *a, const int *lda);
void dsymm_(const char*, const char*, const int*, const int*,
            const double*, const double*, const int*,
            const double*, const int*,
            const double*, double*, const int*);

//void dsyrk_
void zgerc_(const int *m, const int *n,
            const std::complex<double> *alpha, const std::complex<double> *x, const int *incx,
            const std::complex<double> *y, const int *incy,
            std::complex<double> *a, const int *lda);
void zgemv_(const char*, const int*, const int*,
            const std::complex<double>*, const std::complex<double>*, const int*,
            const std::complex<double>*, const int*,
            const std::complex<double>*, std::complex<double>*, const int*);
void zgemm_(const char*, const char*,
            const int*, const int*, const int*,
            const std::complex<double>*, const std::complex<double>*, const int*,
            const std::complex<double>*, const int*,
            const std::complex<double>*, std::complex<double>*, const int*);

void CINTdaxpy2v(const int n, const double a,
                 const double *x, const double *y, double *v);
void CINTdmat_transpose(double *a_t, const double *a, const int m, const int n);
void CINTzmat_transpose(std::complex<double> *a_t, const std::complex<double> *a,
                        const int m, const int n);
void CINTzmat_dagger(std::complex<double> *a_c, const std::complex<double> *a,
                     const int m, const int n);

