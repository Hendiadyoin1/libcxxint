/*
 * Copyright (C) 2013  Qiming Sun <osirpt.sun@gmail.com>
 *
 * basic functions
 */
#pragma once
#include <stdint.h>
#include "config.hpp"
#include "fblas.hpp"

#define MIN(X,Y) ((X)<(Y)?(X):(Y))
#define MAX(X,Y) ((X)>(Y)?(X):(Y))
#define SQUARE(r) ((r)[0]*(r)[0] + (r)[1]*(r)[1] + (r)[2]*(r)[2])



void CINTdcmplx_re(const int n, std::complex<double> *z, const double *re);
void CINTdcmplx_im(const int n, std::complex<double> *z, const double *im);
void CINTdcmplx_pp(const int n, std::complex<double> *z, const double *re, const double *im);
void CINTdcmplx_pn(const int n, std::complex<double> *z, const double *re, const double *im);
void CINTdcmplx_np(const int n, std::complex<double> *z, const double *re, const double *im);
void CINTdcmplx_nn(const int n, std::complex<double> *z, const double *re, const double *im);

double CINTsquare_dist(const double *r1, const double *r2);

double CINTgto_norm(int n, double a);

#ifdef WITH_CINT2_INTERFACE
#define ALL_CINT(NAME) \
int c##NAME##_cart(double *out, int *shls, int *atm, int natm, \
            int *bas, int nbas, double *env, CINTOpt *opt) { \
        return NAME##_cart(out, NULL, shls, atm, natm, bas, nbas, env, opt); \
} \
void c##NAME##_cart_optimizer(CINTOpt **opt, int *atm, int natm, \
                         int *bas, int nbas, double *env) { \
        NAME##_optimizer(opt, atm, natm, bas, nbas, env); \
} \
int c##NAME##_sph(double *out, int *shls, int *atm, int natm, \
            int *bas, int nbas, double *env, CINTOpt *opt) { \
        return NAME##_sph(out, NULL, shls, atm, natm, bas, nbas, env, opt); \
} \
void c##NAME##_sph_optimizer(CINTOpt **opt, int *atm, int natm, \
                         int *bas, int nbas, double *env) { \
        NAME##_optimizer(opt, atm, natm, bas, nbas, env); \
} \
int c##NAME(double *out, int *shls, int *atm, int natm, \
            int *bas, int nbas, double *env, CINTOpt *opt) { \
        return NAME##_spinor((std::complex<double> *)out, NULL, shls, \
                             atm, natm, bas, nbas, env, opt); \
} \
void c##NAME##_optimizer(CINTOpt **opt, int *atm, int natm, \
                         int *bas, int nbas, double *env) { \
        NAME##_optimizer(opt, atm, natm, bas, nbas, env); \
}


#define ALL_CINT1E(NAME) \
int c##NAME##_cart(double *out, int *shls, int *atm, int natm, \
            int *bas, int nbas, double *env) { \
        return NAME##_cart(out, NULL, shls, atm, natm, bas, nbas, env, NULL); \
} \
int c##NAME##_sph(double *out, int *shls, int *atm, int natm, \
            int *bas, int nbas, double *env) { \
        return NAME##_sph(out, NULL, shls, atm, natm, bas, nbas, env, NULL); \
} \
int c##NAME(double *out, int *shls, int *atm, int natm, \
            int *bas, int nbas, double *env) { \
        return NAME##_spinor((std::complex<double> *)out, NULL, shls, \
                             atm, natm, bas, nbas, env, NULL); \
}

#else

#define ALL_CINT(NAME)
#define ALL_CINT1E(NAME)

#endif  // WITH_CINT2_INTERFACE
