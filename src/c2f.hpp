/*
 * Copyright (C) 2013  Qiming Sun <osirpt.sun@gmail.com>
 *
 */
#pragma once
#ifdef WITH_FORTRAN
#include "config.hpp"

#define ALL_CINT_FORTRAN_(NAME) \
int c##NAME##_sph_(double *out, int *shls, int *atm, int *natm, \
                    int *bas, int *nbas, double *env, size_t optptr_as_integer8) { \
        CINTOpt **opt = (CINTOpt **)optptr_as_integer8; \
        return NAME##_sph(out, NULL, shls, \
                          atm, *natm, bas, *nbas, env, *opt); \
} \
void c##NAME##_sph_optimizer_(size_t optptr_as_integer8, int *atm, int *natm, \
                              int *bas, int *nbas, double *env) { \
        CINTOpt **opt = (CINTOpt **)optptr_as_integer8; \
        NAME##_optimizer(opt, atm, *natm, bas, *nbas, env); \
} \
int c##NAME##_cart_(double *out, int *shls, int *atm, int *natm, \
                     int *bas, int *nbas, double *env, size_t optptr_as_integer8) { \
        CINTOpt **opt = (CINTOpt **)optptr_as_integer8; \
        return NAME##_cart(out, NULL, shls, \
                           atm, *natm, bas, *nbas, env, *opt); \
} \
void c##NAME##_cart_optimizer_(CINTOpt **opt, int *atm, int *natm, \
                               int *bas, int *nbas, double *env) { \
        NAME##_optimizer(opt, atm, *natm, bas, *nbas, env); \
} \
int c##NAME##_(double *out, int *shls, int *atm, int *natm, \
                int *bas, int *nbas, double *env, size_t optptr_as_integer8) { \
        CINTOpt **opt = (CINTOpt **)optptr_as_integer8; \
        return NAME##_spinor((std::complex<double> *)out, NULL, shls, \
                             atm, *natm, bas, *nbas, env, *opt); \
} \
void c##NAME##_optimizer_(size_t optptr_as_integer8, int *atm, int *natm, \
                         int *bas, int *nbas, double *env) { \
        CINTOpt **opt = (CINTOpt **)optptr_as_integer8; \
        NAME##_optimizer(opt, atm, *natm, bas, *nbas, env); \
}

#define ALL_CINT1E_FORTRAN_(NAME) \
int c##NAME##_sph_(double *out, int *shls, int *atm, int *natm, \
                    int *bas, int *nbas, double *env) { \
        return NAME##_sph(out, NULL, shls, atm, *natm, bas, *nbas, env, NULL); \
} \
int c##NAME##_cart_(double *out, int *shls, int *atm, int *natm, \
                     int *bas, int *nbas, double *env) { \
        return NAME##_cart(out, NULL, shls, \
                           atm, *natm, bas, *nbas, env, NULL); \
} \
int c##NAME##_(double *out, int *shls, int *atm, int *natm, \
                int *bas, int *nbas, double *env) { \
        return NAME##_spinor((std::complex<double> *)out, NULL, shls, \
                             atm, *natm, bas, *nbas, env, NULL); \
}

#else

#define ALL_CINT_FORTRAN_(NAME)
#define ALL_CINT1E_FORTRAN_(NAME)

#endif
