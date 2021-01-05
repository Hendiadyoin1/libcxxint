/*
 * Copyright (C) 2013  Qiming Sun <osirpt.sun@gmail.com>
 *
 */
#pragma once
#include <complex>
using namespace std::literals;
#include "config.hpp"

int CINT1e_loop(double *gctr, CINTEnvVars *envs);

int CINT1e_nuc_loop(double *gctr, CINTEnvVars *envs, double fac, int nuc_id);

int CINT1e_drv(double *out, int *dims, CINTEnvVars *envs,
void (*f_c2s)(double *opij, double *gctr, int *dims,
                 CINTEnvVars *envs), int int1e_type);

int CINT1e_spinor_drv(std::complex<double> *out, int *dims, CINTEnvVars *envs,
void (*f_c2s)(std::complex<double> *opij, double *gctr, int *dims,
               CINTEnvVars *envs), int int1e_type);

double CINTnuc_mod(double aij, int nuc_id, int *atm, double *env);

int CINT3c1e_spheric_drv(double *out, int *dims, CINTEnvVars *envs, CINTOpt *opt,
void (*f_e1_c2s)(double *out, double *gctr, int *dims,
                  CINTEnvVars *envs), int int_type, int is_ssc);
int CINT3c1e_cart_drv(double *out, int *dims, CINTEnvVars *envs, CINTOpt *opt,
int int_type);
int CINT3c1e_spinor_drv(std::complex<double> *out, int *dims, CINTEnvVars *envs, CINTOpt *opt,
void (*f_e1_c2s)(std::complex<double> *opijk, double *gctr, int *dims,
                  CINTEnvVars *envs), int int_type, int is_ssc);

#define INT1E_TYPE_OVLP 0
#define INT1E_TYPE_RINV 1
#define INT1E_TYPE_NUC  2

