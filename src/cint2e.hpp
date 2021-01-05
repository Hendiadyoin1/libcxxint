/*
 * Copyright (C) 2013-  Qiming Sun <osirpt.sun@gmail.com>
 *
 */
#pragma once
#include <complex>
using namespace std::literals;
#include "g1e.hpp"
#include "config.hpp"

void CINTgout2e(double *g, double *gout, int *idx,
                CINTEnvVars *envs, int gout_empty);

int CINT2e_loop(double *gctr, CINTEnvVars *envs, CINTOpt *opt);

int CINT2e_cart_drv(double *out, int *dims, CINTEnvVars *envs, CINTOpt *opt);
int CINT2e_spheric_drv(double *out, int *dims, CINTEnvVars *envs, CINTOpt *opt);
int CINT2e_spinor_drv(std::complex<double> *out, int *dims, CINTEnvVars *envs, CINTOpt *opt,
void (*f_e1_c2s)(std::complex<double> *opij, double *gctr, int *dims,
                 CINTEnvVars *envs), void (*f_e2_c2s)(std::complex<double> *fijkl, std::complex<double> *opij, int *dims,
                 CINTEnvVars *envs));

int CINT3c2e_cart_drv(double *out, int *dims, CINTEnvVars *envs, CINTOpt *opt);
int CINT3c2e_spheric_drv(double *out, int *dims, CINTEnvVars *envs, CINTOpt *opt,
void (*f_e1_c2s)(double *bufijk, double *gctr, int *dims,
                   CINTEnvVars *envs), int is_ssc);
int CINT3c2e_spinor_drv(std::complex<double> *out, int *dims, CINTEnvVars *envs, CINTOpt *opt,
void (*f_e1_c2s)(std::complex<double> *opijk, double *gctr, int *dims,
                  CINTEnvVars *envs), int is_ssc);
int CINT2c2e_cart_drv(double *out, int *dims, CINTEnvVars *envs, CINTOpt *opt);
int CINT2c2e_spheric_drv(double *out, int *dims, CINTEnvVars *envs, CINTOpt *opt);
int CINT2c2e_spinor_drv(std::complex<double> *out, int *dims, CINTEnvVars *envs, CINTOpt *opt,
    void (*f_e1_c2s)(std::complex<double> *opij, double *gctr, int *dims,CINTEnvVars *envs));
