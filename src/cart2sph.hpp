/*
 * Copyright (C) 2013-  Qiming Sun <osirpt.sun@gmail.com>
 *
 * Cartisen GTO to spheric or spinor GTO transformation
 */

/*************************************************
 *
 * transform matrix
 *
 *************************************************/
#pragma once
#include <complex>
using namespace std::literals;
#include "g1e.hpp"

void c2s_sph_1e(double *opij, double *gctr, int *dims, CINTEnvVars *envs);
void c2s_sph_2e1(double *fijkl, double *gctr, int *dims, CINTEnvVars *envs);
void c2s_sph_2e2();

void c2s_cart_1e(double *opij, double *gctr, int *dims, CINTEnvVars *envs);
void c2s_cart_2e1(double *fijkl, double *gctr, int *dims, CINTEnvVars *envs);
void c2s_cart_2e2();

void c2s_sf_1e(std::complex<double> *opij, double *gctr, int *dims, CINTEnvVars *envs);
void c2s_sf_1ei(std::complex<double> *opij, double *gctr, int *dims, CINTEnvVars *envs);

void c2s_si_1e(std::complex<double> *opij, double *gctr, int *dims, CINTEnvVars *envs);
void c2s_si_1ei(std::complex<double> *opij, double *gctr, int *dims, CINTEnvVars *envs);

void c2s_sf_2e1(std::complex<double> *opij, double *gctr, int *dims, CINTEnvVars *envs);
void c2s_sf_2e1i(std::complex<double> *opij, double *gctr, int *dims, CINTEnvVars *envs);

void c2s_sf_2e2(std::complex<double> *fijkl, std::complex<double> *opij, int *dims, CINTEnvVars *envs);
void c2s_sf_2e2i(std::complex<double> *fijkl, std::complex<double> *opij, int *dims, CINTEnvVars *envs);

void c2s_si_2e1(std::complex<double> *opij, double *gctr, int *dims, CINTEnvVars *envs);
void c2s_si_2e1i(std::complex<double> *opij, double *gctr, int *dims, CINTEnvVars *envs);

void c2s_si_2e2(std::complex<double> *fijkl, std::complex<double> *opij, int *dims, CINTEnvVars *envs);
void c2s_si_2e2i(std::complex<double> *fijkl, std::complex<double> *opij, int *dims, CINTEnvVars *envs);

void c2s_sph_3c2e1(double *fijkl, double *gctr, int *dims, CINTEnvVars *envs);
void c2s_cart_3c2e1(double *fijkl, double *gctr, int *dims, CINTEnvVars *envs);
void c2s_sph_3c2e1_ssc(double *fijkl, double *gctr, int *dims, CINTEnvVars *envs);

void c2s_sf_3c2e1(std::complex<double> *opijk, double *gctr, int *dims, CINTEnvVars *envs);
void c2s_sf_3c2e1i(std::complex<double> *opijk, double *gctr, int *dims, CINTEnvVars *envs);
void c2s_si_3c2e1(std::complex<double> *opijk, double *gctr, int *dims, CINTEnvVars *envs);
void c2s_si_3c2e1i(std::complex<double> *opijk, double *gctr, int *dims, CINTEnvVars *envs);
void c2s_sf_3c2e1_ssc(std::complex<double> *opijk, double *gctr, int *dims, CINTEnvVars *envs);
void c2s_sf_3c2e1i_ssc(std::complex<double> *opijk, double *gctr, int *dims, CINTEnvVars *envs);
void c2s_si_3c2e1_ssc(std::complex<double> *opijk, double *gctr, int *dims, CINTEnvVars *envs);
void c2s_si_3c2e1i_ssc(std::complex<double> *opijk, double *gctr, int *dims, CINTEnvVars *envs);

void c2s_sph_3c1e(double *fijkl, double *gctr, int *dims, CINTEnvVars *envs);
void c2s_cart_3c1e(double *fijkl, double *gctr, int *dims, CINTEnvVars *envs);

void c2s_dset0(double *out, int *dims, int *counts);
void c2s_zset0(std::complex<double> *out, int *dims, int *counts);

/*************************************************
 *
 * transform vectors
 *
 *************************************************/
void c2s_sph_vec(double *sph, double *cart, int l, int nvec);
