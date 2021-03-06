/*
 * Copyright (C) 2013-  Qiming Sun <osirpt.sun@gmail.com>
 *
 * Provide the intermediate variable g(nroots,i,j,k,l,[xyz])
 */
#pragma once
#include "cint_const.hpp"
#include "g1e.hpp"

#ifndef HAVE_BC
#define HAVE_BC
struct _BC {
        double c00[MXRYSROOTS*3];
        double c0p[MXRYSROOTS*3];
        double b01[MXRYSROOTS];
        double b00[MXRYSROOTS];
        double b10[MXRYSROOTS];
};
#endif

void CINTg2e_index_xyz(int *idx, const CINTEnvVars *envs);

void CINTinit_int2e_EnvVars(CINTEnvVars *envs, int *ng, int *shls,
                            int *atm, int natm, int *bas, int nbas, double *env);
void CINTinit_int3c2e_EnvVars(CINTEnvVars *envs, int *ng, int *shls,
                              int *atm, int natm, int *bas, int nbas, double *env);
void CINTinit_int2c2e_EnvVars(CINTEnvVars *envs, int *ng, int *shls,
                              int *atm, int natm, int *bas, int nbas, double *env);

int CINTg0_2e(double *g, const double fac, const CINTEnvVars *envs);
void CINTg0_2e_2d(double *g, struct _BC *bc, const CINTEnvVars *envs);
void CINTg0_2e_lj2d4d(double *g, struct _BC *bc, const CINTEnvVars *envs);
void CINTg0_2e_kj2d4d(double *g, struct _BC *bc, const CINTEnvVars *envs);
void CINTg0_2e_il2d4d(double *g, struct _BC *bc, const CINTEnvVars *envs);
void CINTg0_2e_ik2d4d(double *g, struct _BC *bc, const CINTEnvVars *envs);

void CINTg0_lj2d_4d(double *g, struct _BC *bc, const CINTEnvVars *envs);
void CINTg0_kj2d_4d(double *g, struct _BC *bc, const CINTEnvVars *envs);
void CINTg0_il2d_4d(double *g, struct _BC *bc, const CINTEnvVars *envs);
void CINTg0_ik2d_4d(double *g, struct _BC *bc, const CINTEnvVars *envs);

void CINTnabla1i_2e(double *f, const double *g,
                    const int li, const int lj, const int lk, const int ll,
                    const CINTEnvVars *envs);

void CINTnabla1j_2e(double *f, const double *g,
                    const int li, const int lj, const int lk, const int ll,
                    const CINTEnvVars *envs);

void CINTnabla1k_2e(double *f, const double *g,
                    const int li, const int lj, const int lk, const int ll,
                    const CINTEnvVars *envs);

void CINTnabla1l_2e(double *f, const double *g,
                    const int li, const int lj, const int lk, const int ll,
                    const CINTEnvVars *envs);

void CINTx1i_2e(double *f, const double *g, const double *ri,
                const int li, const int lj, const int lk, const int ll,
                const CINTEnvVars *envs);

void CINTx1j_2e(double *f, const double *g, const double *rj,
                const int li, const int lj, const int lk, const int ll,
                const CINTEnvVars *envs);

void CINTx1k_2e(double *f, const double *g, const double *rk,
                const int li, const int lj, const int lk, const int ll,
                const CINTEnvVars *envs);

void CINTx1l_2e(double *f, const double *g, const double *rl,
                const int li, const int lj, const int lk, const int ll,
                const CINTEnvVars *envs);

#ifdef WITH_F12
void CINTinit_int2e_stg_EnvVars(CINTEnvVars *envs, int *ng, int *shls,
                           int *atm, int natm, int *bas, int nbas, double *env);
void CINTinit_int2e_yp_EnvVars(CINTEnvVars *envs, int *ng, int *shls,
                           int *atm, int natm, int *bas, int nbas, double *env);
#endif

#ifdef WITH_GTG
void CINTinit_int2e_gtg_EnvVars(CINTEnvVars *envs, int *ng, int *shls,
                                int *atm, int natm, int *bas, int nbas, double *env);
void CINTinit_int3c2e_gtg_EnvVars(CINTEnvVars *envs, int *ng, int *shls,
                                  int *atm, int natm, int *bas, int nbas, double *env);
void CINTinit_int2c2e_gtg_EnvVars(CINTEnvVars *envs, int *ng, int *shls,
                                  int *atm, int natm, int *bas, int nbas, double *env);
#endif


#define G2E_D_I(f, g, li, lj, lk, ll)   CINTnabla1i_2e(f, g, li, lj, lk, ll, envs)
#define G2E_D_J(f, g, li, lj, lk, ll)   CINTnabla1j_2e(f, g, li, lj, lk, ll, envs)
#define G2E_D_K(f, g, li, lj, lk, ll)   CINTnabla1k_2e(f, g, li, lj, lk, ll, envs)
#define G2E_D_L(f, g, li, lj, lk, ll)   CINTnabla1l_2e(f, g, li, lj, lk, ll, envs)
/* r-R_0, R_0 is (0,0,0) */
#define G2E_R0I(f, g, li, lj, lk, ll)   CINTx1i_2e(f, g, envs->ri, li, lj, lk, ll, envs)
#define G2E_R0J(f, g, li, lj, lk, ll)   CINTx1j_2e(f, g, envs->rj, li, lj, lk, ll, envs)
#define G2E_R0K(f, g, li, lj, lk, ll)   CINTx1k_2e(f, g, envs->rk, li, lj, lk, ll, envs)
#define G2E_R0L(f, g, li, lj, lk, ll)   CINTx1l_2e(f, g, envs->rl, li, lj, lk, ll, envs)
/* r-R_C, R_C is common origin */
#define G2E_RCI(f, g, li, lj, lk, ll)   CINTx1i_2e(f, g, dri, li, lj, lk, ll, envs)
#define G2E_RCJ(f, g, li, lj, lk, ll)   CINTx1j_2e(f, g, drj, li, lj, lk, ll, envs)
#define G2E_RCK(f, g, li, lj, lk, ll)   CINTx1k_2e(f, g, drk, li, lj, lk, ll, envs)
#define G2E_RCL(f, g, li, lj, lk, ll)   CINTx1l_2e(f, g, drl, li, lj, lk, ll, envs)
/* origin from center of each basis
 * x1[ijkl]_2e(f, g, ng, li, lj, lk, ll, 0d0) */
#define G2E_R_I(f, g, li, lj, lk, ll)   f = g + envs->g_stride_i
#define G2E_R_K(f, g, li, lj, lk, ll)   f = g + envs->g_stride_k
#define G2E_R_L(f, g, li, lj, lk, ll)   f = g + envs->g_stride_l
#define G2E_R_J(f, g, li, lj, lk, ll)   f = g + envs->g_stride_j
