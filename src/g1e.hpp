/*
 * Copyright (C) 2013  Qiming Sun <osirpt.sun@gmail.com>
 *
 */
#pragma once
#include "config.hpp"

#if !defined HAVE_DEFINED_CINTENVVARS_H
#define HAVE_DEFINED_CINTENVVARS_H
// ref to CINTinit_int1e_EnvVars, CINTinit_int2e_EnvVars
struct CINTEnvVars{
        int *atm;
        int *bas;
        double *env;
        int *shls;
        int natm;
        int nbas;

        int i_l;
        int j_l;
        int k_l;
        int l_l;
        int nfi;  // number of cartesion components
        int nfj;
        int nfk;
        int nfl;
        int nf;  // = nfi*nfj*nfk*nfl;
        int _padding;
        int x_ctr[4];

        int gbits;
        int ncomp_e1; // = 1 if spin free, = 4 when spin included, it
        int ncomp_e2; // corresponds to POSX,POSY,POSZ,POS1, see cint_const.h
        int ncomp_tensor; // e.g. = 3 for gradients

        /* values may diff based on the g0_2d4d algorithm */
        int li_ceil; // power of x, == i_l if nabla is involved, otherwise == i_l
        int lj_ceil;
        int lk_ceil;
        int ll_ceil;
        int g_stride_i; // nrys_roots * shift of (i++,k,l,j)
        int g_stride_k; // nrys_roots * shift of (i,k++,l,j)
        int g_stride_l; // nrys_roots * shift of (i,k,l++,j)
        int g_stride_j; // nrys_roots * shift of (i,k,l,j++)
        int nrys_roots;
        int g_size;  // ref to cint2e.c g = malloc(sizeof(double)*g_size)

        int g2d_ijmax;
        int g2d_klmax;
        double common_factor;
        double expcutoff;
        double rirj[3]; // diff by sign in different g0_2d4d algorithm
        double rkrl[3];
        double *rx_in_rijrx;
        double *rx_in_rklrx;

        double *ri;
        double *rj;
        double *rk;
        double *rl;

        int  (*f_g0_2e)(double *g, const double fac, const CINTEnvVars *envs);
        void (*f_g0_2d4d)(double *g, struct _BC *bc, const CINTEnvVars *envs);
        void (*f_gout)(double *gout, double *g, int *idx, CINTEnvVars *envs, int gout_empty);

        /* values are assigned during calculation */
        int *idx;
        double ai;
        double aj;
        double ak;
        double al;
        double rij[3];
        double rijrx[3];
        double aij;
        double rkl[3];
        double rklrx[3];
        double akl;
};
#endif

void CINTinit_int1e_EnvVars(CINTEnvVars *envs, int *ng, int *shls,
                            int *atm, int natm, int *bas, int nbas, double *env);
void CINTinit_int3c1e_EnvVars(CINTEnvVars *envs, int *ng, int *shls,
                              int *atm, int natm, int *bas, int nbas, double *env);

void CINTg1e_index_xyz(int *idx,const CINTEnvVars *envs);

void CINTg_ovlp(double *g, double ai, double aj, double fac, CINTEnvVars *envs);

void CINTg_nuc(double *g, double aij, double *rij,
               double *cr, double t2, double fac, CINTEnvVars *envs);

void CINTnabla1i_1e(double *f, double *g,
                    int li, int lj, int lk, CINTEnvVars *envs);

void CINTnabla1j_1e(double *f, double *g,
                    int li, int lj, int lk, CINTEnvVars *envs);

void CINTnabla1k_1e(double *f, double *g,
                    int li, int lj, int lk, CINTEnvVars *envs);

void CINTx1i_1e(double *f, double *g, double ri[3],
                int li, int lj, int lk, CINTEnvVars *envs);

void CINTx1j_1e(double *f, double *g, double rj[3],
                int li, int lj, int lk, CINTEnvVars *envs);

void CINTx1k_1e(double *f, double *g, double rk[3],
                int li, int lj, int lk, CINTEnvVars *envs);

void CINTprim_to_ctr(double *gc, int nf, double *gp,
                     int inc, int nprim,
                     int nctr, double *pcoeff);

double CINTcommon_fac_sp(int l);

void CINTprim_to_ctr_0(double *gc, double *gp, double *coeff, int nf,
                       int nprim, int nctr, int non0ctr, int *sortedidx); // TODO change last attribute back to int*
void CINTprim_to_ctr_1(double *gc, double *gp, double *coeff, int nf,
                       int nprim, int nctr, int non0ctr, int *sortedidx); // TODO change last attribute back to int*

#define G1E_D_I(f, g, li, lj, lk)   CINTnabla1i_1e(f, g, li, lj, lk, envs)
#define G1E_D_J(f, g, li, lj, lk)   CINTnabla1j_1e(f, g, li, lj, lk, envs)
#define G1E_D_K(f, g, li, lj, lk)   CINTnabla1k_1e(f, g, li, lj, lk, envs)
/* r-R_0, R_0 is (0,0,0) */
#define G1E_R0I(f, g, li, lj, lk)   CINTx1i_1e(f, g, envs->ri, li, lj, lk, envs)
#define G1E_R0J(f, g, li, lj, lk)   CINTx1j_1e(f, g, envs->rj, li, lj, lk, envs)
#define G1E_R0K(f, g, li, lj, lk)   CINTx1k_1e(f, g, envs->rk, li, lj, lk, envs)
/* r-R_C, R_C is common origin */
#define G1E_RCI(f, g, li, lj, lk)   CINTx1i_1e(f, g, dri, li, lj, lk, envs)
#define G1E_RCJ(f, g, li, lj, lk)   CINTx1j_1e(f, g, drj, li, lj, lk, envs)
#define G1E_RCK(f, g, li, lj, lk)   CINTx1k_1e(f, g, drk, li, lj, lk, envs)
/* origin from center of each basis
 * x1[ij]_1e(f, g, ng, li, lj, 0d0) */
#define G1E_R_I(f, g, li, lj, lk)   f = g + envs->g_stride_i
#define G1E_R_J(f, g, li, lj, lk)   f = g + envs->g_stride_j
#define G1E_R_K(f, g, li, lj, lk)   f = g + envs->g_stride_k
