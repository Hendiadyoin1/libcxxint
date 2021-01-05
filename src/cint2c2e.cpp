/*
 * Copyright (C) 2013-  Qiming Sun <osirpt.sun@gmail.com>
 *
 * 2-center 2-electron integrals
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cint_bas.hpp"
#include "g2e.hpp"
#include "optimizer.hpp"
#include "cint2e.hpp"
#include "misc.hpp"
#include "cart2sph.hpp"
#include "c2f.hpp"

#define PRIM2CTR0(ctrsymb, gp, ngp) \
        if (ctrsymb##_ctr > 1) {\
                if (*ctrsymb##empty) { \
                        CINTprim_to_ctr_0(gctr##ctrsymb, gp, c##ctrsymb+ctrsymb##p, \
                                          ngp, ctrsymb##_prim, ctrsymb##_ctr, \
                                          non0ctr##ctrsymb[ctrsymb##p], \
                                          non0idx##ctrsymb+ctrsymb##p*ctrsymb##_ctr); \
                } else { \
                        CINTprim_to_ctr_1(gctr##ctrsymb, gp, c##ctrsymb+ctrsymb##p, \
                                          ngp, ctrsymb##_prim, ctrsymb##_ctr, \
                                          non0ctr##ctrsymb[ctrsymb##p], \
                                          non0idx##ctrsymb+ctrsymb##p*ctrsymb##_ctr); \
                } \
        } \
        *ctrsymb##empty = 0


int CINT2c2e_loop_nopt(double *gctr, CINTEnvVars *envs)
{
        int *shls  = envs->shls;
        int *bas = envs->bas;
        double *env = envs->env;
        int i_sh = shls[0];
        int k_sh = shls[1];
        int i_ctr = envs->x_ctr[0];
        int k_ctr = envs->x_ctr[1];
        int i_prim = bas(NPRIM_OF, i_sh);
        int k_prim = bas(NPRIM_OF, k_sh);
        double *ai = env + bas(PTR_EXP, i_sh);
        double *ak = env + bas(PTR_EXP, k_sh);
        double *ci = env + bas(PTR_COEFF, i_sh);
        double *ck = env + bas(PTR_COEFF, k_sh);
        int n_comp = envs->ncomp_tensor;
        double fac1i, fac1k;
        int ip, kp;
        int empty[3] = {1, 1, 1};
        int *iempty = empty + 0;
        int *kempty = empty + 1;
        int *gempty = empty + 2;
        /* COMMON_ENVS_AND_DECLARE end */
        const int nc = i_ctr * k_ctr;
        const int leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const int lenk = envs->nf * nc * n_comp; // gctrk
        const int leni = envs->nf * i_ctr * n_comp; // gctri
        const int len0 = envs->nf * n_comp; // gout
        const int len = leng + lenk + leni + len0;
        double *g = new double[len]; // only one allocation and then the pointer gets split up
        double *g1 = g + leng;
        double *gout, *gctri, *gctrk;

        if (n_comp == 1) {
                gctrk = gctr;
        } else {
                gctrk = g1;
                g1 += lenk;
        }
        if (k_ctr == 1) {
                gctri = gctrk;
                iempty = kempty;
        } else {
                gctri = g1;
                g1 += leni;
        }
        if (i_ctr == 1) {
                gout = gctri;
                gempty = iempty;
        } else {
                gout = g1;
        }

        envs->idx = new int[envs->nf * 3];
        CINTg1e_index_xyz(envs->idx, envs);
        int *non0idxi = new int[i_prim*i_ctr];
        int *non0idxk = new int[k_prim*k_ctr];
        int *non0ctri = new int[i_prim];
        int *non0ctrk = new int[k_prim];
        if (i_ctr > 1) {
                CINTOpt_non0coeff_byshell(non0idxi, non0ctri, ci, i_prim, i_ctr);
        }
        if (k_ctr > 1) {
                CINTOpt_non0coeff_byshell(non0idxk, non0ctrk, ck, k_prim, k_ctr);
        }

        *kempty = 1;
        for (kp = 0; kp < k_prim; kp++) {
                envs->ak = ak[kp];
                envs->akl = ak[kp]; // to use CINTg0_2e
                if (k_ctr == 1) {
                        fac1k = envs->common_factor * ck[kp];
                } else {
                        fac1k = envs->common_factor;
                        *iempty = 1;
                }
                for (ip = 0; ip < i_prim; ip++) {
                        envs->ai = ai[ip];
                        envs->aij = ai[ip];
                        if (i_ctr == 1) {
                                fac1i = fac1k*ci[ip];
                        } else {
                                fac1i = fac1k;
                        }
                        if ((*envs->f_g0_2e)(g, fac1i, envs)) {
                                (*envs->f_gout)(gout, g, envs->idx, envs, *gempty);
                                PRIM2CTR0(i, gout, envs->nf*n_comp);
                        }
                } // end loop i_prim
                if (!*iempty) {
                        PRIM2CTR0(k, gctri, envs->nf*i_ctr*n_comp);
                }
        } // end loop k_prim

        if (n_comp > 1 && !*kempty) {
                CINTdmat_transpose(gctr, gctrk, envs->nf*nc, n_comp);
        }
        delete[] envs->idx; envs->idx = nullptr;
        return !*kempty;
}


#define COMMON_ENVS_AND_DECLARE \
        int *shls = envs->shls; \
        int *bas = envs->bas; \
        double *env = envs->env; \
        int i_sh = shls[0]; \
        int k_sh = shls[1]; \
        int i_ctr = envs->x_ctr[0]; \
        int k_ctr = envs->x_ctr[1]; \
        int i_prim = bas(NPRIM_OF, i_sh); \
        int k_prim = bas(NPRIM_OF, k_sh); \
        double *ai = env + bas(PTR_EXP, i_sh); \
        double *ak = env + bas(PTR_EXP, k_sh); \
        double *ci = env + bas(PTR_COEFF, i_sh); \
        double *ck = env + bas(PTR_COEFF, k_sh); \
        int n_comp = envs->ncomp_tensor; \
        double fac1i, fac1k; \
        int ip, kp; \
        int empty[3] = {1, 1, 1}; \
        int *iempty = empty + 0; \
        int *kempty = empty + 1; \
        int *gempty = empty + 2; \
        int *non0idxi = new int[i_prim*i_ctr]; \
        int *non0idxk = new int[k_prim*k_ctr]; \
        int *non0ctri = new int[i_prim]; \
        int *non0ctrk = new int[k_prim]; \
        if (i_ctr > 1) { \
                CINTOpt_non0coeff_byshell(non0idxi, non0ctri, ci, i_prim, i_ctr); \
        } \
        if (k_ctr > 1) { \
                CINTOpt_non0coeff_byshell(non0idxk, non0ctrk, ck, k_prim, k_ctr); \
        }

#define USE_OPT \
        envs->idx = opt->index_xyz_array[envs->i_l*LMAX1+envs->k_l]

#define PRIM2CTR(ctrsymb, gp, ngp) \
        if (ctrsymb##_ctr > 1) {\
                if (*ctrsymb##empty) { \
                        CINTprim_to_ctr_0(gctr##ctrsymb, gp, c##ctrsymb+ctrsymb##p, \
                                          ngp, ctrsymb##_prim, ctrsymb##_ctr, \
                                          non0ctr##ctrsymb[ctrsymb##p], \
                                          non0idx##ctrsymb+ctrsymb##p*ctrsymb##_ctr); \
                } else { \
                        CINTprim_to_ctr_1(gctr##ctrsymb, gp, c##ctrsymb+ctrsymb##p, \
                                          ngp, ctrsymb##_prim, ctrsymb##_ctr, \
                                          non0ctr##ctrsymb[ctrsymb##p], \
                                          non0idx##ctrsymb+ctrsymb##p*ctrsymb##_ctr); \
                } \
        } \
        *ctrsymb##empty = 0

// i_ctr = k_ctr = 1;
int CINT2c2e_11_loop(double *gctr, CINTEnvVars *envs, const CINTOpt *opt)
{
        COMMON_ENVS_AND_DECLARE;
        const int nc = 1;
        const int leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const int len0 = envs->nf * n_comp;
        const int len = leng + len0;
        double *g = new double[len];
        double *gout;
        if (n_comp == 1) {
                gout = gctr;
        } else {
                gout = g + leng;
        }

        USE_OPT;

        for (kp = 0; kp < k_prim; kp++) {
                envs->ak = ak[kp];
                envs->akl = ak[kp];
                fac1k = envs->common_factor * ck[kp];
                for (ip = 0; ip < i_prim; ip++) {
                        envs->ai = ai[ip];
                        envs->aij = ai[ip];
                        fac1i = fac1k*ci[ip];
                        if ((*envs->f_g0_2e)(g, fac1i, envs)) {
                                (*envs->f_gout)(gout, g, envs->idx, envs, *empty);
                                *empty = 0;
                        }
                } // end loop i_prim
        } // end loop k_prim

        if (n_comp > 1 && !*empty) {
                CINTdmat_transpose(gctr, gout, envs->nf*nc, n_comp);
        }
        return !*empty;
}

// i_ctr = n; k_ctr = 1;
int CINT2c2e_n1_loop(double *gctr, CINTEnvVars *envs, const CINTOpt *opt)
{
        COMMON_ENVS_AND_DECLARE;

        const int nc = i_ctr;
        const int leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const int leni = envs->nf * i_ctr * n_comp; // gctri
        const int len0 = envs->nf * n_comp; // gout
        const int len = leng + leni + len0;
        double *g= new double[len];
        double *g1 = g + leng;
        double *gout, *gctri;
        if (n_comp == 1) {
                gctri = gctr;
        } else {
                gctri = g1;
                g1 += leni;
        }
        gout = g1;

        USE_OPT;

        for (kp = 0; kp < k_prim; kp++) {
                envs->ak = ak[kp];
                envs->akl = ak[kp];
                fac1k = envs->common_factor * ck[kp];
                for (ip = 0; ip < i_prim; ip++) {
                        envs->ai = ai[ip];
                        envs->aij = ai[ip];
                        fac1i = fac1k;
                        if ((*envs->f_g0_2e)(g, fac1i, envs)) {
                                (*envs->f_gout)(gout, g, envs->idx, envs, 1);
                                PRIM2CTR(i, gout, envs->nf*n_comp);
                        }
                } // end loop i_prim
        } // end loop k_prim

        if (n_comp > 1 && !*iempty) {
                CINTdmat_transpose(gctr, gctri, envs->nf*nc, n_comp);
        }
        return !*iempty;
}

// k_ctr = n; i_ctr = 1;
int CINT2c2e_1n_loop(double *gctr, CINTEnvVars *envs, const CINTOpt *opt)
{
        COMMON_ENVS_AND_DECLARE;

        const int nc = k_ctr;
        const int leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const int lenk = envs->nf * k_ctr * n_comp; // gctrk
        const int len0 = envs->nf * n_comp; // gout
        const int len = leng + lenk + len0;
        double *g= new double[len];
        double *g1 = g + leng;
        double *gout, *gctrk;
        if (n_comp == 1) {
                gctrk = gctr;
        } else {
                gctrk = g1;
                g1 += lenk;
        }
        gout = g1;

        USE_OPT;

        for (kp = 0; kp < k_prim; kp++) {
                envs->ak = ak[kp];
                envs->akl = ak[kp];
                fac1k = envs->common_factor;
                *iempty = 1;
                for (ip = 0; ip < i_prim; ip++) {
                        envs->ai = ai[ip];
                        envs->aij = ai[ip];
                        fac1i = fac1k*ci[ip];
                        if ((*envs->f_g0_2e)(g, fac1i, envs)) {
                                (*envs->f_gout)(gout, g, envs->idx, envs, *iempty);
                                *iempty = 0;
                        }
                } // end loop i_prim
                if (!*iempty) {
                        PRIM2CTR(k, gout,envs->nf*n_comp);
                }
        } // end loop k_prim

        if (n_comp > 1 && !*kempty) {
                CINTdmat_transpose(gctr, gctrk, envs->nf*nc, n_comp);
        }
        return !*kempty;
}


int CINT2c2e_loop(double *gctr, CINTEnvVars *envs, const CINTOpt *opt)
{
        COMMON_ENVS_AND_DECLARE;
        const int nc = i_ctr * k_ctr;
        const int leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const int lenk = envs->nf * nc * n_comp; // gctrk
        const int leni = envs->nf * i_ctr * n_comp; // gctri
        const int len0 = envs->nf * n_comp; // gout
        const int len = leng + lenk + leni + len0;
        double *g = new double[len];
        double *g1 = g + leng;
        double *gout, *gctri, *gctrk;

        if (n_comp == 1) {
                gctrk = gctr;
        } else {
                gctrk = g1;
                g1 += lenk;
        }
        if (k_ctr == 1) {
                gctri = gctrk;
                iempty = kempty;
        } else {
                gctri = g1;
                g1 += leni;
        }
        if (i_ctr == 1) {
                gout = gctri;
                gempty = iempty;
        } else {
                gout = g1;
        }

        /* USE_OPT */
        envs->idx = opt->index_xyz_array[envs->i_l*LMAX1+envs->k_l];
        /* USE_OPT end */

        *kempty = 1;
        for (kp = 0; kp < k_prim; kp++) {
                envs->ak = ak[kp];
                envs->akl = ak[kp];
                if (k_ctr == 1) {
                        fac1k = envs->common_factor * ck[kp];
                } else {
                        fac1k = envs->common_factor;
                        *iempty = 1;
                }
                for (ip = 0; ip < i_prim; ip++) {
                        envs->ai = ai[ip];
                        envs->aij = ai[ip];
                        if (i_ctr == 1) {
                                fac1i = fac1k*ci[ip];
                        } else {
                                fac1i = fac1k;
                        }
                        if ((*envs->f_g0_2e)(g, fac1i, envs)) {
                                (*envs->f_gout)(gout, g, envs->idx, envs, *gempty);
                                PRIM2CTR(i, gout, envs->nf*n_comp);
                        }
                } // end loop i_prim
                if (!*iempty) {
                        PRIM2CTR(k, gctri, envs->nf*i_ctr*n_comp);
                }
        } // end loop k_prim

        if (n_comp > 1 && !*kempty) {
                CINTdmat_transpose(gctr, gctrk, envs->nf*nc, n_comp);
        }
        return !*kempty;
}

static int (*CINTf_2c2e_loop[8])(double *gctr, CINTEnvVars *envs, const CINTOpt *opt) = {
        CINT2c2e_loop,
        CINT2c2e_n1_loop,
        CINT2c2e_1n_loop,
        CINT2c2e_11_loop,
};

int CINT2c2e_cart_drv(double *out, int *dims, CINTEnvVars *envs, CINTOpt *opt)
{
        int *x_ctr = envs->x_ctr;
        int nc = envs->nf * x_ctr[0] * x_ctr[1];
        int n_comp = envs->ncomp_e1 * envs->ncomp_e2 * envs->ncomp_tensor;
        double *gctr=new double[nc*n_comp];

        int n;
        int has_value;

        if (opt != NULL) {
                n = ((envs->x_ctr[0]==1) << 1) + (envs->x_ctr[1]==1);
                has_value = CINTf_2c2e_loop[n](gctr, envs, opt);
        } else {
                has_value = CINT2c2e_loop_nopt(gctr, envs);
        }

        int counts[4];
        counts[0] = envs->nfi * x_ctr[0];
        counts[1] = envs->nfk * x_ctr[1];
        counts[2] = 1;
        counts[3] = 1;
        if (dims == NULL) {
                dims = counts;
        }
        int nout = dims[0] * dims[1];
        if (has_value) {
                for (n = 0; n < n_comp; n++) {
                        c2s_cart_1e(out+nout*n, gctr+nc*n, dims, envs);
                }
        } else {
                for (n = 0; n < n_comp; n++) {
                        c2s_dset0(out+nout*n, dims, counts);
                }
        }
        return has_value;
}
int CINT2c2e_spheric_drv(double *out, int *dims, CINTEnvVars *envs, CINTOpt *opt)
{
        int *x_ctr = envs->x_ctr;
        int nc = envs->nf * x_ctr[0] * x_ctr[1];
        int n_comp = envs->ncomp_e1 * envs->ncomp_e2 * envs->ncomp_tensor;
        double *gctr = new double[nc*n_comp];

        int n;
        int has_value;

        if (opt != NULL) {
                n = ((envs->x_ctr[0]==1) << 1) + (envs->x_ctr[1]==1);
                has_value = CINTf_2c2e_loop[n](gctr, envs, opt);
        } else {
                has_value = CINT2c2e_loop_nopt(gctr, envs);
        }

        int counts[4];
        counts[0] = (envs->i_l*2+1) * x_ctr[0];
        counts[1] = (envs->k_l*2+1) * x_ctr[1];
        counts[2] = 1;
        counts[3] = 1;
        if (dims == NULL) {
                dims = counts;
        }
        int nout = dims[0] * dims[1];
        if (has_value) {
                for (n = 0; n < n_comp; n++) {
                        c2s_sph_1e(out+nout*n, gctr+nc*n, dims, envs);
                }
        } else {
                for (n = 0; n < n_comp; n++) {
                        c2s_dset0(out+nout*n, dims, counts);
                }
        }
        return has_value;
}
// (spinor|spinor)
int CINT2c2e_spinor_drv(std::complex<double> *out, int *dims, CINTEnvVars *envs, CINTOpt *opt,
void (*f_e1_c2s)(std::complex<double> *opij, double *gctr, int *dims,CINTEnvVars *envs))
{
        if (envs->ncomp_e1 > 1 || envs->ncomp_e2 > 1) {
                fprintf(stderr, "CINT2c2e_spinor_drv not implemented\n");
                exit(1);
        }
        int *x_ctr = envs->x_ctr;
        int counts[4];
        counts[0] = CINTcgto_spinor(envs->shls[0], envs->bas);
        counts[1] = CINTcgto_spinor(envs->shls[1], envs->bas);
        counts[2] = 1;
        counts[3] = 1;
        int nc = envs->nf * x_ctr[0] * x_ctr[1];
        int n_comp = envs->ncomp_e1 * envs->ncomp_e2 * envs->ncomp_tensor;
        double *gctr = new double[nc*n_comp];

        int n, has_value;

        if (opt != NULL) {
                n = ((envs->x_ctr[0]==1) << 1) + (envs->x_ctr[1]==1);
                has_value = CINTf_2c2e_loop[n](gctr, envs, opt);
        } else {
                has_value = CINT2c2e_loop_nopt(gctr, envs);
        }

        if (dims == NULL) {
                dims = counts;
        }
        int nout = dims[0] * dims[1];
        if (has_value) {
                for (n = 0; n < n_comp; n++) {
                        (*f_e1_c2s)(out+nout*n, gctr, dims, envs);
                        gctr += nc;
                }
        } else {
                for (n = 0; n < n_comp; n++) {
                        c2s_zset0(out+nout*n, dims, counts);
                }
        }
        return has_value;
}


int int2c2e_sph(double *out, int *dims, int *shls, int *atm, int natm,
                int *bas, int nbas, double *env, CINTOpt *opt)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2c2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        return CINT2c2e_spheric_drv(out, dims, &envs, opt);
}
void int2c2e_optimizer(CINTOpt **opt, int *atm, int natm,
                       int *bas, int nbas, double *env)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTall_2c2e_optimizer(opt, ng, atm, natm, bas, nbas, env);
}

int int2c2e_cart(double *out, int *dims, int *shls, int *atm, int natm,
                 int *bas, int nbas, double *env, CINTOpt *opt)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2c2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        return CINT2c2e_cart_drv(out, dims, &envs, opt);
}
 
int int2c2e_spinor(std::complex<double> *out, int *dims, int *shls, int *atm, int natm,
                   int *bas, int nbas, double *env, CINTOpt *opt)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2c2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        return CINT2c2e_spinor_drv(out, dims, &envs, opt, &c2s_sf_1e);
}


ALL_CINT(int2c2e)
ALL_CINT_FORTRAN_(int2c2e)

