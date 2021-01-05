/*
 * Copyright (C) 2013-  Qiming Sun <osirpt.sun@gmail.com>
 *
 * 3-center 2-electron integrals
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cint_bas.hpp"
#include "optimizer.hpp"
#include "g2e.hpp"
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


int CINT3c2e_loop_nopt(double *gctr, CINTEnvVars *envs)
{
        int *shls  = envs->shls;
        int *bas = envs->bas;
        double *env = envs->env;
        int i_sh = shls[0];
        int j_sh = shls[1];
        int k_sh = shls[2];
        int i_ctr = envs->x_ctr[0];
        int j_ctr = envs->x_ctr[1];
        int k_ctr = envs->x_ctr[2];
        int i_prim = bas(NPRIM_OF, i_sh);
        int j_prim = bas(NPRIM_OF, j_sh);
        int k_prim = bas(NPRIM_OF, k_sh);
        //double *ri = envs->ri;
        //double *rj = envs->rj;
        double *ai = env + bas(PTR_EXP, i_sh);
        double *aj = env + bas(PTR_EXP, j_sh);
        double *ak = env + bas(PTR_EXP, k_sh);
        double *ci = env + bas(PTR_COEFF, i_sh);
        double *cj = env + bas(PTR_COEFF, j_sh);
        double *ck = env + bas(PTR_COEFF, k_sh);

        double expcutoff = envs->expcutoff;
        double *log_maxci, *log_maxcj;
        PairData *pdata_base, *pdata_ij;
        log_maxci = new double[i_prim+j_prim];
        pdata_base = new PairData[i_prim*j_prim];
        log_maxcj = log_maxci + i_prim;
        CINTOpt_log_max_pgto_coeff(log_maxci, ci, i_prim, i_ctr);
        CINTOpt_log_max_pgto_coeff(log_maxcj, cj, j_prim, j_ctr);
        if (CINTset_pairdata(pdata_base, ai, aj, envs->ri, envs->rj,
                             log_maxci, log_maxcj, envs->li_ceil, envs->lj_ceil,
                             i_prim, j_prim, SQUARE(envs->rirj), expcutoff)) {
                return 0;
        }

        int n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
        double fac1i, fac1j, fac1k;
        int ip, jp, kp;
        int empty[4] = {1, 1, 1, 1};
        int *iempty = empty + 0;
        int *jempty = empty + 1;
        int *kempty = empty + 2;
        int *gempty = empty + 3;
        /* COMMON_ENVS_AND_DECLARE end */

        double expij;
        double *rij;
        //const double dist_ij = SQUARE(envs->rirj);

        int *idx = new int[envs->nf * 3];
        CINTg2e_index_xyz(idx, envs);

        int *non0ctri, *non0ctrj, *non0ctrk;
        int *non0idxi, *non0idxj, *non0idxk;
        non0ctri = new int[i_prim+j_prim+k_prim+i_prim*i_ctr+j_prim*j_ctr+k_prim*k_ctr];
        non0ctrj = non0ctri + i_prim;
        non0ctrk = non0ctrj + j_prim;
        non0idxi = non0ctrk + k_prim;
        non0idxj = non0idxi + i_prim*i_ctr;
        non0idxk = non0idxj + j_prim*j_ctr;
        CINTOpt_non0coeff_byshell(non0idxi, non0ctri, ci, i_prim, i_ctr);
        CINTOpt_non0coeff_byshell(non0idxj, non0ctrj, cj, j_prim, j_ctr);
        CINTOpt_non0coeff_byshell(non0idxk, non0ctrk, ck, k_prim, k_ctr);

        const int nc = i_ctr * j_ctr * k_ctr;
        // (irys,i,j,k,l,coord,0:1); +1 for nabla-r12
        const int leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const int lenk = envs->nf * nc * n_comp; // gctrk
        const int lenj = envs->nf * i_ctr * j_ctr * n_comp; // gctrj
        const int leni = envs->nf * i_ctr * n_comp; // gctri
        const int len0 = envs->nf * n_comp; // gout
        const int len = leng + lenk + lenj + leni + len0;
        double *g = new double[len];  // must be allocated last in this function
        double *g1 = g + leng;
        double *gout, *gctri, *gctrj, *gctrk;

        if (n_comp == 1) {
                gctrk = gctr;
        } else {
                gctrk = g1;
                g1 += lenk;
        }
        if (k_ctr == 1) {
                gctrj = gctrk;
                jempty = kempty;
        } else {
                gctrj = g1;
                g1 += lenj;
        }
        if (j_ctr == 1) {
                gctri = gctrj;
                iempty = jempty;
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

        *kempty = 1;
        for (kp = 0; kp < k_prim; kp++) {
                envs->ak = ak[kp];
                envs->akl = ak[kp];
                if (k_ctr == 1) {
                        fac1k = envs->common_factor * ck[kp];
                } else {
                        fac1k = envs->common_factor;
                        *jempty = 1;
                }

                pdata_ij = pdata_base;
                for (jp = 0; jp < j_prim; jp++) {
                        envs->aj = aj[jp];
                        if (j_ctr == 1) {
                                fac1j = fac1k * cj[jp];
                        } else {
                                fac1j = fac1k;
                                *iempty = 1;
                        }
                        for (ip = 0; ip < i_prim; ip++, pdata_ij++) {
                                if (pdata_ij->cceij > expcutoff) {
                                        goto i_contracted;
                                }
                                envs->ai = ai[ip];
                                envs->aij = ai[ip] + aj[jp];
                                expij = pdata_ij->eij;
                                rij = pdata_ij->rij;
                                envs->rij[0] = rij[0];
                                envs->rij[1] = rij[1];
                                envs->rij[2] = rij[2];
                                envs->rijrx[0] = rij[0] - envs->rx_in_rijrx[0];
                                envs->rijrx[1] = rij[1] - envs->rx_in_rijrx[1];
                                envs->rijrx[2] = rij[2] - envs->rx_in_rijrx[2];
                                if (i_ctr == 1) {
                                        fac1i = fac1j*ci[ip]*expij;
                                } else {
                                        fac1i = fac1j*expij;
                                }
                                if ((*envs->f_g0_2e)(g, fac1i, envs)) {
                                        (*envs->f_gout)(gout, g, idx, envs, *gempty);
                                        PRIM2CTR0(i, gout, envs->nf*n_comp);
                                }
i_contracted: ;
                        } // end loop i_prim
                        if (!*iempty) {
                                PRIM2CTR0(j, gctri, envs->nf*i_ctr*n_comp);
                        }
                } // end loop j_prim
                if (!*jempty) {
                        PRIM2CTR0(k, gctrj, envs->nf*i_ctr*j_ctr*n_comp);
                }
        } // end loop k_prim

        if (n_comp > 1 && !*kempty) {
                CINTdmat_transpose(gctr, gctrk, envs->nf*nc, n_comp);
        }
        return !*kempty;
}


#define COMMON_ENVS_AND_DECLARE \
        int *shls = envs->shls; \
        int *bas = envs->bas; \
        double *env = envs->env; \
        int i_sh = shls[0]; \
        int j_sh = shls[1]; \
        if (opt->pairdata != NULL && \
            opt->pairdata[i_sh*opt->nbas+j_sh] == NOVALUE) { \
                return 0; \
        } \
        int k_sh = shls[2]; \
        int i_ctr = envs->x_ctr[0]; \
        int j_ctr = envs->x_ctr[1]; \
        int k_ctr = envs->x_ctr[2]; \
        int i_prim = bas(NPRIM_OF, i_sh); \
        int j_prim = bas(NPRIM_OF, j_sh); \
        int k_prim = bas(NPRIM_OF, k_sh); \
        double *ai = env + bas(PTR_EXP, i_sh); \
        double *aj = env + bas(PTR_EXP, j_sh); \
        double *ak = env + bas(PTR_EXP, k_sh); \
        double *ci = env + bas(PTR_COEFF, i_sh); \
        double *cj = env + bas(PTR_COEFF, j_sh); \
        double *ck = env + bas(PTR_COEFF, k_sh); \
        double expcutoff = envs->expcutoff; \
        PairData *pdata_base, *pdata_ij; \
        if (opt->pairdata != NULL) { \
                pdata_base = opt->pairdata[i_sh*opt->nbas+j_sh]; \
        } else { \
                double *log_maxci = opt->log_max_coeff[i_sh]; \
                double *log_maxcj = opt->log_max_coeff[j_sh]; \
                pdata_base = new PairData[i_prim*j_prim]; \
                if (CINTset_pairdata(pdata_base, ai, aj, envs->ri, envs->rj, \
                                     log_maxci, log_maxcj, envs->li_ceil, envs->lj_ceil, \
                                     i_prim, j_prim, SQUARE(envs->rirj), expcutoff)) { \
                        return 0; \
                } \
        } \
        int n_comp = envs->ncomp_e1 * envs->ncomp_tensor; \
        double fac1i, fac1j, fac1k; \
        int ip, jp, kp; \
        int empty[4] = {1, 1, 1, 1}; \
        int *iempty = empty + 0; \
        int *jempty = empty + 1; \
        int *kempty = empty + 2; \
        int *gempty = empty + 3; \
        int *non0ctri = opt->non0ctr[i_sh]; \
        int *non0ctrj = opt->non0ctr[j_sh]; \
        int *non0idxi = opt->sortedidx[i_sh]; \
        int *non0idxj = opt->sortedidx[j_sh]; \
        int *non0ctrk, *non0idxk; \
        non0ctrk = new int[k_prim+k_prim*k_ctr]; \
        non0idxk = non0ctrk + k_prim; \
        CINTOpt_non0coeff_byshell(non0idxk, non0ctrk, ck, k_prim, k_ctr); \
        double expij; \
        double *rij; \
        int *idx = opt->index_xyz_array[envs->i_l*LMAX1*LMAX1 \
                                        +envs->j_l*LMAX1+envs->k_l]; \
        if (idx == NULL) { \
                idx = new int[envs->nf * 3]; \
                CINTg2e_index_xyz(idx, envs); \
        }

#define SET_RIJ    \
        if (pdata_ij->cceij > expcutoff) { \
                goto i_contracted; \
        } \
        envs->ai  = ai[ip]; \
        envs->aij = ai[ip] + aj[jp]; \
        expij = pdata_ij->eij; \
        rij = pdata_ij->rij; \
        envs->rij[0] = rij[0]; \
        envs->rij[1] = rij[1]; \
        envs->rij[2] = rij[2]; \
        envs->rijrx[0] = rij[0] - envs->rx_in_rijrx[0]; \
        envs->rijrx[1] = rij[1] - envs->rx_in_rijrx[1]; \
        envs->rijrx[2] = rij[2] - envs->rx_in_rijrx[2]

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


// i_ctr = j_ctr = k_ctr = 1;
int CINT3c2e_111_loop(double *gctr, CINTEnvVars *envs, const CINTOpt *opt)
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

        for (kp = 0; kp < k_prim; kp++) {
                envs->ak = ak[kp];
                envs->akl = ak[kp];
                fac1k = envs->common_factor * ck[kp];

                pdata_ij = pdata_base;
                for (jp = 0; jp < j_prim; jp++) {
                        envs->aj = aj[jp];
                        fac1j = fac1k * cj[jp];
                        for (ip = 0; ip < i_prim; ip++, pdata_ij++) {
                                SET_RIJ;
                                fac1i = fac1j*ci[ip]*expij;
                                if ((*envs->f_g0_2e)(g, fac1i, envs)) {
                                        (*envs->f_gout)(gout, g, idx, envs, *empty);
                                        *empty = 0;
                                }
i_contracted: ;
                        } // end loop i_prim
                } // end loop j_prim
        } // end loop k_prim

        if (n_comp > 1 && !*empty) {
                CINTdmat_transpose(gctr, gout, envs->nf*nc, n_comp);
        }
        return !*empty;
}

// i_ctr = n; j_ctr = k_ctr = 1;
int CINT3c2e_n11_loop(double *gctr, CINTEnvVars *envs, const CINTOpt *opt)
{
        COMMON_ENVS_AND_DECLARE;

        const int nc = i_ctr;
        const int leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const int leni = envs->nf * i_ctr * n_comp; // gctri
        const int len0 = envs->nf * n_comp; // gout
        const int len = leng + leni + len0;
        double *g = new double[len];
        double *g1 = g + leng;
        double *gout, *gctri;
        if (n_comp == 1) {
                gctri = gctr;
        } else {
                gctri = g1;
                g1 += leni;
        }
        gout = g1;

        for (kp = 0; kp < k_prim; kp++) {
                envs->ak = ak[kp];
                envs->akl = ak[kp];
                fac1k = envs->common_factor * ck[kp];
                pdata_ij = pdata_base;
                for (jp = 0; jp < j_prim; jp++) {
                        envs->aj = aj[jp];
                        fac1j = fac1k * cj[jp];
                        for (ip = 0; ip < i_prim; ip++, pdata_ij++) {
                                SET_RIJ;
                                fac1i = fac1j*expij;
                                if ((*envs->f_g0_2e)(g, fac1i, envs)) {
                                        (*envs->f_gout)(gout, g, idx, envs, 1);
                                        PRIM2CTR(i, gout,envs->nf*n_comp);
                                }
i_contracted: ;
                        } // end loop i_prim
                } // end loop j_prim
        } // end loop k_prim

        if (n_comp > 1 && !*iempty) {
                CINTdmat_transpose(gctr, gctri, envs->nf*nc, n_comp);
        }
        return !*iempty;
}

// j_ctr = n; i_ctr = k_ctr = 1;
int CINT3c2e_1n1_loop(double *gctr, CINTEnvVars *envs, const CINTOpt *opt)
{
        COMMON_ENVS_AND_DECLARE;

        const int nc = j_ctr;
        const int leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const int lenj = envs->nf * j_ctr * n_comp; // gctrj
        const int len0 = envs->nf * n_comp; // gout
        const int len = leng + lenj + len0;
        double *g = new double[len];
        double *g1 = g + leng;
        double *gout, *gctrj;
        if (n_comp == 1) {
                gctrj = gctr;
        } else {
                gctrj = g1;
                g1 += lenj;
        }
        gout = g1;

        for (kp = 0; kp < k_prim; kp++) {
                envs->ak = ak[kp];
                envs->akl = ak[kp];
                fac1k = envs->common_factor * ck[kp];
                pdata_ij = pdata_base;
                for (jp = 0; jp < j_prim; jp++) {
                        envs->aj = aj[jp];
                        fac1j = fac1k;
                        *iempty = 1;
                        for (ip = 0; ip < i_prim; ip++, pdata_ij++) {
                                SET_RIJ;
                                fac1i = fac1j*ci[ip]*expij;
                                if ((*envs->f_g0_2e)(g, fac1i, envs)) {
                                        (*envs->f_gout)(gout, g, idx, envs, *iempty);
                                        *iempty = 0;
                                }
i_contracted: ;
                        } // end loop i_prim
                        if (!*iempty) {
                                PRIM2CTR(j, gout,envs->nf*n_comp);
                        }
                } // end loop j_prim
        } // end loop k_prim

        if (n_comp > 1 && !*jempty) {
                CINTdmat_transpose(gctr, gctrj, envs->nf*nc, n_comp);
        }
        return !*jempty;
}


int CINT3c2e_loop(double *gctr, CINTEnvVars *envs, const CINTOpt *opt)
{
        COMMON_ENVS_AND_DECLARE;
        const int nc = i_ctr * j_ctr * k_ctr;
        // (irys,i,j,k,coord,0:1); +1 for nabla-r12
        const int leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const int lenk = envs->nf * nc * n_comp; // gctrk
        const int lenj = envs->nf * i_ctr * j_ctr * n_comp; // gctrj
        const int leni = envs->nf * i_ctr * n_comp; // gctri
        const int len0 = envs->nf * n_comp; // gout
        const int len = leng + lenk + lenj + leni + len0;
        double *g = new double[len];
        double *g1 = g + leng;
        double *gout, *gctri, *gctrj, *gctrk;

        if (n_comp == 1) {
                gctrk = gctr;
        } else {
                gctrk = g1;
                g1 += lenk;
        }
        if (k_ctr == 1) {
                gctrj = gctrk;
                jempty = kempty;
        } else {
                gctrj = g1;
                g1 += lenj;
        }
        if (j_ctr == 1) {
                gctri = gctrj;
                iempty = jempty;
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

        *kempty = 1;
        for (kp = 0; kp < k_prim; kp++) {
                envs->ak = ak[kp];
                envs->akl = ak[kp];
                if (k_ctr == 1) {
                        fac1k = envs->common_factor * ck[kp];
                } else {
                        fac1k = envs->common_factor;
                        *jempty = 1;
                }
                pdata_ij = pdata_base;
                for (jp = 0; jp < j_prim; jp++) {
                        envs->aj = aj[jp];
                        if (j_ctr == 1) {
                                fac1j = fac1k * cj[jp];
                        } else {
                                fac1j = fac1k;
                                *iempty = 1;
                        }
                        for (ip = 0; ip < i_prim; ip++, pdata_ij++) {
                                SET_RIJ;
                                if (i_ctr == 1) {
                                        fac1i = fac1j*ci[ip]*expij;
                                } else {
                                        fac1i = fac1j*expij;
                                }
                                if ((*envs->f_g0_2e)(g, fac1i, envs)) {
                                        (*envs->f_gout)(gout, g, idx, envs, *gempty);
                                        PRIM2CTR(i, gout, envs->nf*n_comp);
                                }
i_contracted: ;
                        } // end loop i_prim
                        if (!*iempty) {
                                PRIM2CTR(j, gctri, envs->nf*i_ctr*n_comp);
                        }
                } // end loop j_prim
                if (!*jempty) {
                        PRIM2CTR0(k, gctrj, envs->nf*i_ctr*j_ctr*n_comp);
                }
        } // end loop k_prim

        if (n_comp > 1 && !*kempty) {
                CINTdmat_transpose(gctr, gctrk, envs->nf*nc, n_comp);
        }
        return !*kempty;
}

static int (*CINTf_3c2e_loop[8])(double *gctr, CINTEnvVars *envs, const CINTOpt *opt) = {
        CINT3c2e_loop,
        CINT3c2e_loop,
        CINT3c2e_loop,
        CINT3c2e_n11_loop,
        CINT3c2e_loop,
        CINT3c2e_1n1_loop,
        CINT3c2e_loop,
        CINT3c2e_111_loop,
};

#define PAIRDATA_NON0IDX_SIZE(ps) \
                int *bas = envs->bas; \
                int *shls  = envs->shls; \
                int i_prim = bas(NPRIM_OF, shls[0]); \
                int j_prim = bas(NPRIM_OF, shls[1]); \
                int k_prim = bas(NPRIM_OF, shls[2]); \
                int ps = (i_prim*j_prim * 5 \
                           + i_prim * x_ctr[0] \
                           + j_prim * x_ctr[1] \
                           + k_prim * x_ctr[2] \
                           +(i_prim+j_prim)*2 + k_prim + envs->nf*3);

int CINT3c2e_cart_drv(double *out, int *dims, CINTEnvVars *envs, CINTOpt *opt)
{
        int *x_ctr = envs->x_ctr;
        int nc = envs->nf * x_ctr[0] * x_ctr[1] * x_ctr[2];
        int n_comp = envs->ncomp_e1 * envs->ncomp_tensor;

        double *gctr = new double[nc*n_comp];

        int n;
        int has_value;

        if (opt != NULL) {
                n = ((envs->x_ctr[0]==1) << 2) + ((envs->x_ctr[1]==1) << 1) + (envs->x_ctr[2]==1);
                has_value = CINTf_3c2e_loop[n](gctr, envs, opt);
        } else {
                has_value = CINT3c2e_loop_nopt(gctr, envs);
        }

        int counts[4];
        counts[0] = envs->nfi * x_ctr[0];
        counts[1] = envs->nfj * x_ctr[1];
        counts[2] = envs->nfk * x_ctr[2];
        counts[3] = 1;
        if (dims == NULL) {
                dims = counts;
        }
        int nout = dims[0] * dims[1] * dims[2];
        if (has_value) {
                for (n = 0; n < n_comp; n++) {
                        c2s_cart_3c2e1(out+nout*n, gctr+nc*n, dims, envs);
                }
        } else {
                for (n = 0; n < n_comp; n++) {
                        c2s_dset0(out+nout*n, dims, counts);
                }
        }
        return has_value;
}
int CINT3c2e_spheric_drv(double *out, int *dims, CINTEnvVars *envs, CINTOpt *opt,
void (*f_e1_c2s)(double *bufijk, double *gctr, int *dims,
                   CINTEnvVars *envs), int is_ssc)
{
        int *x_ctr = envs->x_ctr;
        int nc = envs->nf * x_ctr[0] * x_ctr[1] * x_ctr[2];
        int n_comp = envs->ncomp_e1 * envs->ncomp_tensor;

        double *gctr = new double[nc*n_comp];

        int n;
        int has_value;

        if (opt != NULL) {
                n = ((envs->x_ctr[0]==1) << 2) + ((envs->x_ctr[1]==1) << 1) + (envs->x_ctr[2]==1);
                has_value = CINTf_3c2e_loop[n](gctr, envs, opt);
        } else {
                has_value = CINT3c2e_loop_nopt(gctr, envs);
        }

        int counts[4];
        counts[0] = (envs->i_l*2+1) * x_ctr[0];
        counts[1] = (envs->j_l*2+1) * x_ctr[1];
        if (is_ssc) {
                counts[2] = envs->nfk * x_ctr[2];
        } else {
                counts[2] = (envs->k_l*2+1) * x_ctr[2];
        }
        counts[3] = 1;
        if (dims == NULL) {
                dims = counts;
        }
        int nout = dims[0] * dims[1] * dims[2];
        if (has_value) {
                for (n = 0; n < n_comp; n++) {
                        (*f_e1_c2s)(out+nout*n, gctr+nc*n, dims, envs);
                }
        } else {
                for (n = 0; n < n_comp; n++) {
                        c2s_dset0(out+nout*n, dims, counts);
                }
        }
        return has_value;
}
int CINT3c2e_spinor_drv(std::complex<double> *out, int *dims, CINTEnvVars *envs, CINTOpt *opt,
void (*f_e1_c2s)(std::complex<double> *opijk, double *gctr, int *dims,
                  CINTEnvVars *envs), int is_ssc)
{
        int *x_ctr = envs->x_ctr;
        int counts[4];
        counts[0] = CINTcgto_spinor(envs->shls[0], envs->bas);
        counts[1] = CINTcgto_spinor(envs->shls[1], envs->bas);
        if (is_ssc) {
                counts[2] = envs->nfk * x_ctr[2];
        } else {
                counts[2] = (envs->k_l*2+1) * x_ctr[2];
        }
        counts[3] = 1;
        int nc = envs->nf * x_ctr[0] * x_ctr[1] * x_ctr[2];
        int n_comp = envs->ncomp_e1 * envs->ncomp_e2 * envs->ncomp_tensor;

        double *gctr = new double[nc*n_comp];

        int n;
        int has_value;

        if (opt != NULL) {
                n = ((envs->x_ctr[0]==1) << 2) + ((envs->x_ctr[1]==1) << 1) + (envs->x_ctr[2]==1);
                has_value = CINTf_3c2e_loop[n](gctr, envs, opt);
        } else {
                has_value = CINT3c2e_loop_nopt(gctr, envs);
        }

        if (dims == NULL) {
                dims = counts;
        }
        int nout = dims[0] * dims[1] * dims[2];
        if (has_value) {
                for (n = 0; n < envs->ncomp_e2 * envs->ncomp_tensor; n++) {
                        (*f_e1_c2s)(out+nout*n, gctr, dims, envs);
                        gctr += nc * envs->ncomp_e1;
                }
        } else {
                for (n = 0; n < envs->ncomp_e2 * envs->ncomp_tensor; n++) {
                        c2s_zset0(out+nout*n, dims, counts);
                }
        }
        return has_value;
}


int int3c2e_sph(double *out, int *dims, int *shls, int *atm, int natm,
                int *bas, int nbas, double *env, CINTOpt *opt)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        return CINT3c2e_spheric_drv(out, dims, &envs, opt, &c2s_sph_3c2e1, 0);
}
void int3c2e_optimizer(CINTOpt **opt, int *atm, int natm,
                       int *bas, int nbas, double *env)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTall_3c2e_optimizer(opt, ng, atm, natm, bas, nbas, env);
}

int int3c2e_cart(double *out, int *dims, int *shls, int *atm, int natm,
                 int *bas, int nbas, double *env, CINTOpt *opt)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        return CINT3c2e_cart_drv(out, dims, &envs, opt);
}

int int3c2e_spinor(std::complex<double> *out, int *dims, int *shls, int *atm, int natm,
                   int *bas, int nbas, double *env, CINTOpt *opt)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        return CINT3c2e_spinor_drv(out, dims, &envs, opt, &c2s_sf_3c2e1, 0);
}

int int3c2e_sph_ssc(double *out, int *dims, int *shls, int *atm, int natm,
                    int *bas, int nbas, double *env, CINTOpt *opt)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        return CINT3c2e_spheric_drv(out, dims, &envs, opt, &c2s_sph_3c2e1_ssc, 1);
}
int int3c2e_spinor_ssc(std::complex<double> *out, int *dims, int *shls, int *atm, int natm,
                       int *bas, int nbas, double *env, CINTOpt *opt)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        return CINT3c2e_spinor_drv(out, dims, &envs, opt, &c2s_sf_3c2e1_ssc, 1);
}

void CINTgout2e_int3c2e_spsp1(double *g,
double *gout, int *idx, CINTEnvVars *envs, int gout_empty);
int int3c2e_spsp1_spinor_ssc(std::complex<double> *out, int *dims, int *shls, int *atm, int natm,
                             int *bas, int nbas, double *env, CINTOpt *opt)
{
        int ng[] = {1, 1, 0, 0, 2, 4, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e_int3c2e_spsp1;
        return CINT3c2e_spinor_drv(out, dims, &envs, opt, &c2s_si_3c2e1_ssc, 1);
}
void int3c2e_ssc_optimizer(CINTOpt **opt, int *atm, int natm,
                           int *bas, int nbas, double *env)
{
        int3c2e_ssc_optimizer(opt, atm, natm, bas, nbas, env);
}



ALL_CINT(int3c2e)
ALL_CINT_FORTRAN_(int3c2e)

