/*
 * Copyright (C) 2013-  Qiming Sun <osirpt.sun@gmail.com>
 *
 * basic cGTO integrals
 */

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

int CINT2e_loop_nopt(double *gctr, CINTEnvVars *envs)
{
        /* COMMON_ENVS_AND_DECLARE */
        int *shls  = envs->shls;
        int *bas = envs->bas;
        double *env = envs->env;
        int i_sh = shls[0];
        int j_sh = shls[1];
        int k_sh = shls[2];
        int l_sh = shls[3];
        int i_ctr = envs->x_ctr[0];
        int j_ctr = envs->x_ctr[1];
        int k_ctr = envs->x_ctr[2];
        int l_ctr = envs->x_ctr[3];
        int i_prim = bas(NPRIM_OF, i_sh);
        int j_prim = bas(NPRIM_OF, j_sh);
        int k_prim = bas(NPRIM_OF, k_sh);
        int l_prim = bas(NPRIM_OF, l_sh);
        //double *ri = envs->ri;
        //double *rj = envs->rj;
        double *rk = envs->rk;
        double *rl = envs->rl;
        double *ai = env + bas(PTR_EXP, i_sh);
        double *aj = env + bas(PTR_EXP, j_sh);
        double *ak = env + bas(PTR_EXP, k_sh);
        double *al = env + bas(PTR_EXP, l_sh);
        double *ci = env + bas(PTR_COEFF, i_sh);
        double *cj = env + bas(PTR_COEFF, j_sh);
        double *ck = env + bas(PTR_COEFF, k_sh);
        double *cl = env + bas(PTR_COEFF, l_sh);
        double expcutoff = envs->expcutoff;
        double *log_maxci, *log_maxcj, *log_maxck, *log_maxcl;
        PairData *pdata_base, *pdata_ij;
        log_maxci = new double[i_prim+j_prim+k_prim+l_prim];
        pdata_base = new PairData[i_prim*j_prim];
        log_maxcj = log_maxci + i_prim;
        log_maxck = log_maxcj + j_prim;
        log_maxcl = log_maxck + k_prim;
        CINTOpt_log_max_pgto_coeff(log_maxci, ci, i_prim, i_ctr);
        CINTOpt_log_max_pgto_coeff(log_maxcj, cj, j_prim, j_ctr);
        if (CINTset_pairdata(pdata_base, ai, aj, envs->ri, envs->rj,
                             log_maxci, log_maxcj, envs->li_ceil, envs->lj_ceil,
                             i_prim, j_prim, SQUARE(envs->rirj), expcutoff)) {
                return 0;
        }
        CINTOpt_log_max_pgto_coeff(log_maxck, ck, k_prim, k_ctr);
        CINTOpt_log_max_pgto_coeff(log_maxcl, cl, l_prim, l_ctr);

        int n_comp = envs->ncomp_e1 * envs->ncomp_e2 * envs->ncomp_tensor;
        double fac1i, fac1j, fac1k, fac1l;
        int ip, jp, kp, lp;
        int empty[5] = {1, 1, 1, 1, 1};
        int *iempty = empty + 0;
        int *jempty = empty + 1;
        int *kempty = empty + 2;
        int *lempty = empty + 3;
        int *gempty = empty + 4;
        /* COMMON_ENVS_AND_DECLARE end */

        double rr_kl = SQUARE(envs->rkrl);
        double log_rr_kl = (envs->lk_ceil+envs->ll_ceil+1)*log(rr_kl+1)/2;
        double akl, ekl, expijkl, ccekl;
        double *rij;
        //const double dist_ij = SQUARE(envs->rirj);
        const double dist_kl = SQUARE(envs->rkrl);

        int *idx = new int[envs->nf * 3];
        CINTg2e_index_xyz(idx, envs);

        int *non0ctri, *non0ctrj, *non0ctrk, *non0ctrl;
        int *non0idxi, *non0idxj, *non0idxk, *non0idxl;
        non0ctri = new int[i_prim+j_prim+k_prim+l_prim+i_prim*i_ctr+j_prim*j_ctr+k_prim*k_ctr+l_prim*l_ctr]; // single allocation is faster than 5
        non0ctrj = non0ctri + i_prim;
        non0ctrk = non0ctrj + j_prim;
        non0ctrl = non0ctrk + k_prim;
        non0idxi = non0ctrl + l_prim;
        non0idxj = non0idxi + i_prim*i_ctr;
        non0idxk = non0idxj + j_prim*j_ctr;
        non0idxl = non0idxk + k_prim*k_ctr;
        CINTOpt_non0coeff_byshell(non0idxi, non0ctri, ci, i_prim, i_ctr);
        CINTOpt_non0coeff_byshell(non0idxj, non0ctrj, cj, j_prim, j_ctr);
        CINTOpt_non0coeff_byshell(non0idxk, non0ctrk, ck, k_prim, k_ctr);
        CINTOpt_non0coeff_byshell(non0idxl, non0ctrl, cl, l_prim, l_ctr);

        const int nc = i_ctr * j_ctr * k_ctr * l_ctr;
        // (irys,i,j,k,l,coord,0:1); +1 for nabla-r12
        const int leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const int lenl = envs->nf * nc * n_comp; // gctrl
        const int lenk = envs->nf * i_ctr * j_ctr * k_ctr * n_comp; // gctrk
        const int lenj = envs->nf * i_ctr * j_ctr * n_comp; // gctrj
        const int leni = envs->nf * i_ctr * n_comp; // gctri
        const int len0 = envs->nf * n_comp; // gout
        const int len = leng + lenl + lenk + lenj + leni + len0;
        double *g = new double[len];  // must be allocated last in this function
        double *g1 = g + leng;
        double *gout, *gctri, *gctrj, *gctrk, *gctrl;
        double eijcutoff;

        if (n_comp == 1) {
                gctrl = gctr;
        } else {
                gctrl = g1;
                g1 += lenl;
        }
        if (l_ctr == 1) {
                gctrk = gctrl;
                kempty = lempty;
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

        *lempty = 1;
        for (lp = 0; lp < l_prim; lp++) {
                envs->al = al[lp];
                if (l_ctr == 1) {
                        fac1l = envs->common_factor * cl[lp];
                } else {
                        fac1l = envs->common_factor;
                        *kempty = 1;
                }
                for (kp = 0; kp < k_prim; kp++) {
                        akl = ak[kp] + al[lp];
                        ekl = dist_kl * ak[kp] * al[lp] / akl;
                        ccekl = ekl - log_rr_kl - log_maxck[kp] - log_maxcl[lp];
                        // ccekl is almost the overlap of |k> |l>. For typical
                        // chemistry systems, use overlap to prescreen eri
                        // almost 100% works.
                        // The largest error may appear for two Gaussians with
                        // dist_kl ~4^2 and ak=al ~2.5. If in the middle of |k>
                        // and |l> it happens to exist steep functions |i>, |j>.
                        // The error ~ 2*\sqrt{2*ak/pi} ~ 3. So in the worst case,
                        // an integral ~= 3*cutoff may be incorrectly dropped.
                        // Increasing expcutoff by ln(3) can guarantee to get
                        // the required accuracy in any circumstance.
                        if (ccekl > expcutoff) {
                                goto k_contracted;
                        }
                        envs->ak = ak[kp];
                        envs->akl = akl;
                        envs->rkl[0] = (ak[kp]*rk[0] + al[lp]*rl[0]) / envs->akl;
                        envs->rkl[1] = (ak[kp]*rk[1] + al[lp]*rl[1]) / envs->akl;
                        envs->rkl[2] = (ak[kp]*rk[2] + al[lp]*rl[2]) / envs->akl;
                        envs->rklrx[0] = envs->rkl[0] - envs->rx_in_rklrx[0];
                        envs->rklrx[1] = envs->rkl[1] - envs->rx_in_rklrx[1];
                        envs->rklrx[2] = envs->rkl[2] - envs->rx_in_rklrx[2];
                        eijcutoff = expcutoff - MAX(ccekl, 0);
                        ekl = exp(-ekl);

                        if (k_ctr == 1) {
                                fac1k = fac1l * ck[kp];
                        } else {
                                fac1k = fac1l;
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
                                        if (pdata_ij->cceij > eijcutoff) {
                                                goto i_contracted;
                                        }
                                        envs->ai = ai[ip];
                                        envs->aij = ai[ip] + aj[jp];
                                        rij = pdata_ij->rij;
                                        envs->rij[0] = rij[0];
                                        envs->rij[1] = rij[1];
                                        envs->rij[2] = rij[2];
                                        envs->rijrx[0] = rij[0] - envs->rx_in_rijrx[0];
                                        envs->rijrx[1] = rij[1] - envs->rx_in_rijrx[1];
                                        envs->rijrx[2] = rij[2] - envs->rx_in_rijrx[2];
                                        expijkl = pdata_ij->eij * ekl;
                                        if (i_ctr == 1) {
                                                fac1i = fac1j*ci[ip]*expijkl;
                                        } else {
                                                fac1i = fac1j*expijkl;
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
                                PRIM2CTR(k, gctrj,envs->nf*i_ctr*j_ctr*n_comp);
                        }
k_contracted: ;
                } // end loop k_prim
                if (!*kempty) {
                        PRIM2CTR(l, gctrk, envs->nf*i_ctr*j_ctr*k_ctr*n_comp);
                }
        } // end loop l_prim

        if (n_comp > 1 && !*lempty) {
                CINTdmat_transpose(gctr, gctrl, envs->nf*nc, n_comp);
        }
        return !*lempty;
}


#define COMMON_ENVS_AND_DECLARE \
        int *shls = envs->shls; \
        int *bas = envs->bas; \
        double *env = envs->env; \
        int i_sh = shls[0]; \
        int j_sh = shls[1]; \
        int k_sh = shls[2]; \
        int l_sh = shls[3]; \
        if (opt->pairdata != NULL && \
            ((opt->pairdata[i_sh*opt->nbas+j_sh] == NOVALUE) || \
             (opt->pairdata[k_sh*opt->nbas+l_sh] == NOVALUE))) { \
                return 0; \
        } \
        int i_ctr = envs->x_ctr[0]; \
        int j_ctr = envs->x_ctr[1]; \
        int k_ctr = envs->x_ctr[2]; \
        int l_ctr = envs->x_ctr[3]; \
        int i_prim = bas(NPRIM_OF, i_sh); \
        int j_prim = bas(NPRIM_OF, j_sh); \
        int k_prim = bas(NPRIM_OF, k_sh); \
        int l_prim = bas(NPRIM_OF, l_sh); \
        double *ai = env + bas(PTR_EXP, i_sh); \
        double *aj = env + bas(PTR_EXP, j_sh); \
        double *ak = env + bas(PTR_EXP, k_sh); \
        double *al = env + bas(PTR_EXP, l_sh); \
        double *ci = env + bas(PTR_COEFF, i_sh); \
        double *cj = env + bas(PTR_COEFF, j_sh); \
        double *ck = env + bas(PTR_COEFF, k_sh); \
        double *cl = env + bas(PTR_COEFF, l_sh); \
        double expcutoff = envs->expcutoff; \
        PairData *_pdata_ij, *_pdata_kl, *pdata_ij, *pdata_kl; \
        if (opt->pairdata != NULL) { \
                _pdata_ij = opt->pairdata[i_sh*opt->nbas+j_sh]; \
                _pdata_kl = opt->pairdata[k_sh*opt->nbas+l_sh]; \
        } else { \
                double *log_maxci = opt->log_max_coeff[i_sh]; \
                double *log_maxcj = opt->log_max_coeff[j_sh]; \
                _pdata_ij = new PairData[i_prim*j_prim + k_prim*l_prim]; \
                if (CINTset_pairdata(_pdata_ij, ai, aj, envs->ri, envs->rj, \
                                     log_maxci, log_maxcj, envs->li_ceil, envs->lj_ceil, \
                                     i_prim, j_prim, SQUARE(envs->rirj), expcutoff)) { \
                        return 0; \
                } \
                double *log_maxck = opt->log_max_coeff[k_sh]; \
                double *log_maxcl = opt->log_max_coeff[l_sh]; \
                _pdata_kl = _pdata_ij + i_prim*j_prim; \
                if (CINTset_pairdata(_pdata_kl, ak, al, envs->rk, envs->rl, \
                                     log_maxck, log_maxcl, envs->lk_ceil, envs->ll_ceil, \
                                     k_prim, l_prim, SQUARE(envs->rkrl), expcutoff)) { \
                        return 0; \
                } \
        } \
        int n_comp = envs->ncomp_e1 * envs->ncomp_e2 * envs->ncomp_tensor; \
        double fac1i, fac1j, fac1k, fac1l; \
        int ip, jp, kp, lp; \
        int empty[5] = {1, 1, 1, 1, 1}; \
        int *iempty = empty + 0; \
        int *jempty = empty + 1; \
        int *kempty = empty + 2; \
        int *lempty = empty + 3; \
        int *gempty = empty + 4; \
        int *non0ctri = opt->non0ctr[i_sh]; \
        int *non0ctrj = opt->non0ctr[j_sh]; \
        int *non0ctrk = opt->non0ctr[k_sh]; \
        int *non0ctrl = opt->non0ctr[l_sh]; \
        int *non0idxi = opt->sortedidx[i_sh]; \
        int *non0idxj = opt->sortedidx[j_sh]; \
        int *non0idxk = opt->sortedidx[k_sh]; \
        int *non0idxl = opt->sortedidx[l_sh]; \
        double expij, expkl, eijcutoff, eklcutoff; \
        eklcutoff = expcutoff; \
        double *rij, *rkl; \
        int *idx = opt->index_xyz_array[envs->i_l*LMAX1*LMAX1*LMAX1 \
                                        +envs->j_l*LMAX1*LMAX1 \
                                        +envs->k_l*LMAX1 \
                                        +envs->l_l]; \
        if (idx == NULL) { \
                idx = new int[envs->nf * 3]; \
                CINTg2e_index_xyz(idx, envs); \
        }

#define SET_RIJ(I,J)    \
        if (pdata_##I##J->cceij > e##I##J##cutoff) { \
                goto I##_contracted; } \
        envs->a##I = a##I[I##p]; \
        envs->a##I##J = a##I[I##p] + a##J[J##p]; \
        exp##I##J = pdata_##I##J->eij; \
        r##I##J = pdata_##I##J->rij; \
        envs->r##I##J[0] = r##I##J[0]; \
        envs->r##I##J[1] = r##I##J[1]; \
        envs->r##I##J[2] = r##I##J[2]; \
        envs->r##I##J##rx[0] = r##I##J[0] - envs->rx_in_r##I##J##rx[0]; \
        envs->r##I##J##rx[1] = r##I##J[1] - envs->rx_in_r##I##J##rx[1]; \
        envs->r##I##J##rx[2] = r##I##J[2] - envs->rx_in_r##I##J##rx[2];

// i_ctr = j_ctr = k_ctr = l_ctr = 1;
int CINT2e_1111_loop(double *gctr, CINTEnvVars *envs, CINTOpt *opt)
{
        COMMON_ENVS_AND_DECLARE;
        const int nc = 1;
        const int leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const int len0 = envs->nf * n_comp;
        const int len = leng + len0;
        double *gout;
        double *g  = new double[len];
        if (n_comp == 1) {
                gout = gctr;
        } else {
                gout = g + leng;
        }

        pdata_kl = _pdata_kl;
        for (lp = 0; lp < l_prim; lp++) {
                envs->al = al[lp];
                fac1l = envs->common_factor * cl[lp];
                for (kp = 0; kp < k_prim; kp++, pdata_kl++) {
                        SET_RIJ(k, l);
                        fac1k = fac1l * ck[kp];
                        eijcutoff = eklcutoff - MAX(pdata_kl->cceij, 0);
                        pdata_ij = _pdata_ij;
                        for (jp = 0; jp < j_prim; jp++) {
                                envs->aj = aj[jp];
                                fac1j = fac1k * cj[jp];
                                for (ip = 0; ip < i_prim; ip++, pdata_ij++) {
                                        SET_RIJ(i, j);
                                        fac1i = fac1j*ci[ip]*expij*expkl;
                                        if ((*envs->f_g0_2e)(g, fac1i, envs)) {
                                                (*envs->f_gout)(gout, g, idx, envs, *empty);
                                                *empty = 0;
                                        }
i_contracted: ;
                                } // end loop i_prim
                        } // end loop j_prim
k_contracted: ;
                } // end loop k_prim
        } // end loop l_prim

        if (n_comp > 1 && !*empty) {
                CINTdmat_transpose(gctr, gout, envs->nf*nc, n_comp);
        }
        return !*empty;
}

// i_ctr = n; j_ctr = k_ctr = l_ctr = 1;
int CINT2e_n111_loop(double *gctr, CINTEnvVars *envs, CINTOpt *opt)
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

        pdata_kl = _pdata_kl;
        for (lp = 0; lp < l_prim; lp++) {
                envs->al = al[lp];
                fac1l = envs->common_factor * cl[lp];
                for (kp = 0; kp < k_prim; kp++, pdata_kl++) {
                        SET_RIJ(k, l);
                        fac1k = fac1l * ck[kp];
                        eijcutoff = eklcutoff - MAX(pdata_kl->cceij, 0);
                        pdata_ij = _pdata_ij;
                        for (jp = 0; jp < j_prim; jp++) {
                                envs->aj = aj[jp];
                                fac1j = fac1k * cj[jp];
                                for (ip = 0; ip < i_prim; ip++, pdata_ij++) {
                                        if (pdata_ij->cceij > eijcutoff) {
                                                goto i_contracted;
                                        }
                                        SET_RIJ(i, j);
                                        fac1i = fac1j*expij*expkl;
                                        if ((*envs->f_g0_2e)(g, fac1i, envs)) {
                                                (*envs->f_gout)(gout, g, idx, envs, 1);
                                                PRIM2CTR(i, gout,envs->nf*n_comp);
                                        }
i_contracted: ;
                                } // end loop i_prim
                        } // end loop j_prim
k_contracted: ;
                } // end loop k_prim
        } // end loop l_prim

        if (n_comp > 1 && !*iempty) {
                CINTdmat_transpose(gctr, gctri, envs->nf*nc, n_comp);
        }
        return !*iempty;
}

// j_ctr = n; i_ctr = k_ctr = l_ctr = 1;
int CINT2e_1n11_loop(double *gctr, CINTEnvVars *envs, CINTOpt *opt)
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

        pdata_kl = _pdata_kl;
        for (lp = 0; lp < l_prim; lp++) {
                envs->al = al[lp];
                fac1l = envs->common_factor * cl[lp];
                for (kp = 0; kp < k_prim; kp++, pdata_kl++) {
                        SET_RIJ(k, l);
                        fac1k = fac1l * ck[kp];
                        eijcutoff = eklcutoff - MAX(pdata_kl->cceij, 0);
                        pdata_ij = _pdata_ij;
                        for (jp = 0; jp < j_prim; jp++) {
                                envs->aj = aj[jp];
                                fac1j = fac1k;
                                *iempty = 1;
                                for (ip = 0; ip < i_prim; ip++, pdata_ij++) {
                                        SET_RIJ(i, j);
                                        fac1i = fac1j*ci[ip]*expij*expkl;
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
k_contracted: ;
                } // end loop k_prim
        } // end loop l_prim

        if (n_comp > 1 && !*jempty) {
                CINTdmat_transpose(gctr, gctrj, envs->nf*nc, n_comp);
        }
        return !*jempty;
}

// k_ctr = n; i_ctr = j_ctr = l_ctr = 1;
int CINT2e_11n1_loop(double *gctr, CINTEnvVars *envs, CINTOpt *opt)
{
        COMMON_ENVS_AND_DECLARE;

        const int nc = k_ctr;
        const int leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const int lenk = envs->nf * k_ctr * n_comp; // gctrk
        const int len0 = envs->nf * n_comp; // gout
        const int len = leng + lenk + len0;
        double *g = new double[len];
        double *g1 = g + leng;
        double *gout, *gctrk;
        if (n_comp == 1) {
                gctrk = gctr;
        } else {
                gctrk = g1;
                g1 += lenk;
        }
        gout = g1;

        pdata_kl = _pdata_kl;
        for (lp = 0; lp < l_prim; lp++) {
                envs->al = al[lp];
                fac1l = envs->common_factor * cl[lp];
                for (kp = 0; kp < k_prim; kp++, pdata_kl++) {
                        SET_RIJ(k, l);
                        fac1k = fac1l;
                        eijcutoff = eklcutoff - MAX(pdata_kl->cceij, 0);
                        pdata_ij = _pdata_ij;
                        *jempty = 1;
                        for (jp = 0; jp < j_prim; jp++) {
                                envs->aj = aj[jp];
                                fac1j = fac1k * cj[jp];
                                for (ip = 0; ip < i_prim; ip++, pdata_ij++) {
                                        SET_RIJ(i, j);
                                        fac1i = fac1j*ci[ip]*expij*expkl;
                                        if ((*envs->f_g0_2e)(g, fac1i, envs)) {
                                                (*envs->f_gout)(gout, g, idx, envs, *jempty);
                                                *jempty = 0;
                                        }
i_contracted: ;
                                } // end loop i_prim
                        } // end loop j_prim
                        if (!*jempty) {
                                PRIM2CTR(k, gout,envs->nf*n_comp);
                        }
k_contracted: ;
                } // end loop k_prim
        } // end loop l_prim

        if (n_comp > 1 && !*kempty) {
                CINTdmat_transpose(gctr, gctrk, envs->nf*nc, n_comp);
        }
        return !*kempty;
}

// l_ctr = n; i_ctr = j_ctr = k_ctr = 1;
int CINT2e_111n_loop(double *gctr, CINTEnvVars *envs, CINTOpt *opt)
{
        COMMON_ENVS_AND_DECLARE;

        const int nc = l_ctr;
        const int leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const int lenl = envs->nf * l_ctr * n_comp; // gctrl
        const int len0 = envs->nf * n_comp; // gout
        const int len = leng + lenl + len0;
        double *g = new double[len];
        double *g1 = g + leng;
        double *gout, *gctrl;
        if (n_comp == 1) {
                gctrl = gctr;
        } else {
                gctrl = g1;
                g1 += lenl;
        }
        gout = g1;

        pdata_kl = _pdata_kl;
        for (lp = 0; lp < l_prim; lp++) {
                envs->al = al[lp];
                fac1l = envs->common_factor;
                *kempty = 1;
                for (kp = 0; kp < k_prim; kp++, pdata_kl++) {
                        SET_RIJ(k, l);
                        fac1k = fac1l * ck[kp];
                        eijcutoff = eklcutoff - MAX(pdata_kl->cceij, 0);
                        pdata_ij = _pdata_ij;
                        for (jp = 0; jp < j_prim; jp++) {
                                envs->aj = aj[jp];
                                fac1j = fac1k * cj[jp];
                                for (ip = 0; ip < i_prim; ip++, pdata_ij++) {
                                        SET_RIJ(i, j);
                                        fac1i = fac1j*ci[ip]*expij*expkl;
                                        if ((*envs->f_g0_2e)(g, fac1i, envs)) {
                                                (*envs->f_gout)(gout, g, idx, envs, *kempty);
                                                *kempty = 0;
                                        }
i_contracted: ;
                                } // end loop i_prim
                        } // end loop j_prim
k_contracted: ;
                } // end loop k_prim
                if (!*kempty) {
                        PRIM2CTR(l, gout,envs->nf*n_comp);
                }
        } // end loop l_prim

        if (n_comp > 1 && !*lempty) {
                CINTdmat_transpose(gctr, gctrl, envs->nf*nc, n_comp);
        }
        return !*lempty;
}


int CINT2e_loop(double *gctr, CINTEnvVars *envs, CINTOpt *opt)
{
        COMMON_ENVS_AND_DECLARE;
        const int nc = i_ctr * j_ctr * k_ctr * l_ctr;
        const int leng = envs->g_size * 3 * ((1<<envs->gbits)+1); // (irys,i,j,k,l,coord,0:1);
        const int lenl = envs->nf * nc * n_comp; // gctrl
        const int lenk = envs->nf * i_ctr * j_ctr * k_ctr * n_comp; // gctrk
        const int lenj = envs->nf * i_ctr * j_ctr * n_comp; // gctrj
        const int leni = envs->nf * i_ctr * n_comp; // gctri
        const int len0 = envs->nf * n_comp; // gout
        const int len = leng + lenl + lenk + lenj + leni + len0;
        double *g = new double[len];
        double *g1 = g + leng;
        double *gout, *gctri, *gctrj, *gctrk, *gctrl;

        if (n_comp == 1) {
                gctrl = gctr;
        } else {
                gctrl = g1;
                g1 += lenl;
        }
        if (l_ctr == 1) {
                gctrk = gctrl;
                kempty = lempty;
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

        pdata_kl = _pdata_kl;
        for (lp = 0; lp < l_prim; lp++) {
                envs->al = al[lp];
                if (l_ctr == 1) {
                        fac1l = envs->common_factor * cl[lp];
                } else {
                        fac1l = envs->common_factor;
                        *kempty = 1;
                }
                for (kp = 0; kp < k_prim; kp++, pdata_kl++) {
                        /* SET_RIJ(k, l); */
                        if (pdata_kl->cceij > eklcutoff) {
                                goto k_contracted;
                        }
                        envs->ak = ak[kp];
                        envs->akl = ak[kp] + al[lp];
                        expkl = pdata_kl->eij;
                        rkl = pdata_kl->rij;
                        envs->rkl[0] = rkl[0];
                        envs->rkl[1] = rkl[1];
                        envs->rkl[2] = rkl[2];
                        envs->rklrx[0] = rkl[0] - envs->rx_in_rklrx[0];
                        envs->rklrx[1] = rkl[1] - envs->rx_in_rklrx[1];
                        envs->rklrx[2] = rkl[2] - envs->rx_in_rklrx[2];
                        eijcutoff = eklcutoff - MAX(pdata_kl->cceij, 0);
                        /* SET_RIJ(k, l); end */
                        if (k_ctr == 1) {
                                fac1k = fac1l * ck[kp];
                        } else {
                                fac1k = fac1l;
                                *jempty = 1;
                        }

                        pdata_ij = _pdata_ij;
                        for (jp = 0; jp < j_prim; jp++) {
                                envs->aj = aj[jp];
                                if (j_ctr == 1) {
                                        fac1j = fac1k * cj[jp];
                                } else {
                                        fac1j = fac1k;
                                        *iempty = 1;
                                }
                                for (ip = 0; ip < i_prim; ip++, pdata_ij++) {
                                        /* SET_RIJ(i, j); */
                                        if (pdata_ij->cceij > eijcutoff) {
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
                                        /* SET_RIJ(i, j); end */
                                        if (i_ctr == 1) {
                                                fac1i = fac1j*ci[ip] * expij*expkl;
                                        } else {
                                                fac1i = fac1j * expij*expkl;
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
                                PRIM2CTR(k, gctrj, envs->nf*i_ctr*j_ctr*n_comp);
                        }
k_contracted: ;
                } // end loop k_prim
                if (!*kempty) {
//TODO: merge this contraction with COPY_AND_CLOSING for n_comp>1
                        PRIM2CTR(l, gctrk, envs->nf*i_ctr*j_ctr*k_ctr*n_comp);
                }
        } // end loop l_prim

        if (n_comp > 1 && !*lempty) {
                CINTdmat_transpose(gctr, gctrl, envs->nf*nc, n_comp);
        }
        return !*lempty;
}


static int (*CINTf_2e_loop[16])(double *gctr, CINTEnvVars *envs, CINTOpt *opt) = {
        CINT2e_loop,
        CINT2e_loop,
        CINT2e_loop,
        CINT2e_loop,
        CINT2e_loop,
        CINT2e_loop,
        CINT2e_loop,
        CINT2e_n111_loop,
        CINT2e_loop,
        CINT2e_loop,
        CINT2e_loop,
        CINT2e_1n11_loop,
        CINT2e_loop,
        CINT2e_11n1_loop,
        CINT2e_111n_loop,
        CINT2e_1111_loop,
};

#define PAIRDATA_NON0IDX_SIZE(ps) \
                int *bas = envs->bas; \
                int *shls  = envs->shls; \
                int i_prim = bas(NPRIM_OF, shls[0]); \
                int j_prim = bas(NPRIM_OF, shls[1]); \
                int k_prim = bas(NPRIM_OF, shls[2]); \
                int l_prim = bas(NPRIM_OF, shls[3]); \
                int ps = ((i_prim*j_prim + k_prim*l_prim) * 5 \
                           + i_prim * x_ctr[0] \
                           + j_prim * x_ctr[1] \
                           + k_prim * x_ctr[2] \
                           + l_prim * x_ctr[3] \
                           +(i_prim+j_prim+k_prim+l_prim)*2 + envs->nf*3);

int CINT2e_cart_drv(double *out, int *dims, CINTEnvVars *envs, CINTOpt *opt)
{
        int *x_ctr = envs->x_ctr;
        int nc = envs->nf * x_ctr[0] * x_ctr[1] * x_ctr[2] * x_ctr[3];
        int n_comp = envs->ncomp_e1 * envs->ncomp_e2 * envs->ncomp_tensor;

        double *gctr = new double[nc*n_comp];

        int n, has_value;
        if (opt != NULL) {
                n = ((x_ctr[0]==1) << 3) + ((x_ctr[1]==1) << 2)
                  + ((x_ctr[2]==1) << 1) +  (x_ctr[3]==1);
                has_value = CINTf_2e_loop[n](gctr, envs, opt);
        } else {
                has_value = CINT2e_loop_nopt(gctr, envs);
        }

        int counts[4];
        counts[0] = envs->nfi * x_ctr[0];
        counts[1] = envs->nfj * x_ctr[1];
        counts[2] = envs->nfk * x_ctr[2];
        counts[3] = envs->nfl * x_ctr[3];
        if (dims == NULL) {
                dims = counts;
        }
        int nout = dims[0] * dims[1] * dims[2] * dims[3];
        if (has_value) {
                for (n = 0; n < n_comp; n++) {
                        c2s_cart_2e1(out+nout*n, gctr+nc*n, dims, envs);
                }
        } else {
                for (n = 0; n < n_comp; n++) {
                        c2s_dset0(out+nout*n, dims, counts);
                }
        }
        return has_value;
}
int CINT2e_spheric_drv(double *out, int *dims, CINTEnvVars *envs, CINTOpt *opt)
{
        int *x_ctr = envs->x_ctr;
        int nc = envs->nf * x_ctr[0] * x_ctr[1] * x_ctr[2] * x_ctr[3];
        int n_comp = envs->ncomp_e1 * envs->ncomp_e2 * envs->ncomp_tensor;

        double *gctr = new double[nc*n_comp];

        int n, has_value;
        if (opt != NULL) {
                n = ((x_ctr[0]==1) << 3) + ((x_ctr[1]==1) << 2)
                  + ((x_ctr[2]==1) << 1) +  (x_ctr[3]==1);
                has_value = CINTf_2e_loop[n](gctr, envs, opt);
        } else {
                has_value = CINT2e_loop_nopt(gctr, envs);
        }

        int counts[4];
        counts[0] = (envs->i_l*2+1) * x_ctr[0];
        counts[1] = (envs->j_l*2+1) * x_ctr[1];
        counts[2] = (envs->k_l*2+1) * x_ctr[2];
        counts[3] = (envs->l_l*2+1) * x_ctr[3];
        if (dims == NULL) {
                dims = counts;
        }
        int nout = dims[0] * dims[1] * dims[2] * dims[3];
        if (has_value) {
                for (n = 0; n < n_comp; n++) {
                        c2s_sph_2e1(out+nout*n, gctr+nc*n, dims, envs);
                }
        } else {
                for (n = 0; n < n_comp; n++) {
                        c2s_dset0(out+nout*n, dims, counts);
                }
        }
        return has_value;
}
int CINT2e_spinor_drv(std::complex<double> *out, int *dims, CINTEnvVars *envs, CINTOpt *opt,
void (*f_e1_c2s)(std::complex<double> *opij, double *gctr, int *dims,
                 CINTEnvVars *envs), void (*f_e2_c2s)(std::complex<double> *fijkl, std::complex<double> *opij, int *dims,
                 CINTEnvVars *envs))
{
        int *shls = envs->shls;
        int *bas = envs->bas;
        int counts[4];
        counts[0] = CINTcgto_spinor(shls[0], bas);
        counts[1] = CINTcgto_spinor(shls[1], bas);
        counts[2] = CINTcgto_spinor(shls[2], bas);
        counts[3] = CINTcgto_spinor(shls[3], bas);
        int *x_ctr = envs->x_ctr;
        int nc = envs->nf * x_ctr[0] * x_ctr[1] * x_ctr[2] * x_ctr[3];
        int n_comp = envs->ncomp_e1 * envs->ncomp_e2 * envs->ncomp_tensor;
        int n1 = counts[0] * envs->nfk * x_ctr[2]
                           * envs->nfl * x_ctr[3] * counts[1];

        double *gctr = new double[nc*n_comp];

        int n, m, has_value;
        if (opt != NULL) {
                n = ((x_ctr[0]==1) << 3) + ((x_ctr[1]==1) << 2)
                  + ((x_ctr[2]==1) << 1) +  (x_ctr[3]==1);
                has_value = CINTf_2e_loop[n](gctr, envs, opt);
        } else {
                has_value = CINT2e_loop_nopt(gctr, envs);
        }

        if (dims == NULL) {
                dims = counts;
        }
        int nout = dims[0] * dims[1] * dims[2] * dims[3];
        if (has_value) {
                std::complex<double> *opij = new std::complex<double>[n1*envs->ncomp_e2];
                for (n = 0; n < envs->ncomp_tensor; n++) {
                        for (m = 0; m < envs->ncomp_e2; m++) {
                                (*f_e1_c2s)(opij+n1*m, gctr, dims, envs);
                                gctr += nc * envs->ncomp_e1;
                        }
                        (*f_e2_c2s)(out+nout*n, opij, dims, envs);
                }
        } else {
                for (n = 0; n < envs->ncomp_tensor; n++) {
                        c2s_zset0(out+nout*n, dims, counts);
                }
        }

        return has_value;
}


/*
 * <ki|jl> = (ij|kl); i,j\in electron 1; k,l\in electron 2
 */
void CINTgout2e(double *gout, double *g, int *idx,
                CINTEnvVars *envs, int gout_empty)
{
        int nf = envs->nf;
        int i, ix, iy, iz, n;
        double s;

        if (gout_empty) {
                switch (envs->nrys_roots) {
                        case 1:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] = g[ix] * g[iy] * g[iz];
                                }
                                break;
                        case 2:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] = g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1];
                                }
                                break;
                        case 3:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] = g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1]
                                                + g[ix+2] * g[iy+2] * g[iz+2];
                                }
                                break;
                        case 4:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] = g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1]
                                                + g[ix+2] * g[iy+2] * g[iz+2]
                                                + g[ix+3] * g[iy+3] * g[iz+3];
                                }
                                break;
                        case 5:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] = g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1]
                                                + g[ix+2] * g[iy+2] * g[iz+2]
                                                + g[ix+3] * g[iy+3] * g[iz+3]
                                                + g[ix+4] * g[iy+4] * g[iz+4];
                                }
                                break;
                        case 6:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] = g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1]
                                                + g[ix+2] * g[iy+2] * g[iz+2]
                                                + g[ix+3] * g[iy+3] * g[iz+3]
                                                + g[ix+4] * g[iy+4] * g[iz+4]
                                                + g[ix+5] * g[iy+5] * g[iz+5];
                                }
                                break;
                        case 7:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] = g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1]
                                                + g[ix+2] * g[iy+2] * g[iz+2]
                                                + g[ix+3] * g[iy+3] * g[iz+3]
                                                + g[ix+4] * g[iy+4] * g[iz+4]
                                                + g[ix+5] * g[iy+5] * g[iz+5]
                                                + g[ix+6] * g[iy+6] * g[iz+6];
                                }
                                break;
                        case 8:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] = g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1]
                                                + g[ix+2] * g[iy+2] * g[iz+2]
                                                + g[ix+3] * g[iy+3] * g[iz+3]
                                                + g[ix+4] * g[iy+4] * g[iz+4]
                                                + g[ix+5] * g[iy+5] * g[iz+5]
                                                + g[ix+6] * g[iy+6] * g[iz+6]
                                                + g[ix+7] * g[iy+7] * g[iz+7];
                                }
                                break;
                        default:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        s = 0;
                                        for (i = 0; i < envs->nrys_roots; i++) {
                                                s += g[ix+i] * g[iy+i] * g[iz+i];
                                        }
                                        gout[n] = s;
                                }
                                break;
                } // end switch nroots
        } else { // not flag_acc
                switch (envs->nrys_roots) {
                        case 1:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] += g[ix] * g[iy] * g[iz];
                                }
                                break;
                        case 2:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] +=g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1];
                                }
                                break;
                        case 3:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] +=g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1]
                                                + g[ix+2] * g[iy+2] * g[iz+2];
                                }
                                break;
                        case 4:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] +=g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1]
                                                + g[ix+2] * g[iy+2] * g[iz+2]
                                                + g[ix+3] * g[iy+3] * g[iz+3];
                                }
                                break;
                        case 5:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] +=g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1]
                                                + g[ix+2] * g[iy+2] * g[iz+2]
                                                + g[ix+3] * g[iy+3] * g[iz+3]
                                                + g[ix+4] * g[iy+4] * g[iz+4];
                                }
                                break;
                        case 6:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] +=g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1]
                                                + g[ix+2] * g[iy+2] * g[iz+2]
                                                + g[ix+3] * g[iy+3] * g[iz+3]
                                                + g[ix+4] * g[iy+4] * g[iz+4]
                                                + g[ix+5] * g[iy+5] * g[iz+5];
                                }
                                break;
                        case 7:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] +=g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1]
                                                + g[ix+2] * g[iy+2] * g[iz+2]
                                                + g[ix+3] * g[iy+3] * g[iz+3]
                                                + g[ix+4] * g[iy+4] * g[iz+4]
                                                + g[ix+5] * g[iy+5] * g[iz+5]
                                                + g[ix+6] * g[iy+6] * g[iz+6];
                                }
                                break;
                        case 8:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] +=g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1]
                                                + g[ix+2] * g[iy+2] * g[iz+2]
                                                + g[ix+3] * g[iy+3] * g[iz+3]
                                                + g[ix+4] * g[iy+4] * g[iz+4]
                                                + g[ix+5] * g[iy+5] * g[iz+5]
                                                + g[ix+6] * g[iy+6] * g[iz+6]
                                                + g[ix+7] * g[iy+7] * g[iz+7];
                                }
                                break;
                        default:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        s = 0;
                                        for (i = 0; i < envs->nrys_roots; i++) {
                                                s += g[ix+i] * g[iy+i] * g[iz+i];
                                        }
                                        gout[n] += s;
                                }
                                break;
                } // end switch nroots
        }
}

int int2e_sph(double *out, int *dims, int *shls, int *atm, int natm,
              int *bas, int nbas, double *env, CINTOpt *opt)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        return CINT2e_spheric_drv(out, dims, &envs, opt);
}
void int2e_optimizer(CINTOpt **opt, int *atm, int natm,
                     int *bas, int nbas, double *env)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTall_2e_optimizer(opt, ng, atm, natm, bas, nbas, env);
}

int int2e_cart(double *out, int *dims, int *shls, int *atm, int natm,
               int *bas, int nbas, double *env, CINTOpt *opt)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        return CINT2e_cart_drv(out, dims, &envs, opt);
}

/*
 * spinor <ki|jl> = (ij|kl); i,j\in electron 1; k,l\in electron 2
 */
int int2e_spinor(std::complex<double> *out, int *dims, int *shls, int *atm, int natm,
                 int *bas, int nbas, double *env, CINTOpt *opt)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        return CINT2e_spinor_drv(out, dims, &envs, opt,
                                 &c2s_sf_2e1, &c2s_sf_2e2);
}


ALL_CINT(int2e)
ALL_CINT_FORTRAN_(int2e)

