/*
 * Copyright (C) 2013-  Qiming Sun <osirpt.sun@gmail.com>
 *
 * 4-center 1-electron integrals
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cint_bas.hpp"
#include "g1e.hpp"
#include "optimizer.hpp"
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

void CINTg4c1e_ovlp(double *g, double fac, const CINTEnvVars *envs);
void CINTg4c1e_index_xyz(int *idx, const CINTEnvVars *envs);
int CINTinit_int4c1e_EnvVars(CINTEnvVars *envs, int *ng, int *shls,
                             int *atm, int natm,
                             int *bas, int nbas, double *env);
void CINTgout3c1e(double *g, double *gout, int *idx,
                  CINTEnvVars *envs, int gout_empty);

int CINT4c1e_loop_nopt(double *gctr, CINTEnvVars *envs)
{
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
        double *ri = envs->ri;
        double *rj = envs->rj;
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

        int *idx = new int[envs->nf * 3];
        CINTg4c1e_index_xyz(idx, envs);

        int *non0ctri, *non0ctrj, *non0ctrk, *non0ctrl;
        int *non0idxi, *non0idxj, *non0idxk, *non0idxl;
        non0ctri = new int[i_prim+j_prim+k_prim+l_prim+i_prim*i_ctr+j_prim*j_ctr+k_prim*k_ctr+l_prim*l_ctr];
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

        int nc = i_ctr * j_ctr * k_ctr * l_ctr;
        // (irys,i,j,k,l,coord,0:1); +1 for nabla-r12
        int leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        int lenl = envs->nf * nc * n_comp; // gctrl
        int lenk = envs->nf * i_ctr * j_ctr * k_ctr * n_comp; // gctrk
        int lenj = envs->nf * i_ctr * j_ctr * n_comp; // gctrj
        int leni = envs->nf * i_ctr * n_comp; // gctri
        int len0 = envs->nf * n_comp; // gout
        int len = leng + lenl + lenk + lenj + leni + len0;
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

        double eij, ekl, expijkl, eijkl, a0;
        double rijrkl[3];
        double dist_ij = SQUARE(envs->rirj);
        double dist_kl = SQUARE(envs->rkrl);
        double common_fac = envs->common_factor * SQRTPI * M_PI *
                CINTcommon_fac_sp(envs->i_l) * CINTcommon_fac_sp(envs->j_l) *
                CINTcommon_fac_sp(envs->k_l) * CINTcommon_fac_sp(envs->l_l);

        *lempty = 1;
        for (lp = 0; lp < l_prim; lp++) {
                envs->al = al[lp];
                if (l_ctr == 1) {
                        fac1l = common_fac * cl[lp];
                } else {
                        fac1l = common_fac;
                        *kempty = 1;
                }
                for (kp = 0; kp < k_prim; kp++) {
                        envs->ak = ak[kp];
                        envs->akl = ak[kp] + al[lp];
                        ekl = dist_kl * ak[kp] * al[lp] / envs->akl;
                        if (ekl > EXPCUTOFF) {
                                goto k_contracted;
                        }
                        envs->rkl[0] = (ak[kp]*rk[0] + al[lp]*rl[0]) / envs->akl;
                        envs->rkl[1] = (ak[kp]*rk[1] + al[lp]*rl[1]) / envs->akl;
                        envs->rkl[2] = (ak[kp]*rk[2] + al[lp]*rl[2]) / envs->akl;
                        envs->rklrx[0] = envs->rkl[0] - envs->rx_in_rklrx[0];
                        envs->rklrx[1] = envs->rkl[1] - envs->rx_in_rklrx[1];
                        envs->rklrx[2] = envs->rkl[2] - envs->rx_in_rklrx[2];
                        if (k_ctr == 1) {
                                fac1k = fac1l * ck[kp];
                        } else {
                                fac1k = fac1l;
                                *jempty = 1;
                        }

                        for (jp = 0; jp < j_prim; jp++) {
                                envs->aj = aj[jp];
                                if (j_ctr == 1) {
                                        fac1j = fac1k * cj[jp];
                                } else {
                                        fac1j = fac1k;
                                        *iempty = 1;
                                }
                                for (ip = 0; ip < i_prim; ip++) {
                                        envs->ai = ai[ip];
                                        envs->aij = ai[ip] + aj[jp];
                                        eij = dist_ij * ai[ip] * aj[jp] / envs->aij;
                                        if (eij+ekl > EXPCUTOFF) {
                                                goto i_contracted;
                                        }
                                        envs->rij[0] = (ai[ip]*ri[0] + aj[jp]*rj[0]) / envs->aij;
                                        envs->rij[1] = (ai[ip]*ri[1] + aj[jp]*rj[1]) / envs->aij;
                                        envs->rij[2] = (ai[ip]*ri[2] + aj[jp]*rj[2]) / envs->aij;
                                        envs->rijrx[0] = envs->rij[0] - envs->rx_in_rijrx[0];
                                        envs->rijrx[1] = envs->rij[1] - envs->rx_in_rijrx[1];
                                        envs->rijrx[2] = envs->rij[2] - envs->rx_in_rijrx[2];
                                        rijrkl[0] = envs->rij[0] - envs->rkl[0];
                                        rijrkl[1] = envs->rij[1] - envs->rkl[1];
                                        rijrkl[2] = envs->rij[2] - envs->rkl[2];
                                        a0 = envs->aij * envs->akl / (envs->aij + envs->akl);
                                        eijkl = a0 * SQUARE(rijrkl);
                                        eijkl += eij + ekl;
                                        if (eijkl > EXPCUTOFF) {
                                                goto i_contracted;
                                        }
                                        expijkl = exp(-eijkl);
                                        if (i_ctr == 1) {
                                                fac1i = fac1j*ci[ip]*expijkl;
                                        } else {
                                                fac1i = fac1j*expijkl;
                                        }
                                        CINTg4c1e_ovlp(g, fac1i, envs);
                                        (*envs->f_gout)(gout, g, idx, envs, *gempty);
                                        PRIM2CTR0(i, gout, envs->nf*n_comp);
i_contracted: ;
                                } // end loop i_prim
                                if (!*iempty) {
                                        PRIM2CTR0(j, gctri, envs->nf*i_ctr*n_comp);
                                }
                        } // end loop j_prim
                        if (!*jempty) {
                                PRIM2CTR0(k, gctrj,envs->nf*i_ctr*j_ctr*n_comp);
                        }
k_contracted: ;
                } // end loop k_prim
                if (!*kempty) {
                        PRIM2CTR0(l, gctrk, envs->nf*i_ctr*j_ctr*k_ctr*n_comp);
                }
        } // end loop l_prim

        if (n_comp > 1 && !*lempty) {
                CINTdmat_transpose(gctr, gctrl, envs->nf*nc, n_comp);
        }
        return !*lempty;
}


#define PAIRDATA_NON0IDX_SIZE(ps) \
                int *bas = envs->bas; \
                int *shls  = envs->shls; \
                int i_prim = bas(NPRIM_OF, shls[0]); \
                int j_prim = bas(NPRIM_OF, shls[1]); \
                int k_prim = bas(NPRIM_OF, shls[2]); \
                int l_prim = bas(NPRIM_OF, shls[3]); \
                int ps = (i_prim * x_ctr[0] \
                           + j_prim * x_ctr[1] \
                           + k_prim * x_ctr[2] \
                           + l_prim * x_ctr[3] \
                           + envs->nf*3);

int CINT4c1e_cart_drv(double *out, int *dims, CINTEnvVars *envs, CINTOpt *opt)
{
        int *x_ctr = envs->x_ctr;
        int nc = envs->nf * x_ctr[0] * x_ctr[1] * x_ctr[2] * x_ctr[3];
        int n_comp = envs->ncomp_e1 * envs->ncomp_e2 * envs->ncomp_tensor;

        double *gctr = new double[nc*n_comp];

        int n, has_value;
        has_value = CINT4c1e_loop_nopt(gctr, envs);

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
int CINT4c1e_spheric_drv(double *out, int *dims, CINTEnvVars *envs, CINTOpt *opt)
{
        int *x_ctr = envs->x_ctr;
        int nc = envs->nf * x_ctr[0] * x_ctr[1] * x_ctr[2] * x_ctr[3];
        int n_comp = envs->ncomp_e1 * envs->ncomp_e2 * envs->ncomp_tensor;
        
        double *gctr = new double[nc*n_comp];

        int n, has_value;
        has_value = CINT4c1e_loop_nopt(gctr, envs);

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

int int4c1e_sph(double *out, int *dims, int *shls, int *atm, int natm,
                 int *bas, int nbas, double *env, CINTOpt *opt)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int4c1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout3c1e;
        return CINT4c1e_spheric_drv(out, dims, &envs, opt);
}
void int4c1e_optimizer(CINTOpt **opt, int *atm, int natm,
                       int *bas, int nbas, double *env)
{
        *opt = NULL;
}

int int4c1e_cart(double *out, int *dims, int *shls, int *atm, int natm,
                  int *bas, int nbas, double *env, CINTOpt *opt)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int4c1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout3c1e;
        return CINT4c1e_cart_drv(out, dims, &envs, opt);
}

int int4c1e_spinor(std::complex<double> *out, int *dims, int *shls, int *atm, int natm,
                   int *bas, int nbas, double *env, CINTOpt *opt)
{
        fprintf(stderr, "int4c1e_spinor not implemented\n");
        exit(1);
}

ALL_CINT(int4c1e)
ALL_CINT_FORTRAN_(int4c1e)
