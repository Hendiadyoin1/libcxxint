/*
 * Copyright (C) 2013-  Qiming Sun <osirpt.sun@gmail.com>
 *
 * 3-center 1-electron integrals
 */

#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <math.h>
#include "cint_bas.hpp"
#include "g2e.hpp"
#include "g3c1e.hpp"
#include "optimizer.hpp"
#include "cint1e.hpp"
#include "cint2e.hpp"
#include "misc.hpp"
#include "cart2sph.hpp"
#include "rys_roots.hpp"
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


int CINT3c1e_loop_nopt(double *gctr, CINTEnvVars *envs)
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
        double *ri = envs->ri;
        double *rj = envs->rj;
        double *rk = envs->rk;
        double *ai = env + bas(PTR_EXP, i_sh);
        double *aj = env + bas(PTR_EXP, j_sh);
        double *ak = env + bas(PTR_EXP, k_sh);
        double *ci = env + bas(PTR_COEFF, i_sh);
        double *cj = env + bas(PTR_COEFF, j_sh);
        double *ck = env + bas(PTR_COEFF, k_sh);
        int n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
        double fac1i, fac1j, fac1k;
        int ip, jp, kp;
        int empty[4] = {1, 1, 1, 1};
        int *iempty = empty + 0;
        int *jempty = empty + 1;
        int *kempty = empty + 2;
        int *gempty = empty + 3;
        /* COMMON_ENVS_AND_DECLARE end */

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

        double eijk, dijk, aijk;
        double aiajrr, aiakrr, ajakrr;
        double rirk[3];
        double rjrk[3];
        rirk[0] = ri[0] - rk[0];
        rirk[1] = ri[1] - rk[1];
        rirk[2] = ri[2] - rk[2];
        rjrk[0] = rj[0] - rk[0];
        rjrk[1] = rj[1] - rk[1];
        rjrk[2] = rj[2] - rk[2];
        const double rr_ij = SQUARE(envs->rirj);
        const double rr_ik = SQUARE(      rirk);
        const double rr_jk = SQUARE(      rjrk);

        *kempty = 1;
        for (kp = 0; kp < k_prim; kp++) {
                envs->ak = ak[kp];
                if (k_ctr == 1) {
                        fac1k = envs->common_factor * ck[kp];
                } else {
                        fac1k = envs->common_factor;
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
                        ajakrr = aj[jp] * ak[kp] * rr_jk;
                        for (ip = 0; ip < i_prim; ip++) {
                                envs->ai = ai[ip];
                                aijk = ai[ip] + aj[jp] + ak[kp];
                                aiakrr = ai[ip] * ak[kp] * rr_ik;
                                aiajrr = ai[ip] * aj[jp] * rr_ij;
                                eijk = (aiajrr+aiakrr+ajakrr) / aijk;
                                if (eijk > EXPCUTOFF) {
                                        continue;
                                }

                                if (i_ctr == 1) {
                                        fac1i = fac1j*ci[ip]*exp(-eijk);
                                } else {
                                        fac1i = fac1j*exp(-eijk);
                                }
                                dijk = fac1i / (aijk * sqrt(aijk));
                                CINTg3c1e_ovlp(g, ai[ip], aj[jp], ak[kp], dijk, envs);
                                (*envs->f_gout)(gout, g, idx, envs, *gempty);

                                PRIM2CTR0(i, gout, envs->nf*n_comp);
                        } // end loop i_prim
                        if (!*iempty) {
                                PRIM2CTR0(j, gctri, envs->nf*i_ctr*n_comp);
                        }
                } // end loop j_prim
                if (!*jempty) {
                        PRIM2CTR0(k, gctrj,envs->nf*i_ctr*j_ctr*n_comp);
                }
        } // end loop k_prim

        if (n_comp > 1 && !*kempty) {
                CINTdmat_transpose(gctr, gctrk, envs->nf*nc, n_comp);
        }
        return !*kempty;
}

int CINT3c1e_nuc_loop_nopt(double *gctr, CINTEnvVars *envs,
                            double fac, int nuc_id)
{
        int *shls  = envs->shls;
        int *atm = envs->atm;
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
        double *ri = envs->ri;
        double *rj = envs->rj;
        double *rk = envs->rk;
        double *ai = env + bas(PTR_EXP, i_sh);
        double *aj = env + bas(PTR_EXP, j_sh);
        double *ak = env + bas(PTR_EXP, k_sh);
        double *ci = env + bas(PTR_COEFF, i_sh);
        double *cj = env + bas(PTR_COEFF, j_sh);
        double *ck = env + bas(PTR_COEFF, k_sh);
        int n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
        double fac1i, fac1j, fac1k;
        int i, ip, jp, kp;
        int empty[4] = {1, 1, 1, 1};
        int *iempty = empty + 0;
        int *jempty = empty + 1;
        int *kempty = empty + 2;
        int *gempty = empty + 3;
        int rys_empty;

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

        double *cr;
        double t2, tau, x, u[MXRYSROOTS], w[MXRYSROOTS];
        /* COMMON_ENVS_AND_DECLARE end */
        const int nc = i_ctr * j_ctr * k_ctr;
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

        if (nuc_id < 0) {
                cr = &env[PTR_RINV_ORIG];
        } else {
                cr = &env[atm(PTR_COORD, nuc_id)];
        }

        double eijk, dijk, aijk;
        double aiajrr, aiakrr, ajakrr;
        double rirk[3];
        double rjrk[3];
        double rijk[3];
        rirk[0] = ri[0] - rk[0];
        rirk[1] = ri[1] - rk[1];
        rirk[2] = ri[2] - rk[2];
        rjrk[0] = rj[0] - rk[0];
        rjrk[1] = rj[1] - rk[1];
        rjrk[2] = rj[2] - rk[2];
        const double rr_ij = SQUARE(envs->rirj);
        const double rr_ik = SQUARE(      rirk);
        const double rr_jk = SQUARE(      rjrk);
        fac *= envs->common_factor;

        *kempty = 1;
        for (kp = 0; kp < k_prim; kp++) {
                envs->ak = ak[kp];
                if (k_ctr == 1) {
                        fac1k = fac * ck[kp];
                } else {
                        fac1k = fac;
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
                        ajakrr = aj[jp] * ak[kp] * rr_jk;
                        for (ip = 0; ip < i_prim; ip++) {
                                envs->ai = ai[ip];
                                aijk = ai[ip] + aj[jp] + ak[kp];
                                aiakrr = ai[ip] * ak[kp] * rr_ik;
                                aiajrr = ai[ip] * aj[jp] * rr_ij;
                                eijk = (aiajrr+aiakrr+ajakrr) / aijk;
                                if (eijk > EXPCUTOFF) {
                                        continue;
                                }

                                if (i_ctr == 1) {
                                        fac1i = fac1j*ci[ip]*exp(-eijk);
                                } else {
                                        fac1i = fac1j*exp(-eijk);
                                }
                                dijk = fac1i / aijk;

                                rijk[0] = (ai[ip] * ri[0] + aj[jp] * rj[0] + ak[kp] * rk[0]) / aijk;
                                rijk[1] = (ai[ip] * ri[1] + aj[jp] * rj[1] + ak[kp] * rk[1]) / aijk;
                                rijk[2] = (ai[ip] * ri[2] + aj[jp] * rj[2] + ak[kp] * rk[2]) / aijk;
                                tau = CINTnuc_mod(aijk, nuc_id, atm, env);
                                x = aijk * CINTsquare_dist(rijk, cr) * tau * tau;
                                CINTrys_roots(envs->nrys_roots, x, u, w);
                                rys_empty = *gempty;
                                for (i = 0; i < envs->nrys_roots; i++) {
                                        t2 = u[i] / (1 + u[i]) * tau * tau;
                                        CINTg3c1e_nuc(g, ai[ip], aj[jp], ak[kp], rijk, cr, t2,
                                                      dijk * w[i] * tau, envs);

                                        (*envs->f_gout)(gout, g, idx, envs, rys_empty);
                                        rys_empty = 0;
                                }

                                PRIM2CTR0(i, gout, envs->nf*n_comp);
                        } // end loop i_prim
                        if (!*iempty) {
                                PRIM2CTR0(j, gctri, envs->nf*i_ctr*n_comp);
                        }
                } // end loop j_prim
                if (!*jempty) {
                        PRIM2CTR0(k, gctrj,envs->nf*i_ctr*j_ctr*n_comp);
                }
        } // end loop k_prim

        if (n_comp > 1 && !*kempty) {
                CINTdmat_transpose(gctr, gctrk, envs->nf*nc, n_comp);
        }
        return !*kempty;
}

#define PAIRDATA_NON0IDX_SIZE(ps) \
                int *bas = envs->bas; \
                int *shls  = envs->shls; \
                int i_prim = bas(NPRIM_OF, shls[0]); \
                int j_prim = bas(NPRIM_OF, shls[1]); \
                int k_prim = bas(NPRIM_OF, shls[2]); \
                int ps = (i_prim * x_ctr[0] \
                           + j_prim * x_ctr[1] \
                           + k_prim * x_ctr[2] \
                           + envs->nf*3);

int CINT3c1e_cart_drv(double *out, int *dims, CINTEnvVars *envs, CINTOpt *opt,
int int_type)
{
        int *x_ctr = envs->x_ctr;
        int nc = envs->nf * x_ctr[0] * x_ctr[1] * x_ctr[2];
        int n_comp = envs->ncomp_e1 * envs->ncomp_tensor;

        double *gctr = new double[nc*n_comp];;
        int n;
        int has_value = 0;

        if (int_type == INT1E_TYPE_OVLP) {
                has_value = CINT3c1e_loop_nopt(gctr, envs);
        } else if (int_type == INT1E_TYPE_RINV) {
                has_value = CINT3c1e_nuc_loop_nopt(gctr, envs, 1, -1);
        } else {
                int *atm = envs->atm;
                int i;
                int buf_filled;
                double fac;
                double *buf = new double[nc*n_comp];
                double *env = envs->env;
                memset(gctr,0,nc*n_comp);

                for (n = 0; n < envs->natm; n++) {
                        if (atm(NUC_MOD_OF,n) == FRAC_CHARGE_NUC) {
                                fac = -env[atm(PTR_FRAC_CHARGE,n)];
                        } else if (atm(CHARGE_OF,n) != 0) {
                                fac = -abs(atm(CHARGE_OF,n));
                        } else {
                                fac = 0;
                        }

                        if (fac != 0) {
                                buf_filled = CINT3c1e_nuc_loop_nopt(buf, envs, fac, n);
                                if (buf_filled) {
                                        for (i = 0; i < nc*n_comp; i++) {
                                                gctr[i] += buf[i];
                                        }
                                }
                                has_value |= buf_filled;
                        }
                }
                delete[] buf; buf = nullptr; // release memory for c2s function
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
int CINT3c1e_spheric_drv(double *out, int *dims, CINTEnvVars *envs, CINTOpt *opt,
                void (*f_e1_c2s)(double *out, double *gctr, int *dims,CINTEnvVars *envs),
                int int_type, int is_ssc)
{
        int *x_ctr = envs->x_ctr;
        int nc = envs->nf * x_ctr[0] * x_ctr[1] * x_ctr[2];
        int n_comp = envs->ncomp_e1 * envs->ncomp_tensor;

        double *gctr = new double[nc*n_comp];
        int n;
        int has_value = 0;

        if (int_type == INT1E_TYPE_OVLP) {
                has_value = CINT3c1e_loop_nopt(gctr, envs);
        } else if (int_type == INT1E_TYPE_RINV) {
                has_value = CINT3c1e_nuc_loop_nopt(gctr, envs, 1, -1);
        } else {
                int *atm = envs->atm;
                int i;
                int buf_filled;
                double fac;
                double *buf = new double[nc*n_comp];
                memset(gctr,0,nc*n_comp);

                for (n = 0; n < envs->natm; n++) {
                        if (atm(CHARGE_OF,n) != 0) {
                                fac = -abs(atm(CHARGE_OF,n));
                                buf_filled = CINT3c1e_nuc_loop_nopt(buf, envs, fac, n);
                                if (buf_filled) {
                                        for (i = 0; i < nc*n_comp; i++) {
                                                gctr[i] += buf[i];
                                        }
                                }
                                has_value |= buf_filled;
                        }
                }
                delete[] buf; buf = nullptr; // release memory for c2s function
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

// TODO: ssc type c2s transformation
int CINT3c1e_spinor_drv(std::complex<double> *out, int *dims, CINTEnvVars *envs, CINTOpt *opt,
void (*f_e1_c2s)(std::complex<double> *opijk, double *gctr, int *dims,
                  CINTEnvVars *envs), int int_type, int is_ssc)
{
        fprintf(stderr, "CINT3c1e_spinor_drv not implemented");
        exit(1);
}

void CINTgout3c1e(double *gout, double *g, int *idx, CINTEnvVars *envs, int gout_empty)
{
        int ix, iy, iz, n;

        if (gout_empty) {
                for (n = 0; n < envs->nf; n++, idx+=3) {
                        ix = idx[0];
                        iy = idx[1];
                        iz = idx[2];
                        gout[n] = g[ix] * g[iy] * g[iz];
                }
        } else {
                for (n = 0; n < envs->nf; n++, idx+=3) {
                        ix = idx[0];
                        iy = idx[1];
                        iz = idx[2];
                        gout[n] += g[ix] * g[iy] * g[iz];
                }
        }
}

int int3c1e_sph(double *out, int *dims, int *shls, int *atm, int natm,
                 int *bas, int nbas, double *env, CINTOpt *opt)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout3c1e;
        return CINT3c1e_spheric_drv(out, dims, &envs, opt, &c2s_sph_3c1e, 0, 0);
}
void int3c1e_optimizer(CINTOpt **opt, int *atm, int natm,
                       int *bas, int nbas, double *env)
{
        *opt = NULL;
}

int int3c1e_cart(double *out, int *dims, int *shls, int *atm, int natm,
                  int *bas, int nbas, double *env, CINTOpt *opt)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout3c1e;
        return CINT3c1e_cart_drv(out, dims, &envs, opt, 0);
}

int int3c1e_spinor(std::complex<double> *out, int *dims, int *shls, int *atm, int natm,
                   int *bas, int nbas, double *env, CINTOpt *opt)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout3c1e;
        return CINT3c1e_spinor_drv(out, dims, &envs, opt, &c2s_sf_3c2e1, 0, 0);
}

int int3c1e_rinv_sph(double *out, int *dims, int *shls, int *atm, int natm,
                      int *bas, int nbas, double *env, CINTOpt *opt)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout3c1e;
        return CINT3c1e_spheric_drv(out, dims, &envs, opt, &c2s_sph_3c1e, 1, 0);
}
void int3c1e_rinv_optimizer(CINTOpt **opt, int *atm, int natm,
                       int *bas, int nbas, double *env)
{
        *opt = NULL;
}

int int3c1e_rinv_cart(double *out, int *dims, int *shls, int *atm, int natm,
                       int *bas, int nbas, double *env, CINTOpt *opt)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout3c1e;
        return CINT3c1e_cart_drv(out, dims, &envs, opt, 1);
}

int int3c1e_rinv_spinor(std::complex<double> *out, int *dims, int *shls, int *atm, int natm,
                         int *bas, int nbas, double *env, CINTOpt *opt)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout3c1e;
        return CINT3c1e_spinor_drv(out, dims, &envs, opt, &c2s_sf_3c2e1, 1, 0);
}


ALL_CINT(int3c1e)
ALL_CINT_FORTRAN_(int3c1e)

