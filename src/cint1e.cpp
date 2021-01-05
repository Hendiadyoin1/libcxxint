/*
 * Copyright (C) 2013-  Qiming Sun <osirpt.sun@gmail.com>
 *
 * basic cGTO integrals
 */

#include <stdlib.h>
#include <math.h>
#include <cstring>
#include "cint_bas.hpp"
#include "optimizer.hpp"
#include "g1e.hpp"
#include "cint1e.hpp"
#include "misc.hpp"
#include "cart2sph.hpp"
#include "c2f.hpp"
#include "rys_roots.hpp"


/*
 * 1e GTO integral basic loop for < i|j>, no 1/r
 */
int CINT1e_loop(double *gctr, CINTEnvVars *envs)
{
        int *shls  = envs->shls;
        int *bas = envs->bas;
        double *env = envs->env;
        int i_sh = shls[0];
        int j_sh = shls[1];
        int i_l = envs->i_l;
        int j_l = envs->j_l;
        int i_ctr = envs->x_ctr[0];
        int j_ctr = envs->x_ctr[1];
        int i_prim = bas(NPRIM_OF, i_sh);
        int j_prim = bas(NPRIM_OF, j_sh);
        int n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
        int nf = envs->nf;
        double *ri = envs->ri;
        double *rj = envs->rj;
        double *ai = env + bas(PTR_EXP, i_sh);
        double *aj = env + bas(PTR_EXP, j_sh);
        double *ci = env + bas(PTR_COEFF, i_sh);
        double *cj = env + bas(PTR_COEFF, j_sh);
        int ip, jp, n;
        int has_value = 0;
        int *idx = new int[nf * 3];
        double aij, dij, eij, rrij;
        double 
                *g = new double[envs->g_size * 3 * ((1<<envs->gbits)+1)],// +1 as buffer
                *gout = new double[nf * n_comp],
                *gctri = new double[nf * i_ctr * n_comp];
        CINTg1e_index_xyz(idx, envs);

        rrij = CINTsquare_dist(ri, rj);
        double fac = envs->common_factor * CINTcommon_fac_sp(i_l) * CINTcommon_fac_sp(j_l);
        double expcutoff = envs->expcutoff;

        for (jp = 0; jp < j_prim; jp++) {
                envs->aj = aj[jp];
                n = nf * i_ctr * n_comp;
                memset(gctri, 0, sizeof(double) * n);
                for (ip = 0; ip < i_prim; ip++) {
                        envs->ai = ai[ip];
                        aij = ai[ip] + aj[jp];
                        eij = (ai[ip] * aj[jp] / aij) * rrij;
                        if (eij > expcutoff)
                                continue;
                        has_value = 1;

                        dij = exp(-eij) / (aij * sqrt(aij)) * fac;
                        CINTg_ovlp(g, ai[ip], aj[jp], dij, envs);

                        memset(gout, 0, sizeof(double) * nf * n_comp);
                        (*envs->f_gout)(gout, g, idx, envs, 1);

                        n = nf * n_comp;
                        CINTprim_to_ctr(gctri, n, gout, 1, i_prim, i_ctr, ci+ip);
                }
                n = nf * i_ctr;
                CINTprim_to_ctr(gctr, n, gctri, n_comp, j_prim, j_ctr, cj+jp);
        }
        return has_value;
}

/*
 * For given charge distribution, calculate temporary parameter tau.
 * The charge parameter zeta is defined as    rho(r) = Norm * exp(-zeta*r^2)
 */
double CINTnuc_mod(double aij, int nuc_id, int *atm, double *env)
{
        double zeta;
        if (nuc_id < 0) {
                zeta = env[PTR_RINV_ZETA];
        } else if (atm(NUC_MOD_OF, nuc_id) == GAUSSIAN_NUC) {
                zeta = env[atm(PTR_ZETA, nuc_id)];
        } else {
                zeta = 0;
        }

        if (zeta > 0) {
                return sqrt(zeta / (aij + zeta));
        } else {
                return 1;
        }
}

/*
 * 1e GTO integral basic loop for < i|1/r|j>, no 1/r
 * if nuc_id >= 0: nuclear attraction, use nuclear model
 * if nuc_id <  0: 1/r potential, do not use nuclear model
 */
int CINT1e_nuc_loop(double *gctr, CINTEnvVars *envs, double fac, int nuc_id)
{
        int *shls  = envs->shls;
        int *atm = envs->atm;
        int *bas = envs->bas;
        double *env = envs->env;
        int i_sh = shls[0];
        int j_sh = shls[1];
        int i_l = envs->i_l;
        int j_l = envs->j_l;
        int i_ctr = envs->x_ctr[0];
        int j_ctr = envs->x_ctr[1];
        int i_prim = bas(NPRIM_OF, i_sh);
        int j_prim = bas(NPRIM_OF, j_sh);
        int nf = envs->nf;
        int n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
        double *ri = envs->ri;
        double *rj = envs->rj;
        double *ai = env + bas(PTR_EXP, i_sh);
        double *aj = env + bas(PTR_EXP, j_sh);
        double *ci = env + bas(PTR_COEFF, i_sh);
        double *cj = env + bas(PTR_COEFF, j_sh);
        int ip, jp, i, n;
        int has_value = 0;
        double tau;
        double *cr;
        double x, u[MXRYSROOTS], w[MXRYSROOTS];
        int *idx = new int[nf * 3];
        double rij[3], aij, dij, eij, rrij, t2;
        double 
                *g = new double[envs->g_size * 3 * ((1<<envs->gbits)+1)], // +1 as buffer
                *gout = new double[nf * n_comp],
                *gctri = new double[nf * i_ctr * n_comp];
        double expcutoff = envs->expcutoff;

        if (nuc_id < 0) {
                cr = &env[PTR_RINV_ORIG];
        } else {
                cr = &env[atm(PTR_COORD, nuc_id)];
        }

        CINTg1e_index_xyz(idx, envs);

        rrij = CINTsquare_dist(ri, rj);
        fac *= envs->common_factor * CINTcommon_fac_sp(i_l) * CINTcommon_fac_sp(j_l);

        for (jp = 0; jp < j_prim; jp++) {
                envs->aj = aj[jp];
                n = nf * i_ctr * n_comp;
                memset(gctri, 0, sizeof(double) * n);
                for (ip = 0; ip < i_prim; ip++) {
                        envs->ai = ai[ip];
                        aij = ai[ip] + aj[jp];
                        eij = (ai[ip] * aj[jp] / aij) * rrij;
                        if (eij > expcutoff)
                                continue;
                        has_value = 1;

                        rij[0] = (ai[ip] * ri[0] + aj[jp] * rj[0]) / aij;
                        rij[1] = (ai[ip] * ri[1] + aj[jp] * rj[1]) / aij;
                        rij[2] = (ai[ip] * ri[2] + aj[jp] * rj[2]) / aij;
                        tau = CINTnuc_mod(aij, nuc_id, atm, env);
                        x = aij * CINTsquare_dist(rij, cr) * tau * tau;
                        CINTrys_roots(envs->nrys_roots, x, u, w);

                        dij = exp(-eij) / aij * fac;
                        memset(gout, 0, sizeof(double) * nf * n_comp);
                        for (i = 0; i < envs->nrys_roots; i++) {
                                t2 = u[i] / (1 + u[i]) * tau * tau;
                                CINTg_nuc(g, aij, rij, cr, t2,
                                          dij * w[i] * tau, envs);

                                (*envs->f_gout)(gout, g, idx, envs, 1);
                        }

                        n = nf * n_comp;
                        CINTprim_to_ctr(gctri, n, gout, 1, i_prim, i_ctr, ci+ip);
                }
                n = nf * i_ctr;
                CINTprim_to_ctr(gctr, n, gctri, n_comp, j_prim, j_ctr, cj+jp);
        }
        ;
        return has_value;
}


/*
 * 1e integrals <i|O|j> without 1/r
 */
int CINT1e_drv(double *out, int *dims, CINTEnvVars *envs,
void (*f_c2s)(double *opij, double *gctr, int *dims,
                 CINTEnvVars *envs), int int1e_type)
{
        int *x_ctr = envs->x_ctr;
        int nc = envs->nf * x_ctr[0] * x_ctr[1];
        int n_comp = envs->ncomp_e1 * envs->ncomp_tensor;

        double *gctr = new double[nc*n_comp];

        int nout;
        int n;
        int has_value = 0;
        int *atm = envs->atm;
        double *env = envs->env;
        double charge_fac;

        memset(gctr, 0, sizeof(double) * nc*n_comp);

        memset(gctr,0,nc*n_comp);

        switch (int1e_type) {
        case INT1E_TYPE_OVLP:
                has_value = CINT1e_loop(gctr, envs);
                break;
        case INT1E_TYPE_RINV:
                has_value = CINT1e_nuc_loop(gctr, envs, 1, -1);
                break;
        default:
                for (n = 0; n < envs->natm; n++) {
                        if (atm(NUC_MOD_OF,n) == FRAC_CHARGE_NUC) {
                                charge_fac = -env[atm(PTR_FRAC_CHARGE,n)];
                        } else if (atm(CHARGE_OF,n) != 0) {
                                charge_fac = -abs(atm(CHARGE_OF,n));
                        } else {
                                charge_fac = 0;
                        }
                        if (charge_fac != 0) {
                                has_value = CINT1e_nuc_loop(gctr, envs, charge_fac, n)
                                        || has_value;
                        }
                }
        }

        int counts[4];
        if (dims == NULL) {
                dims = counts;
        }
        if (f_c2s == &c2s_sph_1e) {
                counts[0] = (envs->i_l*2+1) * x_ctr[0];
                counts[1] = (envs->j_l*2+1) * x_ctr[1];
        } else if (f_c2s == &c2s_cart_1e) {
                counts[0] = envs->nfi * x_ctr[0];
                counts[1] = envs->nfj * x_ctr[1];
        }
        nout = dims[0] * dims[1];
        for (n = 0; n < n_comp; n++) {
                (*f_c2s)(out+nout*n, gctr+nc*n, dims, envs);
        }
        return has_value;
}

int CINT1e_spinor_drv(std::complex<double> *out, int *dims, CINTEnvVars *envs,
void (*f_c2s)(std::complex<double> *opij, double *gctr, int *dims,
               CINTEnvVars *envs), int int1e_type)
{
        int *x_ctr = envs->x_ctr;
        int nc = envs->nf * x_ctr[0] * x_ctr[1] * envs->ncomp_e1;
        
        double *gctr = new double[nc*envs->ncomp_tensor];

        int nout;
        int n;
        int has_value = 0;
        int *atm = envs->atm;
        double charge_fac;

        memset(gctr, 0, sizeof(double) * nc*envs->ncomp_tensor);
        switch (int1e_type) {
        case INT1E_TYPE_OVLP:
                has_value = CINT1e_loop(gctr, envs);
                break;
        case INT1E_TYPE_RINV:
                has_value = CINT1e_nuc_loop(gctr, envs, 1, -1);
                break;
        default:
                for (n = 0; n < envs->natm; n++) {
                        if (atm(CHARGE_OF,n) != 0) {
                                charge_fac = -abs(atm(CHARGE_OF,n));
                                has_value = CINT1e_nuc_loop(gctr, envs, charge_fac, n)
                                        || has_value;
                        }
                }
        }

        int counts[4];
        if (dims == NULL) {
                dims = counts;
        }
        counts[0] = CINTcgto_spinor(envs->shls[0], envs->bas);
        counts[1] = CINTcgto_spinor(envs->shls[1], envs->bas);
        nout = dims[0] * dims[1];
        for (n = 0; n < envs->ncomp_tensor; n++) {
                (*f_c2s)(out+nout*n, gctr+nc*n, dims, envs);
        }

        return has_value;
}

void int1e_optimizer(CINTOpt **opt, int *atm, int natm,
                     int *bas, int nbas, double *env)
{
        *opt = NULL;
}

