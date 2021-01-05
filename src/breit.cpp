/*
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 *
 * Breit = Gaunt + gauge
 * Gaunt ~ - \sigma1\dot\sigma2/r12
 * gauge ~  1/2 \sigma1\dot\sigma2/r12 - 1/2 (\sigma1\dot r12) (\sigma2\dot r12)/r12^3
 * Breit ~ -1/2 \sigma1\dot\sigma2/r12 - 1/2 (\sigma1\dot r12) (\sigma2\dot r12)/r12^3
 */

#include <stdlib.h>
#include <complex>
using namespace std::literals;
#include "cint_bas.hpp"
#include "cart2sph.hpp"
#include "g1e.hpp"
#include "g2e.hpp"
#include "optimizer.hpp"
#include "cint1e.hpp"
#include "cint2e.hpp"
#include "misc.hpp"
#include "c2f.hpp"

#define DECLARE(X)      int X(std::complex<double> *out, int *dims, int *shls, \
                              int *atm, int natm, int *bas, int nbas, double *env, \
                              CINTOpt *opt)
#define function_t(name) int (*name)(std::complex<double> *out, int *dims, int *shls,\
                              int *atm, int natm, int *bas, int nbas, double *env,\
                              CINTOpt *opt)

#define BREIT0(X, ncomp_tensor) \
DECLARE(int2e_##X##_spinor); \
DECLARE(int2e_gauge_r1_##X##_spinor); \
DECLARE(int2e_gauge_r2_##X##_spinor); \
void int2e_breit_##X##_optimizer(CINTOpt **opt, int *atm, int natm, \
                                 int *bas, int nbas, double *env) \
{ \
        *opt = NULL; \
} \
int int2e_breit_##X##_spinor(std::complex<double> *out, int *dims, int *shls, \
                             int *atm, int natm, int *bas, int nbas, double *env, \
                             CINTOpt *opt) \
{ \
        return _int2e_breit_drv(out, dims, shls, atm, natm, bas, nbas, env, opt, \
                                ncomp_tensor, &int2e_##X##_spinor, \
                                &int2e_gauge_r1_##X##_spinor, &int2e_gauge_r2_##X##_spinor); \
} \
int cint2e_breit_##X##_spinor(std::complex<double> *out, int *shls, \
                      int *atm, int natm, int *bas, int nbas, double *env, \
                      CINTOpt *opt) \
{ \
        return int2e_breit_##X##_spinor(out, NULL, shls, atm, natm, bas, nbas, env, opt); \
}

static void _copy_to_out(std::complex<double> *out, std::complex<double> *in, int *dims, int *counts)
{
        if (out == in) {
                return;
        }
        int ni = dims[0];
        int nj = dims[1];
        int nk = dims[2];
        int nij = ni * nj;
        int nijk = nij * nk;
        int di = counts[0];
        int dj = counts[1];
        int dk = counts[2];
        int dl = counts[3];
        int dij = di * dj;
        int dijk = dij * dk;
        int i, j, k, l;
        std::complex<double> *pin, *pout;
        for (l = 0; l < dl; l++) {
                for (k = 0; k < dk; k++) {
                        pin  = in  + k * dij;
                        pout = out + k * nij;
                        for (j = 0; j < dj; j++) {
                        for (i = 0; i < di; i++) {
                                pout[j*ni+i] = pin[j*di+i];
                        } }
                }
                in  += dijk;
                out += nijk;
        }
}

static int _int2e_breit_drv(std::complex<double> *out, int *dims, int *shls,
                            int *atm, int natm, int *bas, int nbas, double *env,
                            CINTOpt *opt, int ncomp_tensor,
                            function_t(f_gaunt), function_t(f_gauge_r1), function_t(f_gauge_r2))
{
        int counts[4];
        counts[0] = CINTcgto_spinor(shls[0], bas);
        counts[1] = CINTcgto_spinor(shls[1], bas);
        counts[2] = CINTcgto_spinor(shls[2], bas);
        counts[3] = CINTcgto_spinor(shls[3], bas);
        int nop = counts[0] * counts[1] * counts[2] * counts[3] * ncomp_tensor;
        std::complex<double> *buf = new std::complex<double>[nop*2];
        std::complex<double> *buf1;
        if (dims == NULL) {
                dims = counts;
                buf1 = out;
        } else {
                buf1 = buf + nop;
        }

        int has_value = (*f_gaunt)(buf1, NULL, shls, atm, natm, bas, nbas, env, NULL);

        int i;
        has_value = ((*f_gauge_r1)(buf, NULL, shls, atm, natm, bas, nbas, env, NULL) ||
                     has_value);
        /* [1/2 gaunt] - [1/2 xxx*\sigma1\dot r1] */
        if (has_value) {
                for (i = 0; i < nop; i++) {
                        buf1[i] = -buf1[i] - buf[i];
                }
        }
        /* ... [- 1/2 xxx*\sigma1\dot(-r2)] */
        has_value = ((*f_gauge_r2)(buf, NULL, shls, atm, natm, bas, nbas, env, NULL) ||
                     has_value);
        if (has_value) {
                for (i = 0; i < nop; i++) {
                        buf1[i] = (buf1[i] + buf[i]) * .5;
                }
        }
        _copy_to_out(out, buf1, dims, counts);
        
        return has_value;
}


BREIT0(ssp1ssp2, 1);
BREIT0(ssp1sps2, 1);
BREIT0(sps1ssp2, 1);
BREIT0(sps1sps2, 1);

/* based on
 * '("int2e_breit_r1p2"  ( nabla \, r0 \| dot nabla-r12 \| \, nabla ))
 */
static void CINTgout2e_int2e_breit_r1p2(double *gout, double *g,
                int *idx, CINTEnvVars *envs, int gout_empty) {
        int nf = envs->nf;
        int nrys_roots = envs->nrys_roots;
        int ix, iy, iz, i, n;
        double *g0 = g;
        double *g1 = g0 + envs->g_size * 3;
        double *g2 = g1 + envs->g_size * 3;
        double *g3 = g2 + envs->g_size * 3;
        double *g4 = g3 + envs->g_size * 3;
        double *g5 = g4 + envs->g_size * 3;
        double *g6 = g5 + envs->g_size * 3;
        double *g7 = g6 + envs->g_size * 3;
        double *g8 = g7 + envs->g_size * 3;
        double *g9 = g8 + envs->g_size * 3;
        double *g10 = g9 + envs->g_size * 3;
        double *g11 = g10 + envs->g_size * 3;
        double *g12 = g11 + envs->g_size * 3;
        double *g13 = g12 + envs->g_size * 3;
        double *g14 = g13 + envs->g_size * 3;
        double *g15 = g14 + envs->g_size * 3;
        G2E_D_L(g1, g0, envs->i_l+2, envs->j_l+2, envs->k_l+0, envs->l_l+0);
        G2E_R0J(g3, g1, envs->i_l+1, envs->j_l+0, envs->k_l, envs->l_l);
        G2E_D_J(g4, g0, envs->i_l+1, envs->j_l+1, envs->k_l, envs->l_l);
        G2E_D_I(g5, g0, envs->i_l+1, envs->j_l+1, envs->k_l, envs->l_l);
        for (ix = 0; ix < envs->g_size * 3; ix++) {g4[ix] += g5[ix];}
        G2E_D_J(g5, g1, envs->i_l+1, envs->j_l+1, envs->k_l, envs->l_l);
        G2E_D_I(g6, g1, envs->i_l+1, envs->j_l+1, envs->k_l, envs->l_l);
        for (ix = 0; ix < envs->g_size * 3; ix++) {g5[ix] += g6[ix];}
        G2E_R0J(g7, g5, envs->i_l+1, envs->j_l+0, envs->k_l, envs->l_l);
        G2E_D_I(g12, g4, envs->i_l+0, envs->j_l, envs->k_l, envs->l_l);
        G2E_D_I(g15, g7, envs->i_l+0, envs->j_l, envs->k_l, envs->l_l);
        double s;
        for (n = 0; n < nf; n++) {
                ix = idx[0+n*3];
                iy = idx[1+n*3];
                iz = idx[2+n*3];
                s = 0;
                for (i = 0; i < nrys_roots; i++) {
                        s += g15[ix+i] * g0[iy+i] * g0[iz+i];
                        s += g12[ix+i] * g3[iy+i] * g0[iz+i];
                        s += g12[ix+i] * g0[iy+i] * g3[iz+i];
                        s += g3[ix+i] * g12[iy+i] * g0[iz+i];
                        s += g0[ix+i] * g15[iy+i] * g0[iz+i];
                        s += g0[ix+i] * g12[iy+i] * g3[iz+i];
                        s += g3[ix+i] * g0[iy+i] * g12[iz+i];
                        s += g0[ix+i] * g3[iy+i] * g12[iz+i];
                        s += g0[ix+i] * g0[iy+i] * g15[iz+i];
                }
                if (gout_empty) {
                        gout[n] = s;
                } else {
                        gout[n] += s;
                }
        }
}
void int2e_breit_r1p2_optimizer(CINTOpt **opt, int *atm, int natm, int *bas, int nbas, double *env) {
        int ng[] = {2, 2, 0, 1, 4, 1, 1, 1};
        CINTall_2e_optimizer(opt, ng, atm, natm, bas, nbas, env);
}
int int2e_breit_r1p2_cart(double *out, int *dims, int *shls,
                int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *opt) {
        int ng[] = {2, 2, 0, 1, 4, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e_int2e_breit_r1p2;
        return CINT2e_cart_drv(out, dims, &envs, opt);
} // int2e_breit_r1p2_cart
int int2e_breit_r1p2_sph(double *out, int *dims, int *shls,
                int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *opt) {
        int ng[] = {2, 2, 0, 1, 4, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e_int2e_breit_r1p2;
        return CINT2e_spheric_drv(out, dims, &envs, opt);
} // int2e_breit_r1p2_sph
int int2e_breit_r1p2_spinor(std::complex<double> *out, int *dims, int *shls,
                int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *opt) {
        int ng[] = {2, 2, 0, 1, 4, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e_int2e_breit_r1p2;
        return CINT2e_spinor_drv(out, dims, &envs, opt, &c2s_sf_2e1i, &c2s_sf_2e2i);
} // int2e_breit_r1p2_spinor
ALL_CINT(int2e_breit_r1p2)
ALL_CINT_FORTRAN_(int2e_breit_r1p2)

/* based on
 * '("int2e_breit_r2p2"  ( nabla \, r0 \| dot nabla-r12 \| \, nabla ))
 */
static void CINTgout2e_int2e_breit_r2p2(double *gout,
                double *g, int *idx, CINTEnvVars *envs, int gout_empty) {
        int nf = envs->nf;
        int nrys_roots = envs->nrys_roots;
        int ix, iy, iz, i, n;
        double *g0 = g;
        double *g1 = g0 + envs->g_size * 3;
        double *g2 = g1 + envs->g_size * 3;
        double *g3 = g2 + envs->g_size * 3;
        double *g4 = g3 + envs->g_size * 3;
        double *g5 = g4 + envs->g_size * 3;
        double *g6 = g5 + envs->g_size * 3;
        double *g7 = g6 + envs->g_size * 3;
        double *g8 = g7 + envs->g_size * 3;
        double *g9 = g8 + envs->g_size * 3;
        double *g10 = g9 + envs->g_size * 3;
        double *g11 = g10 + envs->g_size * 3;
        double *g12 = g11 + envs->g_size * 3;
        double *g13 = g12 + envs->g_size * 3;
        double *g14 = g13 + envs->g_size * 3;
        double *g15 = g14 + envs->g_size * 3;
        G2E_R0L(g2, g0, envs->i_l+2, envs->j_l+1, envs->k_l+0, envs->l_l+1);
        G2E_D_L(g3, g2, envs->i_l+2, envs->j_l+1, envs->k_l+0, envs->l_l+0);
        G2E_D_J(g4, g0, envs->i_l+1, envs->j_l+0, envs->k_l, envs->l_l);
        G2E_D_I(g5, g0, envs->i_l+1, envs->j_l+0, envs->k_l, envs->l_l);
        for (ix = 0; ix < envs->g_size * 3; ix++) {g4[ix] += g5[ix];}
        G2E_D_J(g7, g3, envs->i_l+1, envs->j_l+0, envs->k_l, envs->l_l);
        G2E_D_I(g8, g3, envs->i_l+1, envs->j_l+0, envs->k_l, envs->l_l);
        for (ix = 0; ix < envs->g_size * 3; ix++) {g7[ix] += g8[ix];}
        G2E_D_I(g12, g4, envs->i_l+0, envs->j_l, envs->k_l, envs->l_l);
        G2E_D_I(g15, g7, envs->i_l+0, envs->j_l, envs->k_l, envs->l_l);
        double s;
        for (n = 0; n < nf; n++) {
                ix = idx[0+n*3];
                iy = idx[1+n*3];
                iz = idx[2+n*3];
                s = 0;
                for (i = 0; i < nrys_roots; i++) {
                        s += g15[ix+i] * g0[iy+i] * g0[iz+i];
                        s += g12[ix+i] * g3[iy+i] * g0[iz+i];
                        s += g12[ix+i] * g0[iy+i] * g3[iz+i];
                        s += g3[ix+i] * g12[iy+i] * g0[iz+i];
                        s += g0[ix+i] * g15[iy+i] * g0[iz+i];
                        s += g0[ix+i] * g12[iy+i] * g3[iz+i];
                        s += g3[ix+i] * g0[iy+i] * g12[iz+i];
                        s += g0[ix+i] * g3[iy+i] * g12[iz+i];
                        s += g0[ix+i] * g0[iy+i] * g15[iz+i];
                }
                if (gout_empty) {
                        gout[n] = s;
                } else {
                        gout[n] += s;
                }
        }
}
void int2e_breit_r2p2_optimizer(CINTOpt **opt, int *atm, int natm, int *bas, int nbas, double *env) {
        int ng[] = {2, 1, 0, 2, 4, 1, 1, 1};
        CINTall_2e_optimizer(opt, ng, atm, natm, bas, nbas, env);
}
int int2e_breit_r2p2_cart(double *out, int *dims, int *shls,
                int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *opt) {
        int ng[] = {2, 1, 0, 2, 4, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e_int2e_breit_r2p2;
        return CINT2e_cart_drv(out, dims, &envs, opt);
} // int2e_breit_r2p2_cart
int int2e_breit_r2p2_sph(double *out, int *dims, int *shls,
                int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *opt) {
        int ng[] = {2, 1, 0, 2, 4, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e_int2e_breit_r2p2;
        return CINT2e_spheric_drv(out, dims, &envs, opt);
} // int2e_breit_r2p2_sph
int int2e_breit_r2p2_spinor(std::complex<double> *out, int *dims, int *shls,
                int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *opt) {
        int ng[] = {2, 1, 0, 2, 4, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e_int2e_breit_r2p2;
        return CINT2e_spinor_drv(out, dims, &envs, opt, &c2s_sf_2e1i, &c2s_sf_2e2i);
} // int2e_breit_r2p2_spinor
ALL_CINT(int2e_breit_r2p2)
ALL_CINT_FORTRAN_(int2e_breit_r2p2)