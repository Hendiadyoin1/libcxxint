/*
 * Copyright (C) 2013-  Qiming Sun <osirpt.sun@gmail.com>
 */

#include <stdlib.h>
#include "cint_bas.hpp"
#include "g2e.hpp"
#include "optimizer.hpp"
#include "cint2e.hpp"

int int2c2e_gtg_sph(double *out, int *dims, int *shls, int *atm, int natm,
                  int *bas, int nbas, double *env, CINTOpt *opt)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2c2e_gtg_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        return CINT2c2e_spheric_drv(out, dims, &envs, opt);
}
void int2c2e_gtg_optimizer(CINTOpt **opt, int *atm, int natm,
                         int *bas, int nbas, double *env)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTall_2c2e_gtg_optimizer(opt, ng, atm, natm, bas, nbas, env);
}

#define ALL_CINT(NAME) \
int c##NAME##_sph(double *out, int *shls, int *atm, int natm, \
            int *bas, int nbas, double *env, CINTOpt *opt) { \
        return NAME##_sph(out, NULL, shls, atm, natm, bas, nbas, env, opt, NULL); \
} \
void c##NAME##_sph_optimizer(CINTOpt **opt, int *atm, int natm, \
                         int *bas, int nbas, double *env) { \
        NAME##_optimizer(opt, atm, natm, bas, nbas, env); \
}

ALL_CINT(int2c2e_gtg)

