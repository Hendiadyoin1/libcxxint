/*
 * Copyright (C) 2013-  Qiming Sun <osirpt.sun@gmail.com>
 */

#include <stdlib.h>
#include "cint_bas.hpp"
#include "g2e.hpp"
#include "optimizer.hpp"
#include "cint2e.hpp"

int int2e_stg_sph(double *out, int *dims, int *shls, int *atm, int natm,
                  int *bas, int nbas, double *env, CINTOpt *opt)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2e_stg_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        return CINT2e_spheric_drv(out, dims, &envs, opt);
}
void int2e_stg_optimizer(CINTOpt **opt, int *atm, int natm,
                         int *bas, int nbas, double *env)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTall_2e_stg_optimizer(opt, ng, atm, natm, bas, nbas, env);
}

int int2e_yp_sph(double *out, int *dims, int *shls, int *atm, int natm,
                 int *bas, int nbas, double *env, CINTOpt *opt)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2e_yp_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        return CINT2e_spheric_drv(out, dims, &envs, opt);
}
void int2e_yp_optimizer(CINTOpt **opt, int *atm, int natm,
                        int *bas, int nbas, double *env)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTall_2e_stg_optimizer(opt, ng, atm, natm, bas, nbas, env);
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

ALL_CINT(int2e_yp)
ALL_CINT(int2e_stg)



/*
 * ((NABLA i) j| F12 |k l)
 */
void CINTgout2e_int2e_ip1(double *gout,
        double *g, int *idx, CINTEnvVars *envs, int gout_empty);

void int2e_yp_ip1_optimizer(CINTOpt **opt, int *atm, int natm,
                             int *bas, int nbas, double *env) {
        int ng[] = {1, 0, 0, 0, 1, 1, 1, 3};
        CINTall_2e_stg_optimizer(opt, ng, atm, natm, bas, nbas, env);
}
int int2e_yp_ip1_sph(double *out, int *dims, int *shls, int *atm, int natm,
                      int *bas, int nbas, double *env, CINTOpt *opt) {
        int ng[] = {1, 0, 0, 0, 1, 1, 1, 3};
        CINTEnvVars envs;
        CINTinit_int2e_yp_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e_int2e_ip1;
        return CINT2e_spheric_drv(out, dims, &envs, opt);
}
ALL_CINT(int2e_yp_ip1)

void int2e_stg_ip1_optimizer(CINTOpt **opt, int *atm, int natm,
                             int *bas, int nbas, double *env) {
        int ng[] = {1, 0, 0, 0, 1, 1, 1, 3};
        CINTall_2e_stg_optimizer(opt, ng, atm, natm, bas, nbas, env);
}
int int2e_stg_ip1_sph(double *out, int *dims, int *shls, int *atm, int natm,
                      int *bas, int nbas, double *env, CINTOpt *opt) {
        int ng[] = {1, 0, 0, 0, 1, 1, 1, 3};
        CINTEnvVars envs;
        CINTinit_int2e_stg_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e_int2e_ip1;
        return CINT2e_spheric_drv(out, dims, &envs, opt);
}
ALL_CINT(int2e_stg_ip1)

/*
 * ((NABLA NABLA i) j| F12 |k l)
 */
void CINTgout2e_int2e_ipip1(double *gout,
        double *g, int *idx, CINTEnvVars *envs, int gout_empty);

void int2e_yp_ipip1_optimizer(CINTOpt **opt, int *atm, int natm,
                              int *bas, int nbas, double *env) {
        int ng[] = {2, 0, 0, 0, 2, 1, 1, 9};
        CINTall_2e_stg_optimizer(opt, ng, atm, natm, bas, nbas, env);
}
int int2e_yp_ipip1_sph(double *out, int *dims, int *shls, int *atm, int natm,
                       int *bas, int nbas, double *env, CINTOpt *opt) {
        int ng[] = {2, 0, 0, 0, 2, 1, 1, 9};
        CINTEnvVars envs;
        CINTinit_int2e_yp_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e_int2e_ipip1;
        return CINT2e_spheric_drv(out, dims, &envs, opt);
}
ALL_CINT(int2e_yp_ipip1)

void int2e_stg_ipip1_optimizer(CINTOpt **opt, int *atm, int natm,
                               int *bas, int nbas, double *env) {
        int ng[] = {2, 0, 0, 0, 2, 1, 1, 9};
        CINTall_2e_stg_optimizer(opt, ng, atm, natm, bas, nbas, env);
}
int int2e_stg_ipip1_sph(double *out, int *dims, int *shls, int *atm, int natm,
                        int *bas, int nbas, double *env, CINTOpt *opt) {
        int ng[] = {2, 0, 0, 0, 2, 1, 1, 9};
        CINTEnvVars envs;
        CINTinit_int2e_stg_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e_int2e_ipip1;
        return CINT2e_spheric_drv(out, dims, &envs, opt);
}
ALL_CINT(int2e_stg_ipip1)

/*
 * ((NABLA i) (NABLA j)| F12 |k l)
 */
void CINTgout2e_int2e_ipvip1(double *gout,
        double *g, int *idx, CINTEnvVars *envs, int gout_empty);

void int2e_yp_ipvip1_optimizer(CINTOpt **opt, int *atm, int natm,
                               int *bas, int nbas, double *env) {
        int ng[] = {1, 1, 0, 0, 2, 1, 1, 9};
        CINTall_2e_stg_optimizer(opt, ng, atm, natm, bas, nbas, env);
}
int int2e_yp_ipvip1_sph(double *out, int *dims, int *shls, int *atm, int natm,
                        int *bas, int nbas, double *env, CINTOpt *opt) {
        int ng[] = {1, 1, 0, 0, 2, 1, 1, 9};
        CINTEnvVars envs;
        CINTinit_int2e_yp_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e_int2e_ipvip1;
        return CINT2e_spheric_drv(out, dims, &envs, opt);
}
ALL_CINT(int2e_yp_ipvip1)

void int2e_stg_ipvip1_optimizer(CINTOpt **opt, int *atm, int natm,
                                int *bas, int nbas, double *env) {
        int ng[] = {1, 1, 0, 0, 2, 1, 1, 9};
        CINTall_2e_stg_optimizer(opt, ng, atm, natm, bas, nbas, env);
}
int int2e_stg_ipvip1_sph(double *out, int *dims, int *shls, int *atm, int natm,
                         int *bas, int nbas, double *env, CINTOpt *opt) {
        int ng[] = {1, 1, 0, 0, 2, 1, 1, 9};
        CINTEnvVars envs;
        CINTinit_int2e_stg_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e_int2e_ipvip1;
        return CINT2e_spheric_drv(out, dims, &envs, opt);
}
ALL_CINT(int2e_stg_ipvip1)

/*
 * ((NABLA i) j| F12 |(NABLA k) l)
 */
void CINTgout2e_int2e_ip1ip2(double *gout,
        double *g, int *idx, CINTEnvVars *envs, int gout_empty);

void int2e_yp_ip1ip2_optimizer(CINTOpt **opt, int *atm, int natm,
                               int *bas, int nbas, double *env) {
        int ng[] = {1, 0, 1, 0, 2, 1, 1, 9};
        CINTall_2e_stg_optimizer(opt, ng, atm, natm, bas, nbas, env);
}
int int2e_yp_ip1ip2_sph(double *out, int *dims, int *shls, int *atm, int natm,
                        int *bas, int nbas, double *env, CINTOpt *opt) {
        int ng[] = {1, 0, 1, 0, 2, 1, 1, 9};
        CINTEnvVars envs;
        CINTinit_int2e_yp_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e_int2e_ip1ip2;
        return CINT2e_spheric_drv(out, dims, &envs, opt);
}
ALL_CINT(int2e_yp_ip1ip2)

void int2e_stg_ip1ip2_optimizer(CINTOpt **opt, int *atm, int natm,
                                int *bas, int nbas, double *env) {
        int ng[] = {1, 0, 1, 0, 2, 1, 1, 9};
        CINTall_2e_stg_optimizer(opt, ng, atm, natm, bas, nbas, env);
}
int int2e_stg_ip1ip2_sph(double *out, int *dims, int *shls, int *atm, int natm,
                         int *bas, int nbas, double *env, CINTOpt *opt) {
        int ng[] = {1, 0, 1, 0, 2, 1, 1, 9};
        CINTEnvVars envs;
        CINTinit_int2e_stg_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e_int2e_ip1ip2;
        return CINT2e_spheric_drv(out, dims, &envs, opt);
}
ALL_CINT(int2e_stg_ip1ip2)
