/*
 * Copyright (C) 2013  Qiming Sun <osirpt.sun@gmail.com>
 */
#pragma once
#include "config.hpp"

#if !defined HAVE_DEFINED_CINTOPT_H
#define HAVE_DEFINED_CINTOPT_H
struct PairData{
    double rij[3];
    double eij;
    double cceij;
};
struct CINTOpt{
    int **index_xyz_array; // LMAX1**4 pointers to index_xyz
    int **non0ctr;
    int **sortedidx;
    int nbas;
    double **log_max_coeff;
    PairData **pairdata;  // NULL indicates not-initialized, NO_VALUE can be skipped
};
#endif

#define NOVALUE                 ((void *)0xffffffffffffffffuL)
#define MAX_PGTO_FOR_PAIRDATA   2048

void CINTinit_2e_optimizer(CINTOpt **opt, int *atm, int natm,
                           int *bas, int nbas, double *env);
void CINTinit_optimizer(CINTOpt **opt, int *atm, int natm,
                        int *bas, int nbas, double *env);
void CINTdel_2e_optimizer(CINTOpt **opt);
void CINTdel_optimizer(CINTOpt **opt);
void CINTdel_pairdata_optimizer(CINTOpt *cintopt);
void CINTOpt_log_max_pgto_coeff(double *log_maxc, double *coeff, int nprim, int nctr);
void CINTOpt_set_log_maxc(CINTOpt *opt, int *atm, int natm,
                          int *bas, int nbas, double *env);
void CINTOpt_setij(CINTOpt *opt, int *ng,
                   int *atm, int natm, int *bas, int nbas, double *env);
void CINTOpt_non0coeff_byshell(int *sortedidx, int *non0ctr, double *ci,
                               int iprim, int ictr);
void CINTOpt_set_non0coeff(CINTOpt *opt, int *atm, int natm,
                           int *bas, int nbas, double *env);
int CINTset_pairdata(PairData *pairdata, double *ai, double *aj, double *ri, double *rj,
                      double *log_maxci, double *log_maxcj,
                      int li_ceil, int lj_ceil, int iprim, int jprim,
                      double rr_ij, double expcutoff);

void CINTOpt_4cindex_xyz(CINTOpt *opt, int *ng,
                         int *atm, int natm, int *bas, int nbas, double *env);
void CINTOpt_3cindex_xyz(CINTOpt *opt, int *ng,
                         int *atm, int natm, int *bas, int nbas, double *env);
void CINTOpt_2cindex_xyz(CINTOpt *opt, int *ng,
                         int *atm, int natm, int *bas, int nbas, double *env);
void CINTOpt_3c1eindex_xyz(CINTOpt *opt, int *ng,
                           int *atm, int natm, int *bas, int nbas, double *env);

// optimizer examples
void CINTno_optimizer(CINTOpt **opt, int *atm, int natm,
                      int *bas, int nbas, double *env);
void CINTall_1e_optimizer(CINTOpt **opt, int *ng,
                          int *atm, int natm, int *bas, int nbas, double *env);
void CINTall_2e_optimizer(CINTOpt **opt, int *ng,
                          int *atm, int natm, int *bas, int nbas, double *env);
void CINTall_3c2e_optimizer(CINTOpt **opt, int *ng,
                            int *atm, int natm, int *bas, int nbas, double *env);
void CINTall_2c2e_optimizer(CINTOpt **opt, int *ng,
                            int *atm, int natm, int *bas, int nbas, double *env);
void CINTall_3c1e_optimizer(CINTOpt **opt, int *ng,
                            int *atm, int natm, int *bas, int nbas, double *env);

#ifdef WITH_F12
void CINTall_2e_stg_optimizer(CINTOpt **opt, int *ng,
                              int *atm, int natm, int *bas, int nbas, double *env);
#endif

#ifdef WITH_GTG
void CINTall_2e_gtg_optimizer(CINTOpt **opt, int *ng,
                              int *atm, int natm, int *bas, int nbas, double *env);
void CINTall_3c2e_gtg_optimizer(CINTOpt **opt, int *ng,
                                int *atm, int natm, int *bas, int nbas, double *env);
void CINTall_2c2e_gtg_optimizer(CINTOpt **opt, int *ng,
                                int *atm, int natm, int *bas, int nbas, double *env);
#endif
