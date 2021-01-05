#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "rys_roots.hpp"

/*
from mpmath import *
mp.dps = 25
for i in range(14):
    for j in range(14):
        print('   %.19e' % cos(pi*j*(2*i+1)/28)[-26:])
*/
static const double COS_14_14[] = {
 1.                       ,
 9.9371220989324260398e-01,
 9.7492791218182361934e-01,
 9.4388333030836757409e-01,
 9.0096886790241914600e-01,
 8.4672419922828412453e-01,
 7.8183148246802980363e-01,
 7.0710678118654757274e-01,
 6.2348980185873348336e-01,
 5.3203207651533657163e-01,
 4.3388373911755812040e-01,
 3.3027906195516709698e-01,
 2.2252093395631439288e-01,
 1.1196447610330785560e-01,
 1.                       ,
 9.4388333030836757409e-01,
 7.8183148246802980363e-01,
 5.3203207651533657163e-01,
 2.2252093395631439288e-01,
-1.1196447610330785560e-01,
-4.3388373911755812040e-01,
-7.0710678118654757274e-01,
-9.0096886790241914600e-01,
-9.9371220989324260398e-01,
-9.7492791218182361934e-01,
-8.4672419922828412453e-01,
-6.2348980185873348336e-01,
-3.3027906195516709698e-01,
 1.                       ,
 8.4672419922828412453e-01,
 4.3388373911755812040e-01,
-1.1196447610330785560e-01,
-6.2348980185873348336e-01,
-9.4388333030836757409e-01,
-9.7492791218182361934e-01,
-7.0710678118654757274e-01,
-2.2252093395631439288e-01,
 3.3027906195516709698e-01,
 7.8183148246802980363e-01,
 9.9371220989324260398e-01,
 9.0096886790241914600e-01,
 5.3203207651533657163e-01,
 1.                       ,
 7.0710678118654757274e-01,
 0.                       ,
-7.0710678118654757274e-01,
-1.                       ,
-7.0710678118654757274e-01,
 0.                       ,
 7.0710678118654757274e-01,
 1.                       ,
 7.0710678118654757274e-01,
 0.                       ,
-7.0710678118654757274e-01,
-1.                       ,
-7.0710678118654757274e-01,
 1.                       ,
 5.3203207651533657163e-01,
-4.3388373911755812040e-01,
-9.9371220989324260398e-01,
-6.2348980185873348336e-01,
 3.3027906195516709698e-01,
 9.7492791218182361934e-01,
 7.0710678118654757274e-01,
-2.2252093395631439288e-01,
-9.4388333030836757409e-01,
-7.8183148246802980363e-01,
 1.1196447610330785560e-01,
 9.0096886790241914600e-01,
 8.4672419922828412453e-01,
 1.                       ,
 3.3027906195516709698e-01,
-7.8183148246802980363e-01,
-8.4672419922828412453e-01,
 2.2252093395631439288e-01,
 9.9371220989324260398e-01,
 4.3388373911755812040e-01,
-7.0710678118654757274e-01,
-9.0096886790241914600e-01,
 1.1196447610330785560e-01,
 9.7492791218182361934e-01,
 5.3203207651533657163e-01,
-6.2348980185873348336e-01,
-9.4388333030836757409e-01,
 1.                       ,
 1.1196447610330785560e-01,
-9.7492791218182361934e-01,
-3.3027906195516709698e-01,
 9.0096886790241914600e-01,
 5.3203207651533657163e-01,
-7.8183148246802980363e-01,
-7.0710678118654757274e-01,
 6.2348980185873348336e-01,
 8.4672419922828412453e-01,
-4.3388373911755812040e-01,
-9.4388333030836757409e-01,
 2.2252093395631439288e-01,
 9.9371220989324260398e-01,
 1.                       ,
-1.1196447610330785560e-01,
-9.7492791218182361934e-01,
 3.3027906195516709698e-01,
 9.0096886790241914600e-01,
-5.3203207651533657163e-01,
-7.8183148246802980363e-01,
 7.0710678118654757274e-01,
 6.2348980185873348336e-01,
-8.4672419922828412453e-01,
-4.3388373911755812040e-01,
 9.4388333030836757409e-01,
 2.2252093395631439288e-01,
-9.9371220989324260398e-01,
 1.                       ,
-3.3027906195516709698e-01,
-7.8183148246802980363e-01,
 8.4672419922828412453e-01,
 2.2252093395631439288e-01,
-9.9371220989324260398e-01,
 4.3388373911755812040e-01,
 7.0710678118654757274e-01,
-9.0096886790241914600e-01,
-1.1196447610330785560e-01,
 9.7492791218182361934e-01,
-5.3203207651533657163e-01,
-6.2348980185873348336e-01,
 9.4388333030836757409e-01,
 1.                       ,
-5.3203207651533657163e-01,
-4.3388373911755812040e-01,
 9.9371220989324260398e-01,
-6.2348980185873348336e-01,
-3.3027906195516709698e-01,
 9.7492791218182361934e-01,
-7.0710678118654757274e-01,
-2.2252093395631439288e-01,
 9.4388333030836757409e-01,
-7.8183148246802980363e-01,
-1.1196447610330785560e-01,
 9.0096886790241914600e-01,
-8.4672419922828412453e-01,
 1.                       ,
-7.0710678118654757274e-01,
 0.                       ,
 7.0710678118654757274e-01,
-1.                       ,
 7.0710678118654757274e-01,
 0.                       ,
-7.0710678118654757274e-01,
 1.                       ,
-7.0710678118654757274e-01,
 0.                       ,
 7.0710678118654757274e-01,
-1.                       ,
 7.0710678118654757274e-01,
 1.                       ,
-8.4672419922828412453e-01,
 4.3388373911755812040e-01,
 1.1196447610330785560e-01,
-6.2348980185873348336e-01,
 9.4388333030836757409e-01,
-9.7492791218182361934e-01,
 7.0710678118654757274e-01,
-2.2252093395631439288e-01,
-3.3027906195516709698e-01,
 7.8183148246802980363e-01,
-9.9371220989324260398e-01,
 9.0096886790241914600e-01,
-5.3203207651533657163e-01,
 1.                       ,
-9.4388333030836757409e-01,
 7.8183148246802980363e-01,
-5.3203207651533657163e-01,
 2.2252093395631439288e-01,
 1.1196447610330785560e-01,
-4.3388373911755812040e-01,
 7.0710678118654757274e-01,
-9.0096886790241914600e-01,
 9.9371220989324260398e-01,
-9.7492791218182361934e-01,
 8.4672419922828412453e-01,
-6.2348980185873348336e-01,
 3.3027906195516709698e-01,
 1.                       ,
-9.9371220989324260398e-01,
 9.7492791218182361934e-01,
-9.4388333030836757409e-01,
 9.0096886790241914600e-01,
-8.4672419922828412453e-01,
 7.8183148246802980363e-01,
-7.0710678118654757274e-01,
 6.2348980185873348336e-01,
-5.3203207651533657163e-01,
 4.3388373911755812040e-01,
-3.3027906195516709698e-01,
 2.2252093395631439288e-01,
-1.1196447610330785560e-01,
};

void _CINT_clenshaw_dc(double *rr, const double *x, double u, int nroot)
{
        double d[14];
        double g[14];
        double u2 = u * 2.;
        int i, j, k;

        for (i = 0; i < nroot; ++i) {
                for (j = 0; j < 14; j++) {
                        d[j] = 0;
                        g[j] = x[13*14+j];
                }
                for (k = 11; k >= 1; k-=2) {
                        for (j = 0; j < 14; j++) {
                                d[j] = u2 * g[j] - d[j] + x[(k+1)*14+j];
                        }
                        for (j = 0; j < 14; j++){
                                g[j] = u2 * d[j] - g[j] + x[k*14+j];
                        }
                }
                for (j = 0; j < 14; j++) {
                        rr[j+14*i] = u * g[j] - d[j] + x[j] * 0.5;
                }
                x += 196;
        }
}

void _CINT_clenshaw_d1(double *rr, const double *x, double u, int nroot)
{
        int i;
        double d0, g0;
        double u2 = u * 2.;

        for (i = 0; i < nroot; i++) {
                d0 = 0;
                g0 = x[13+14*i];
                d0 = u2 * g0 - d0 + x[12+14*i];
                g0 = u2 * d0 - g0 + x[11+14*i];
                d0 = u2 * g0 - d0 + x[10+14*i];
                g0 = u2 * d0 - g0 + x[9 +14*i];
                d0 = u2 * g0 - d0 + x[8 +14*i];
                g0 = u2 * d0 - g0 + x[7 +14*i];
                d0 = u2 * g0 - d0 + x[6 +14*i];
                g0 = u2 * d0 - g0 + x[5 +14*i];
                d0 = u2 * g0 - d0 + x[4 +14*i];
                g0 = u2 * d0 - g0 + x[3 +14*i];
                d0 = u2 * g0 - d0 + x[2 +14*i];
                g0 = u2 * d0 - g0 + x[1 +14*i];
                rr[i] = u * g0 - d0 + x[14*i] * 0.5;
        }
}

void _CINT_matmul_14_14(double *imc, double *im, int nroot)
{
        double o7 = 0.14285714285714285714;
        double s;
        double d0[14];
        int i, j, k;
        for (i = 0; i < nroot; i++) {
                for (j = 0; j < 14; j++) {
                        d0[j] = 0;
                }
                for (j = 0; j < 14; j++) {
                        s = im[j+14*i];
                        for (k = 0; k < 14; ++k) {
                                d0[k] += s * COS_14_14[j*14+k];
                        }
                }
                for (j = 0; j < 14; j++) {
                        imc[j+14*i] = o7 * d0[j];
                }
        }
}
