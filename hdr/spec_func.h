/*
* Fermi-Dirac function approximations are taken from:
* H. M. Antia, "Rational function approximations for Fermi-Dirac integrals",
* The Astrophysical Journal Supplement Series, 84:101-108, 1993
* The relative accuracy is up to the 1e-12 
*/

#include <gsl/gsl_sf_fermi_dirac.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_errno.h>
#include <stdio.h>
#include <math.h>

double FD_deriv_mhalf(double x);
double FD_mhalf(double x);
double FD_half(double x);
double FD_3half(double x);
double FD_mhalf_squared(double x);
// double Y(double x);
// double DY(double x);

int spec_func_test();