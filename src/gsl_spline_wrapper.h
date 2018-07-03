#ifndef _GSLWRAPPER_HEADER
#define _GSLWRAPPER_HEADER
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>

typedef struct GSL_Spline {
  gsl_spline *spline;
  gsl_interp_accel *xacc;
  double xmin, xmax;
  int allocated;
} GSL_Spline;

GSL_Spline *Create_GSL_Spline(double *x, double *y, int nx);
void Free_GSL_Spline(GSL_Spline *splinecontainer);
double Lookup_GSL_Spline(GSL_Spline *splinecontainer, double x);

#endif
