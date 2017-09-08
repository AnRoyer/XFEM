#include <iostream>
#include <cmath>
#include "bessel.h"

extern "C" {
  void zbesj_(double*, double*, double*, int*, int*, double*,
              double*, int*, int*);
  void zbesy_(double*, double*, double*, int*, int*, double*,
              double*, int*, double*,
              double*, int*);
  void zbesh_(double*, double*, double*, int*, int*, int*, double*,
              double*, int*, int*);
}


int besselJ(double n, int num, double x, double *val)
{
  int nz = 0, ierr = 0, kode = 1;
  double xi = 0.0;
  double* ji = new double[num];
  
  zbesj_(&x, &xi, &n, &kode, &num, val, ji, &nz, &ierr) ;

  delete[] ji;
  
  return 1;
}

int BesselYn(double n, int num, double x, double *val)
{
  int nz = 0, ierr = 0, kode = 1;
  double xi = 0.0;
  double* yi = new double[num];
  double* auxyr = new double[num];
  double* auxyi = new double[num];
  
  zbesy_(&x, &xi, &n, &kode, &num, val, yi, &nz, auxyr, auxyi, &ierr);
  
  delete[] yi;
  delete[] auxyr;
  delete[] auxyi;
  
  return 1;
}
