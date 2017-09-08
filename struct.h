#ifndef struct_h
#define struct_h

#include <complex>
#include <vector>
#include "GEntity.h"

// k = w / c
// Z = rho * c;

enum Method{
  Fem,
  Xfem,
  Lagrange
};

enum Problem2D{
  Square,
  Cylinder
};

typedef struct Param{
  double c_1;
  double c_2;
  double rho_1;
  double rho_2;
  double k_1;
  double k_2;
  double x_bnd;
  double w;
  Method method;
  Problem2D problem2D;
  std::complex<double> wave;
}Param;

typedef struct Physical{
  std::vector<GEntity*> elmDir;
  std::vector<GEntity*> elmInter;
  std::vector<GEntity*> elmInf;
  std::vector<GEntity*> elmOmega1;
  std::vector<GEntity*> elmOmega2;
  std::vector<GEntity*> elmOmega;
}Physical;

#endif /* struct_h */
