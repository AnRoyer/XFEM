#ifndef analytical_h
#define analytical_h

#include <vector>
#include <complex>

#include "GModel.h"

#include "struct.h"

std::vector< std::complex<double> > solveANALYTICAL(GModel* m, int nbNodes, Param param, Physical physical);

#endif /* analytical_h */
