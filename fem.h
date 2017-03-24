#ifndef fem_h
#define fem_h

#include <vector>
#include <complex>

#include "GModel.h"

#include "GEdge.h"

#include "MLine.h"

#include "gmm.h"

#include "struct.h"

std::vector< std::complex<double> > solveFEM(GModel* m, int nbNodes, Param param, Physical physical);
void computeK(gmm::row_matrix< gmm::wsvector< std::complex<double> > > &Ktmp, GModel::eiter eBegin, GModel::eiter eEnd, int nbNodes, Param param);

#endif /* fem_hp*/
