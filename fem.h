#ifndef fem_h
#define fem_h

#include "GModel.h"

#include "gmm.h"
#include "struct.h"

namespace FEM {
    
std::vector< std::complex<double> > solve(GModel* m, int nbNodes, Param param, Physical physical);
void computeK(gmm::row_matrix< gmm::wsvector< std::complex<double> > > &Ktmp, GModel::eiter eBegin, GModel::eiter eEnd, int nbNodes, Param param);
void computeK(gmm::row_matrix< gmm::wsvector< std::complex<double> > > &Ktmp, GModel::fiter fBegin, GModel::fiter fEnd, int nbNodes, Param param);
    
}


#endif /* fem_hp*/
