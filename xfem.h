#ifndef xfem_h
#define xfem_h

#include "GModel.h"
#include "GEdge.h"
#include "gmm.h"
#include "struct.h"

namespace XFEM {
    
std::vector< std::complex<double> > solve(GModel* m, int nbNodes, Param param, Physical physical);
void computeK(gmm::row_matrix< gmm::wsvector< std::complex<double> > > &Ktmp, GModel::eiter eBegin, GModel::eiter eEnd, int nbNodes, Param param);
    
}

#endif /* xfem_h */
