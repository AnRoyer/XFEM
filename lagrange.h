#ifndef lagrange_h
#define lagrange_h

#include "GModel.h"
#include "GEdge.h"
#include "gmm.h"
#include "struct.h"

namespace LAGRANGE {
std::vector< std::complex<double> > solve(GModel* m, Param param, Physical physical);
void computeK(gmm::row_matrix< gmm::wsvector< std::complex<double> > > &Ktmp, std::vector<GEntity*> elms, Param param, std::map<int, std::vector<int> > &enrichedNodes, std::map<int, std::vector<int> > &lagrangianNodes);
    
}

#endif /* lagrange_h */
