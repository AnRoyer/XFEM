#ifndef xfem_h
#define xfem_h

#include "GModel.h"
#include "GEdge.h"
#include "gmm.h"
#include "struct.h"

namespace XFEM {
std::vector< std::complex<double> > solve(GModel* m, Param param, Physical physical);
void computeK(gmm::row_matrix< gmm::wsvector< std::complex<double> > > &Ktmp, std::vector<GEntity*> elms, double k1, double rho1, double c1, double k2, double rho2, double c2, double bnd, int nbNodes);
    
}

#endif /* xfem_h */
