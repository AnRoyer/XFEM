#ifndef fem_h
#define fem_h

#include "GModel.h"

#include "gmm.h"
#include "struct.h"

namespace FEM {
    
std::vector< std::complex<double> > solve(GModel* m, Param param, Physical physical);
void computeK(gmm::row_matrix< gmm::wsvector< std::complex<double> > > &Ktmp, std::vector<GEntity*> elms, double k, double rho, double c);
}


#endif /* fem_hp*/
