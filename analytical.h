#ifndef analytical_h
#define analytical_h

#include "GModel.h"
#include "struct.h"

namespace ANALYTICAL {
    
std::vector< std::complex<double> > solve(GModel* m, int nbNodes, Param param, Physical physical);
    
}

#endif /* analytical_h */
