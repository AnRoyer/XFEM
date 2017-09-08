#ifndef analytical_h
#define analytical_h

#include "GModel.h"
#include "struct.h"

namespace ANALYTICAL {
    
std::vector< std::complex<double> > solve(GModel* m, Param param, Physical physical, bool xfem);
    
}

#endif /* analytical_h */
