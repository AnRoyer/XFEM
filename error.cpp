#include "error.h"

std::vector< std::complex<double> > error(const std::vector< std::complex<double> > ue, std::vector< std::complex<double> > uc)
{
    std::vector< std::complex<double> > e(ue.size());
    if(ue.size() != uc.size())
    {
        return e;
    }
    
    for(unsigned int i = 0; i < ue.size(); i++)
    {
        e[i] = (ue[i] - uc[i])/ue[i];
    }
    
    return e;
}
