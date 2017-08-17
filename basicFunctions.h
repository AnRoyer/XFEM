#ifndef basicFunctions_h
#define basicFunctions_h

#include <vector>

#include "gmm.h"

#include "GEntity.h"
#include "MElement.h"

#include "struct.h"

double BF_order1(unsigned int dim, unsigned int num, double u, double v = 0);
std::vector<double> BFgrad_order1(unsigned int dim, unsigned int num, double u, double v = 0);



#endif /* basicFunctions_h */
