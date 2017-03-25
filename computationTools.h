#ifndef computationTools_h
#define computationTools_h

#include "MElement.h"

double BF_order1(unsigned int dim, unsigned int num, double u, double v = 0);
double BFgrad_order1(unsigned int dim, unsigned int num, double u, double v = 0);

void setLsu(int num, double value);
double F_enrichment(unsigned int dim, double u, double v = 0);
double Fgrad_enrichment(unsigned int dim, double u, double v = 0);

double jacobian1D(MElement* e);

double integraleK(unsigned int dim, unsigned int i, unsigned int j, unsigned int nbPoint);
double integraleM(unsigned int dim, unsigned int i, unsigned int j, unsigned int nbPoint);


#endif /* computationTools_h */
