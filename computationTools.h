#ifndef computationTools_h
#define computationTools_h

#include <vector>

#include "gmm.h"

#include "GEntity.h"
#include "MElement.h"

#include "struct.h"

enum Enrichment{
  None,
  Hat,
  Heaviside
};

std::vector<double> prodMatVec(gmm::dense_matrix<double> &J, std::vector<double> v);

double BFE_order1(Enrichment enrichment, unsigned int dim, unsigned int num, double u, double v = 0);
std::vector<double> BFEgrad_order1(Enrichment enrichment, unsigned int dim, unsigned int num, double u, double v = 0);

void setLsu(int num, double value);
double F_enrichment(unsigned int dim, double u, double v = 0);
std::vector<double> Fgrad_enrichment(unsigned int dim, double u, double v = 0);
double H_enrichment(unsigned int dim, unsigned int num, double u, double v = 0);
std::vector<double> Hgrad_enrichment(unsigned int dim, double u, double v = 0);

gmm::dense_matrix<double> jacobianVol(MElement* e, int dim);
gmm::dense_matrix<double> jacobianSur(MElement* e, int dim);

double integraleK(unsigned int dim, unsigned int i, unsigned int j, unsigned int nbPoint, gmm::dense_matrix<double> J, Enrichment enrichment = None);
double integraleM(unsigned int dim, unsigned int i, unsigned int j, unsigned int nbPoint, gmm::dense_matrix<double> J, Enrichment enrichment = None);
double integraleLD(unsigned int dim, unsigned int i, unsigned int j, unsigned int nbPoint, gmm::dense_matrix<double> J, Enrichment enrichment = None);
double integraleLN(unsigned int dim, unsigned int i, unsigned int j, unsigned int nbPoint, gmm::dense_matrix<double> J, Enrichment enrichment = None);


void sommerfeldCondition(gmm::row_matrix< gmm::wsvector< std::complex<double> > > &Ktmp, std::vector<GEntity*> elmInf, int dim, double k);


#endif /* computationTools_h */
