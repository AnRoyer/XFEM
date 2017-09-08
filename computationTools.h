#ifndef computationTools_h
#define computationTools_h

#include <vector>
#include <map>

#include "gmm.h"

#include "GEntity.h"
#include "MElement.h"
#include "GModel.h"

#include "struct.h"

enum Enrichment{
  None,
  Hat,
  Heaviside
};

std::vector<double> prodMatVec(gmm::dense_matrix<double> &J, std::vector<double> v);

double BFE_order1(Enrichment enrichment, unsigned int dim, unsigned int num, double u, double v = 0);
std::vector<double> BFEgrad_order1(Enrichment enrichment, unsigned int dim, unsigned int num, double u, double v = 0);

void setLst(int num, double value);
double F_enrichment(unsigned int dim, double u, double v = 0);
std::vector<double> Fgrad_enrichment(unsigned int dim, double u, double v = 0);
double H_enrichment(unsigned int dim, unsigned int num, double u, double v = 0);
std::vector<double> Hgrad_enrichment(unsigned int dim, double u, double v = 0);

gmm::dense_matrix<double> jacobianVol(MElement* e, int dim);
gmm::dense_matrix<double> jacobianSur(MElement* e, int dim);
gmm::dense_matrix<double> jacobianHeaviside(int dim, int num);
gmm::dense_matrix<double> jacobianLagrange(int dim, gmm::dense_matrix<double> J);

double integraleK(unsigned int dim, unsigned int i, unsigned int j, unsigned int nbPoint, gmm::dense_matrix<double> J, Enrichment enrichment = None);
double integraleM(unsigned int dim, unsigned int i, unsigned int j, unsigned int nbPoint, gmm::dense_matrix<double> J, Enrichment enrichment = None);
double integraleLD(unsigned int dim, unsigned int i, unsigned int j, unsigned int nbPoint, gmm::dense_matrix<double> J, Enrichment enrichment = None);

void sommerfeldCondition(gmm::row_matrix< gmm::wsvector< std::complex<double> > > &Ktmp, std::vector<GEntity*> elmInf, int dim, double k);

void getTotalNumbersOfNodes(GModel* m, Param param, Physical physical, int *nbTotNodes, std::map<int, std::vector<int> > *enrichedNodes);
void getNumbersOfLagrangian(GModel* m, Param param, Physical physical, int nbNodes, std::map<int, std::vector<int> > &enrichedNodes, std::map<int, std::vector<int> > *lagrangianNodes, int *nbLag);
std::map<int, std::vector<int> > getNeighbor(GModel* m, Param param, std::map<int, std::vector<int> > &enrichedNodes);

#endif /* computationTools_h */
