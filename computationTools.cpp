#include <complex>

#include "computationTools.h"
#include "gauss.h"
#include "basicFunctions.h"

#include "GVertex.h"
#include "GEdge.h"

#include "MPoint.h"
#include "MLine.h"
#include "MTriangle.h"

double lsu[3];

std::vector<double> prodMatVec(gmm::dense_matrix<double> &J, std::vector<double> v)
{
  std::vector<double> result(2, 0.);
  
  for(unsigned int i = 0; i < 2; i++)
  {
    for(unsigned int j = 0; j < 2; j++)
    {
      result[i] += J(i,j)*v[j];
    }
  }
  
  return result;
}

double BFE_order1(Enrichment enrichment, unsigned int dim, unsigned int num, double u, double v)
{
  double value = 0.;
  
  if(enrichment == Hat)
  {
    if(num < 2)
    {
      value = BF_order1(dim, num, u, v);
    }
    else
    {
      value = F_enrichment(dim, u, v)*BF_order1(dim, num%2, u, v);
    }
  }
  else if(enrichment == Heaviside)
  {
    value = H_enrichment(dim, num, u, v)*BF_order1(dim, (int)num/2, u, v);
  }
  
  return value;
}

std::vector<double> BFEgrad_order1(Enrichment enrichment, unsigned int dim, unsigned int num, double u, double v)
{
  std::vector<double> value(2,0.);
  
  if(enrichment == Hat)
  {
    if(num < 2)
    {
      value = BFgrad_order1(dim, num, u, v);
    }
    else
    {
      std::vector<double> BFgrad = BFgrad_order1(dim, num%2, u, v);
      std::vector<double> Fgrad = Fgrad_enrichment(dim, u, v);
      
      for(unsigned int i = 0; i < 2; i++)
      {
        value[i] = F_enrichment(dim, u, v)*BFgrad[i] + Fgrad[i]*BF_order1(dim, num%2, u, v);
      }
    }
  }
  else if(enrichment == Heaviside)
  {
    std::vector<double> BFgrad = BFgrad_order1(dim, (int)num/2, u, v);
    
    for(unsigned int i = 0; i < 2; i++)
    {
      value[i] = H_enrichment(dim, num, u, v)*BFgrad[i];
    }
  }
  
  return value;
}

void setLsu(int num, double value)
{
  lsu[num] = value;
}

//Enrichment function : sum_i(|ls_i|N_i(x)) - |sum_i(ls_i*N_i(x))|
double F_enrichment(unsigned int dim, double u, double v)
{
  double value = 0.;
  
  double term1 = 0.;
  double term2 = 0.;
  
  for(unsigned int i = 0; i < dim+1; i++)
  {
    term1 += std::abs(lsu[i])*BF_order1(dim, i, u, v);
    term2 += lsu[i]*BF_order1(dim, i, u, v);
  }
  
  value = term1 - std::abs(term2);
   
  return value;
}

std::vector<double> Fgrad_enrichment(unsigned int dim, double u, double v)
{
  std::vector<double> value(2, 0.);
  
  std::vector<double> term1(2, 0.);
  std::vector<double> term2(2, 0.);
  double signTerm2 = 0.;
  
  for(unsigned int i = 0; i < dim+1; i++)
  {
    std::vector<double> BFgrad = BFgrad_order1(dim, i, u, v);
    
    term1[0] += std::abs(lsu[i])*BFgrad[0];
    term1[1] += std::abs(lsu[i])*BFgrad[1];
    term2[0] += lsu[i]*BFgrad[0];
    term2[1] += lsu[i]*BFgrad[1];
    signTerm2 += lsu[i]*BF_order1(dim, i, u, v);
  }
  
  if(signTerm2 < 0)
  {
    for(unsigned int i = 0; i < dim+1; i++)
    {
      term2[i] = -term2[i];
    }
  }
  
  for(unsigned int i = 0; i < dim+1; i++)
  {
    value[i] = term1[i] - term2[i];
  }
   
  return value;
}

double H_enrichment(unsigned int dim, unsigned int num, double u, double v)
{
  double value = 0.;
  
  double uBnd = -2.*lsu[0]/(lsu[1]-lsu[0])-1.;
  
  if(num == 0 || num == 2)
  {
    /*  ----------|
     *            |
     *            |
     *            |
     *            |-----------
     */
    
    if(u < uBnd) value = 1.;
    if(u == uBnd) value = 1.;
    if(u > uBnd) value = 0.;
  }
  else if(num == 1 || num == 3)
  {
    /*            |-----------
     *            |
     *            |
     *            |
     *  ----------|
     */
    if(u < uBnd) value = 0.;
    if(u == uBnd) value = 1.;
    if(u > uBnd) value = 1.;
  }
  
  return value;
}

std::vector<double> Hgrad_enrichment(unsigned int dim, double u, double v)
{
  std::vector<double> value(2, 0.);
  
  return value;
}


gmm::dense_matrix<double> jacobianVol(MElement* e, int dim)
{
  gmm::dense_matrix<double> J(2,2);
  if(dim == 0)
  {
    J(0,0) = 1.0;
    J(1,0) = 0.0;
    J(0,1) = 0.0;
    J(1,1) = 1.0;
  }
  else if(dim == 1)
  {
    double x0 = e->getVertex(0)->x();
    double x1 = e->getVertex(1)->x();
    
    J(0,0) = (x1-x0)/2;
    J(1,0) = 0.0;
    J(0,1) = 0.0;
    J(1,1) = 1.0;
  }
  else if(dim == 2)
  {
    double x0 = e->getVertex(0)->x();
    double x1 = e->getVertex(1)->x();
    double x2 = e->getVertex(2)->x();
    
    double y0 = e->getVertex(0)->y();
    double y1 = e->getVertex(1)->y();
    double y2 = e->getVertex(2)->y();
    
    J(0,0) = x1-x0;
    J(1,0) = x2-x0;
    J(0,1) = y1-y0;
    J(1,1) = y2-y0;
  }
  
  return J;
}

gmm::dense_matrix<double> jacobianSur(MElement* e, int dim)
{
  gmm::dense_matrix<double> J(2,2);
  if(dim == 0)
  {
    J(0,0) = 1.0;
    J(1,0) = 0.0;
    J(0,1) = 0.0;
    J(1,1) = 1.0;
  }
  else if(dim == 1)
  {
    double x0 = e->getVertex(0)->x();
    double x1 = e->getVertex(1)->x();
    
    double y0 = e->getVertex(0)->y();
    double y1 = e->getVertex(1)->y();
    
    double value = sqrt((x1-x0)*(x1-x0)/4 + (y1-y0)*(y1-y0)/4);
    
    J(0,0) = value;
    J(1,0) = 0.0;
    J(0,1) = 0.0;
    J(1,1) = 1.0;
  }
  
  return J;
}

double integraleK(unsigned int dim, unsigned int i, unsigned int j, unsigned int nbPoint, gmm::dense_matrix<double> J, Enrichment enrichment)
{
  double value = 0.;
  double detJ = gmm::lu_inverse(J);
  
  switch (dim) {
    case 0:
    {
      std::vector<double>BFgradI = prodMatVec(J, BFgrad_order1(dim, i, px1[0]));
      std::vector<double>BFgradJ = prodMatVec(J, BFgrad_order1(dim, j, px1[0]));
      
      value += pp1[0]*(BFgradI[0]*BFgradJ[0]+BFgradI[1]*BFgradJ[1])*std::abs(detJ);
    }
      break;
    case 1:
    {
      std::vector<double*> lx;
      lx.push_back(lx1);
      lx.push_back(lx2);
      lx.push_back(lx3);
      lx.push_back(lx4);
      lx.push_back(lx5);
      lx.push_back(lx6);
      lx.push_back(lx7);
      lx.push_back(lx8);
      lx.push_back(lx9);
      lx.push_back(lx10);
      lx.push_back(lx11);
      lx.push_back(lx12);
      lx.push_back(lx13);
      lx.push_back(lx14);
      lx.push_back(lx15);
      lx.push_back(lx16);
      lx.push_back(lx17);
      lx.push_back(lx18);
      lx.push_back(lx19);
      lx.push_back(lx20);
    
      std::vector<double*> lp;
      lp.push_back(lp1);
      lp.push_back(lp2);
      lp.push_back(lp3);
      lp.push_back(lp4);
      lp.push_back(lp5);
      lp.push_back(lp6);
      lp.push_back(lp7);
      lp.push_back(lp8);
      lp.push_back(lp9);
      lp.push_back(lp10);
      lp.push_back(lp11);
      lp.push_back(lp12);
      lp.push_back(lp13);
      lp.push_back(lp14);
      lp.push_back(lp15);
      lp.push_back(lp16);
      lp.push_back(lp17);
      lp.push_back(lp18);
      lp.push_back(lp19);
      lp.push_back(lp20);
    
      if(enrichment == None)
      {
        for(unsigned int k = 0; k < nbPoint; k++)
        {
          std::vector<double>BFgradI = prodMatVec(J, BFgrad_order1(dim, i, lx[nbPoint-1][k]));
          std::vector<double>BFgradJ = prodMatVec(J, BFgrad_order1(dim, j, lx[nbPoint-1][k]));
          
          value += lp[nbPoint-1][k]*(BFgradI[0]*BFgradJ[0]+BFgradI[1]*BFgradJ[1])*std::abs(detJ);
        }
      }
      else
      {
        for(unsigned int k = 0; k < nbPoint; k++)
        {
          std::vector<double>BFgradI = prodMatVec(J, BFEgrad_order1(enrichment, dim, i, lx[nbPoint-1][k]));
          std::vector<double>BFgradJ = prodMatVec(J, BFEgrad_order1(enrichment, dim, j, lx[nbPoint-1][k]));
          
          value += lp[nbPoint-1][k]*(BFgradI[0]*BFgradJ[0]+BFgradI[1]*BFgradJ[1])*std::abs(detJ);
        }
      }
    }
    break;
    default:
      break;
  }
  
  return value;
}

double integraleM(unsigned int dim, unsigned int i, unsigned int j, unsigned int nbPoint, gmm::dense_matrix<double> J, Enrichment enrichment)
{
  double detJ = gmm::lu_det(J);
  
  double value = 0.;
  switch (dim) {
    case 0:
      value += pp1[0]*BF_order1(dim, i, px1[0])*BF_order1(dim, j, px1[0])*std::abs(detJ);
      break;
    case 1:
    {
      std::vector<double*> lx;
      lx.push_back(lx1);
      lx.push_back(lx2);
      lx.push_back(lx3);
      lx.push_back(lx4);
      lx.push_back(lx5);
      lx.push_back(lx6);
      lx.push_back(lx7);
      lx.push_back(lx8);
      lx.push_back(lx9);
      lx.push_back(lx10);
      lx.push_back(lx11);
      lx.push_back(lx12);
      lx.push_back(lx13);
      lx.push_back(lx14);
      lx.push_back(lx15);
      lx.push_back(lx16);
      lx.push_back(lx17);
      lx.push_back(lx18);
      lx.push_back(lx19);
      lx.push_back(lx20);
    
      std::vector<double*> lp;
      lp.push_back(lp1);
      lp.push_back(lp2);
      lp.push_back(lp3);
      lp.push_back(lp4);
      lp.push_back(lp5);
      lp.push_back(lp6);
      lp.push_back(lp7);
      lp.push_back(lp8);
      lp.push_back(lp9);
      lp.push_back(lp10);
      lp.push_back(lp11);
      lp.push_back(lp12);
      lp.push_back(lp13);
      lp.push_back(lp14);
      lp.push_back(lp15);
      lp.push_back(lp16);
      lp.push_back(lp17);
      lp.push_back(lp18);
      lp.push_back(lp19);
      lp.push_back(lp20);
    
      if(enrichment == None)
      {
        for(unsigned int k = 0; k < nbPoint; k++)
        {
          double BFEi = BF_order1(dim, i, lx[nbPoint-1][k]);
          double BFEj = BF_order1(dim, j, lx[nbPoint-1][k]);
          
          value += lp[nbPoint-1][k]*BFEi*BFEj*std::abs(detJ);
        }
      }
      else
      {
        for(unsigned int k = 0; k < nbPoint; k++)
        {
          double BFEi = BFE_order1(enrichment, dim, i, lx[nbPoint-1][k]);
          double BFEj = BFE_order1(enrichment, dim, j, lx[nbPoint-1][k]);
              
          value += lp[nbPoint-1][k]*BFEi*BFEj*std::abs(detJ);
        }
      }
    }
      break;
    default:
      break;
  }
  
  return value;
}

double integraleLD(unsigned int dim, unsigned int i, unsigned int j, unsigned int nbPoint, gmm::dense_matrix<double> J, Enrichment enrichment)
{
  double detJ = gmm::lu_det(J);
  
  double value = 0.;
  switch (dim-1) {
    case 0:
    {
      const double uBnd = -2.*lsu[0]/(lsu[1]-lsu[0])-1.;
      value += pp1[0]*BF_order1(dim-1, i, uBnd)*BFE_order1(enrichment, dim, j, uBnd);
    }
      break;
    default:
      break;
  }
  
  return value;
}

double integraleLN(unsigned int dim, unsigned int i, unsigned int j, unsigned int nbPoint, gmm::dense_matrix<double> J, Enrichment enrichment)
{
  double detJ = gmm::lu_det(J);
  
  double value = 0.;
  switch (dim-1) {
    case 0:
    {
      const double uBnd = -2.*lsu[0]/(lsu[1]-lsu[0])-1.;
      value += pp1[0]*BF_order1(dim-1, i, uBnd)*BFEgrad_order1(enrichment, dim, j, uBnd)[0];
    }
      break;
    default:
      break;
  }
  
  return value;
}


void sommerfeldCondition(gmm::row_matrix< gmm::wsvector< std::complex<double> > > &Ktmp, std::vector<GEntity*> elmInf, int dim, double k)
{
  for(unsigned int i = 0; i < elmInf.size(); i++)
  {
    if(dim == 1)
    {
      for(std::vector<MPoint*>::iterator itP = static_cast<GVertex*>(elmInf[i])->points.begin(); itP != static_cast<GVertex*>(elmInf[i])->points.end(); ++itP)
      {
        MPoint* point = *itP;
        int nbN0 = (*itP)->getVertex(0)->getNum()-1;
        
        gmm::dense_matrix<double> J = jacobianSur(point, dim-1);
        
        Ktmp(nbN0, nbN0) += std::complex<double>(0., k*integraleM(dim-1, 0, 0, 1, J));
      }
    }
    else if(dim == 2)
    {
      for(std::vector<MLine*>::iterator itL = static_cast<GEdge*>(elmInf[i])->lines.begin(); itL != static_cast<GEdge*>(elmInf[i])->lines.end(); ++itL)
      {
        MLine* line = *itL;
        int nbN0 = line->getVertex(0)->getNum()-1;
        int nbN1 = line->getVertex(1)->getNum()-1;
        
        gmm::dense_matrix<double> Ke(2,2);
        gmm::dense_matrix<double> J = jacobianSur(line, dim-1);
        
        for(unsigned int i = 0; i < 2; i++)
        {
          for(unsigned int j = 0; j < 2; j++)
          {
            Ke(i,j) = k*integraleM(dim-1, i, j, 3, J);
          }
        }
        
        Ktmp(nbN0, nbN0) += std::complex<double>(0., Ke(0,0));
        Ktmp(nbN0, nbN1) += std::complex<double>(0., Ke(0,1));
        Ktmp(nbN1, nbN0) += std::complex<double>(0., Ke(1,0));
        Ktmp(nbN1, nbN1) += std::complex<double>(0., Ke(1,1));
      }
    }
  }
}


