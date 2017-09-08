#include <complex>
#include <cmath>

#include "computationTools.h"
#include "gauss.h"
#include "basicFunctions.h"

#include "GVertex.h"
#include "GEdge.h"
#include "GFace.h"

#include "MPoint.h"
#include "MLine.h"
#include "MTriangle.h"

double lst[3];

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
      value = F_enrichment(dim, u, v)*BF_order1(dim, num%(dim+1), u, v);
    }
  }
  else if(enrichment == Heaviside)
  {
    value = H_enrichment(dim, num, u, v)*BF_order1(dim, num%(dim+1), u, v);
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
      std::vector<double> BFgrad = BFgrad_order1(dim, num%(dim+1), u, v);
      std::vector<double> Fgrad = Fgrad_enrichment(dim, u, v);
      
      for(unsigned int i = 0; i < 2; i++)
      {
        value[i] = F_enrichment(dim, u, v)*BFgrad[i] + Fgrad[i]*BF_order1(dim, num%(dim+1), u, v);
      }
    }
  }
  else if(enrichment == Heaviside)
  {
    std::vector<double> BFgrad = BFgrad_order1(dim, num%(dim+1), u, v);
    
    for(unsigned int i = 0; i < 2; i++)
    {
      value[i] = H_enrichment(dim, num, u, v)*BFgrad[i];
    }
  }
  
  return value;
}

void setLst(int num, double value)
{
  lst[num] = value;
}

//Enrichment function : sum_i(|ls_i|N_i(x)) - |sum_i(ls_i*N_i(x))|
double F_enrichment(unsigned int dim, double u, double v)
{
  double value = 0.;
  
  double term1 = 0.;
  double term2 = 0.;
  
  for(unsigned int i = 0; i < dim+1; i++)
  {
    term1 += std::abs(lst[i])*BF_order1(dim, i, u, v);
    term2 += lst[i]*BF_order1(dim, i, u, v);
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
    
    term1[0] += std::abs(lst[i])*BFgrad[0];
    term1[1] += std::abs(lst[i])*BFgrad[1];
    term2[0] += lst[i]*BFgrad[0];
    term2[1] += lst[i]*BFgrad[1];
    signTerm2 += lst[i]*BF_order1(dim, i, u, v);
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
  
  if(dim == 1)
  {
    /*  Equation : f(u) = a*u + b
     *
     *  System :
     *   -a + b = lst[0]
     *    a + b = lst[1]
     *
     *    => a = (lst[1] - lst[0])/2;
     *       b = (lst[0] + lst[1])/2;
     *
     *    => bnd = {(u) : (lst[1] - lst[0])*u/2 + (lst[0] + lst[1])/2 = 0}
     *
     */
    const double bnd = (lst[1] - lst[0])*u/2. + (lst[0] + lst[1])/2.;
    //const int order[4] = {0,1,1,0};
    
    if(lst[num%2] < 0)
    {
      if(num/2 == 0)
      {
        if(bnd < 0)  value = 1.;
        if(bnd == 0) value = 0.5;
        if(bnd > 0)  value = 0.;
      }
      else
      {
        if(bnd < 0)  value = 0.;
        if(bnd == 0) value = 0.5;
        if(bnd > 0)  value = 1.;
      }
    }
    else
    {
      if(num/2 == 0)
      {
        if(bnd > 0)  value = 1.;
        if(bnd == 0) value = 0.5;
        if(bnd < 0)  value = 0.;
      }
      else
      {
        if(bnd > 0)  value = 0.;
        if(bnd == 0) value = 0.5;
        if(bnd < 0)  value = 1.;
      }
    }
  }
  else if(dim == 2)
  {
    /*  Equation : f(u) = a*u + b*v + c
     *
     *  System :
     *            c = lst[0]
     *    a     + c = lst[1]
     *        b + c = lst[2]
     *
     *    => a = lst[1] - lst[0];
     *       b = lst[2] - lst[0];
     *       c = lst[0]
     *
     *    => bnd = {(u,v) : (lst[1] - lst[0])*u + (lst[2] - lst[0])*v + lst[0] = 0};
     *
     */
    double bnd = (lst[1] - lst[0])*u + (lst[2] - lst[0])*v + lst[0];
    //const int order[6] = {0,1,2,(lst[0]*lst[1] < 0 ? 1 : 2),(lst[1]*lst[2] < 0 ? 2 : 0),(lst[2]*lst[0] < 0 ? 0 : 1)};
    
    if(lst[num%3] < 0)
    {
      if(num/3 == 0)
      {
        if(bnd < 0)  value = 1.;
        if(bnd == 0) value = 0.5;
        if(bnd > 0)  value = 0.;
      }
      else
      {
        if(bnd < 0)  value = 0.;
        if(bnd == 0) value = 0.5;
        if(bnd > 0)  value = 1.;
      }
    }
    else
    {
      if(num/3 == 0)
      {
        if(bnd > 0)  value = 1.;
        if(bnd == 0) value = 0.5;
        if(bnd < 0)  value = 0.;
      }
      else
      {
        if(bnd > 0)  value = 0.;
        if(bnd == 0) value = 0.5;
        if(bnd < 0)  value = 1.;
      }
    }
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

gmm::dense_matrix<double> jacobianHeaviside(int dim, int num)
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
    /*  The line is split in 2 lines
     *
     *  J(0,0) = u1-u0;
     *  J(1,0) = 0.0;
     *  J(0,1) = 0.0;
     *  J(1,1) = 1.0;
     */
    const double uBnd = -((lst[0] + lst[1])/2.)/((lst[1] - lst[0])/2.);
    if(num == 0)
    {
      J(0,0) = (uBnd+1)/2;
      J(1,0) = 0.0;
      J(0,1) = 0.0;
      J(1,1) = 1.0;
    }
    else if(num == 1)
    {
      J(0,0) = (uBnd-1)/2;
      J(1,0) = 0.0;
      J(0,1) = 0.0;
      J(1,1) = 1.0;
    }
  }
  else if(dim == 2)
  {
    /*  The triangle is split in 2 triangles OR 3 triangles
     *
     *  J(0,0) = u1-u0;
     *  J(1,0) = u2-u0;
     *  J(0,1) = v1-v0;
     *  J(1,1) = v2-v0;
     *
     *  a = lst[1] - lst[0];
     *  b = lst[2] - lst[0];
     *  c = lst[0]
     */
    if(lst[0] == 0)//Only two triangles
    {
      const double u0Bnd = -lst[2]/(lst[1] - lst[2]);
      const double v0Bnd = 1. - u0Bnd;
      if(num == 0)
      {
        J(0,0) = 0;
        J(1,0) = 0;
        J(0,1) = 0;
        J(1,1) = 0;
      }
      else if(num == 1)
      {
        J(0,0) = u0Bnd-1;
        J(1,0) = 0-1;
        J(0,1) = v0Bnd-0;
        J(1,1) = 0-0;
      }
      else if(num == 2)
      {
        J(0,0) = 0-0;
        J(1,0) = u0Bnd-0;
        J(0,1) = 0-1;
        J(1,1) = v0Bnd-1;
      }
    }
    else if(lst[1] == 0)//Only two triangles
    {
      const double u0Bnd = 0.;
      const double v0Bnd = -(lst[0])/(lst[2] - lst[0]);
      if(num == 0)
      {
        J(0,0) = 1-0;
        J(1,0) = u0Bnd-0;
        J(0,1) = 0-0;
        J(1,1) = v0Bnd-0;
      }
      else if(num == 1)
      {
        J(0,0) = 0;
        J(1,0) = 0;
        J(0,1) = 0;
        J(1,1) = 0;
      }
      else if(num == 2)
      {
        J(0,0) = 0-0;
        J(1,0) = u0Bnd-0;
        J(0,1) = 0-1;
        J(1,1) = v0Bnd-1;
      }
    }
    else if(lst[2] == 0)//Only two triangles
    {
      const double u0Bnd = -(lst[0])/(lst[1] - lst[0]);
      const double v0Bnd = 0.;
      if(num == 0)
      {
        J(0,0) = u0Bnd-0;
        J(1,0) = 0-0;
        J(0,1) = v0Bnd-0;
        J(1,1) = 1-0;
      }
      else if(num == 1)
      {
        J(0,0) = 0-1;
        J(1,0) = u0Bnd-1;
        J(0,1) = 1-0;
        J(1,1) = v0Bnd-0;
      }
      else if(num == 2)
      {
        J(0,0) = 0;
        J(1,0) = 0;
        J(0,1) = 0;
        J(1,1) = 0;
      }
    }
    else//Otherwise there are 3 triangles
    {
      double u0Bnd = 0.;
      double v0Bnd = -(lst[0])/(lst[2] - lst[0]);

      double u1Bnd = -lst[2]/(lst[1] - lst[2]);
      double v1Bnd = 1. - u0Bnd;
      
      double u2Bnd = -(lst[0])/(lst[1] - lst[0]);
      double v2Bnd = 0.;
      
      //Check if there are correct
      if(!isfinite(v0Bnd) || v0Bnd < 0 || v0Bnd > 1)
      {
        u0Bnd = INFINITY;
        v0Bnd = INFINITY;
      }
      if(!isfinite(v1Bnd) || v1Bnd < 0 || v1Bnd > 1 || u1Bnd < 0 || u1Bnd > 1)
      {
        u1Bnd = INFINITY;
        v1Bnd = INFINITY;
      }
      if(!isfinite(u2Bnd) || u2Bnd < 0 || u2Bnd > 1)
      {
        u2Bnd = INFINITY;
        v2Bnd = INFINITY;
      }
      
      if(num == 0)
      {
        if(isfinite(u0Bnd) && isfinite(u2Bnd))
        {
          J(0,0) = u0Bnd-0;
          J(1,0) = u2Bnd-0;
          J(0,1) = v0Bnd-0;
          J(1,1) = v2Bnd-0;
        }
        else
        {
          if(isfinite(u0Bnd))
          {
            J(0,0) = u0Bnd-0;
            J(1,0) = u1Bnd-0;
            J(0,1) = v0Bnd-0;
            J(1,1) = v1Bnd-0;
            if(!isfinite(J(0,0)) || !isfinite(J(0,1)) || !isfinite(J(1,0)) || !isfinite(J(1,1)))
            {
              std::cout << J(0,0) << std::endl;
              std::cout << J(1,0) << std::endl;
              std::cout << J(0,1) << std::endl;
              std::cout << J(1,1) << std::endl;
            }
          }
          else//isfinite(u2Bnd)
          {
            J(0,0) = u1Bnd-0;
            J(1,0) = u2Bnd-0;
            J(0,1) = v1Bnd-0;
            J(1,1) = v2Bnd-0;
          }
        }
      }
      else if(num == 1)
      {
        if(isfinite(u0Bnd) && isfinite(u1Bnd))
        {
          J(0,0) = u1Bnd-1;
          J(1,0) = u0Bnd-1;
          J(0,1) = v1Bnd-0;
          J(1,1) = v0Bnd-0;
        }
        else
        {
          if(isfinite(u1Bnd))
          {
            J(0,0) = u1Bnd-1;
            J(1,0) = 0-1;
            J(0,1) = v1Bnd-0;
            J(1,1) = 0-0;
          }
          else//isfinite(u0Bnd)
          {
            J(0,0) = u2Bnd-1;
            J(1,0) = u0Bnd-1;
            J(0,1) = v2Bnd-0;
            J(1,1) = v0Bnd-0;
          }
        }
      }
      else if(num == 2)
      {
        if(isfinite(u1Bnd) && isfinite(u2Bnd))
        {
          J(0,0) = u2Bnd-0;
          J(1,0) = u1Bnd-0;
          J(0,1) = v2Bnd-1;
          J(1,1) = v1Bnd-1;
        }
        else
        {
          if(isfinite(u2Bnd))
          {
            J(0,0) = u2Bnd-0;
            J(1,0) = 1-0;
            J(0,1) = v2Bnd-1;
            J(1,1) = 0-1;
          }
          else//isfinite(u1Bnd)
          {
            J(0,0) = 0-0;
            J(1,0) = u1Bnd-0;
            J(0,1) = 0-1;
            J(1,1) = v1Bnd-1;
          }
        }
      }
    }
  }
  
  return J;
}

gmm::dense_matrix<double> jacobianLagrange(int dim, gmm::dense_matrix<double> J)
{
  gmm::dense_matrix<double> Jl(2,2);
  if(dim == 1)
  {
    double u0Bnd = 0.;
    double v0Bnd = -(lst[0])/(lst[2] - lst[0]);
    
    double u1Bnd = -lst[2]/(lst[1] - lst[2]);
    double v1Bnd = 1. - u0Bnd;
    
    double u2Bnd = -(lst[0])/(lst[1] - lst[0]);
    double v2Bnd = 0.;
    
    double u0 = 0.;
    double v0 = 0.;
    double u1 = 0.;
    double v1 = 0.;
    //Check if there are correct
    if(!isfinite(v0Bnd) || v0Bnd < 0 || v0Bnd > 1)
    {
      u0 = u1Bnd;
      v0 = v1Bnd;
      
      u1 = u2Bnd;
      v1 = v2Bnd;
    }
    if(!isfinite(v1Bnd) || v1Bnd < 0 || v1Bnd > 1 || u1Bnd < 0 || u1Bnd > 1)
    {
      u0 = u2Bnd;
      v0 = v2Bnd;
      
      u1 = u0Bnd;
      v1 = v0Bnd;
    }
    if(!isfinite(u2Bnd) || u2Bnd < 0 || u2Bnd > 1)
    {
      u0 = u0Bnd;
      v0 = v0Bnd;
      
      u1 = u1Bnd;
      v1 = v1Bnd;
    }
    
    std::vector<double> uv0 {u0,v0};
    std::vector<double> uv1 {u1,v1};
    std::vector<double> xy0 = prodMatVec(J, uv0);
    std::vector<double> xy1 = prodMatVec(J, uv1);
    Jl(0,0) = sqrt((xy1[0]-xy0[0])*(xy1[0]-xy0[0])/4 + (xy1[1]-xy0[1])*(xy1[1]-xy0[1])/4);
    Jl(1,0) = 0.0;
    Jl(0,1) = 0.0;
    Jl(1,1) = 1.0;
  }
  
  return Jl;
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
        const double uPiece[2] = {-1.,1};
        for(unsigned int piece = 0; piece < 2; piece++)
        {
          gmm::dense_matrix<double> J2 = jacobianHeaviside(1, piece);
          double detJ2 = gmm::lu_det(J2);
          if(detJ2 != 0)
          {
            for(unsigned int k = 0; k < nbPoint; k++)
            {
              const double u = J2(0,0)*lx[nbPoint-1][k]+uPiece[piece];
          
              std::vector<double>BFgradI = prodMatVec(J, BFEgrad_order1(enrichment, dim, i, u));
              std::vector<double>BFgradJ = prodMatVec(J, BFEgrad_order1(enrichment, dim, j, u));
          
              value += lp[nbPoint-1][k]*(BFgradI[0]*BFgradJ[0]+BFgradI[1]*BFgradJ[1])*std::abs(detJ)*std::abs(detJ2);
            }
          }
        }
      }
    }
    break;
    case 2:
    {
      std::vector<double*> tx;
      tx.push_back(tx1);
      tx.push_back(nullptr);
      tx.push_back(tx3);
      tx.push_back(tx4);
      tx.push_back(nullptr);
      tx.push_back(tx6);
      tx.push_back(tx7);
      tx.push_back(nullptr);
      tx.push_back(nullptr);
      tx.push_back(nullptr);
      tx.push_back(nullptr);
      tx.push_back(tx12);
      tx.push_back(tx13);
      tx.push_back(nullptr);
      tx.push_back(nullptr);
      tx.push_back(tx16);
      
      std::vector<double*> ty;
      ty.push_back(ty1);
      ty.push_back(nullptr);
      ty.push_back(ty3);
      ty.push_back(ty4);
      ty.push_back(nullptr);
      ty.push_back(ty6);
      ty.push_back(ty7);
      ty.push_back(nullptr);
      ty.push_back(nullptr);
      ty.push_back(nullptr);
      ty.push_back(nullptr);
      ty.push_back(ty12);
      ty.push_back(ty13);
      ty.push_back(nullptr);
      ty.push_back(nullptr);
      ty.push_back(ty16);
      
      std::vector<double*> tp;
      tp.push_back(tp1);
      tp.push_back(nullptr);
      tp.push_back(tp3);
      tp.push_back(tp4);
      tp.push_back(nullptr);
      tp.push_back(tp6);
      tp.push_back(tp7);
      tp.push_back(nullptr);
      tp.push_back(nullptr);
      tp.push_back(nullptr);
      tp.push_back(nullptr);
      tp.push_back(tp12);
      tp.push_back(tp13);
      tp.push_back(nullptr);
      tp.push_back(nullptr);
      tp.push_back(tp16);
      
      if(enrichment == None)
      {
        for(unsigned int k = 0; k < nbPoint; k++)
        {
          std::vector<double>BFgradI = prodMatVec(J, BFgrad_order1(dim, i, tx[nbPoint-1][k], ty[nbPoint-1][k]));
          std::vector<double>BFgradJ = prodMatVec(J, BFgrad_order1(dim, j, tx[nbPoint-1][k], ty[nbPoint-1][k]));
          
          value += tp[nbPoint-1][k]*(BFgradI[0]*BFgradJ[0]+BFgradI[1]*BFgradJ[1])*std::abs(detJ);
        }
      }
      else
      {
        const double uPiece[3] = {0.,1.,0.};
        const double vPiece[3] = {0.,0.,1.};
        for(unsigned int piece = 0; piece < 3; piece++)
        {
          gmm::dense_matrix<double> J2 = jacobianHeaviside(dim, piece);
          double detJ2 = gmm::lu_det(J2);
          if(detJ2 != 0)
          {
            for(unsigned int k = 0; k < nbPoint; k++)
            {
              const double u = J2(0,0)*tx[nbPoint-1][k] + J2(0,1)*ty[nbPoint-1][k] + uPiece[piece];
              const double v = J2(1,0)*tx[nbPoint-1][k] + J2(1,1)*ty[nbPoint-1][k] + vPiece[piece];
              
              std::vector<double>BFgradI = prodMatVec(J, BFEgrad_order1(enrichment, dim, i, u, v) );
              std::vector<double>BFgradJ = prodMatVec(J, BFEgrad_order1(enrichment, dim, j, u, v) );
              
              value += tp[nbPoint-1][k]*(BFgradI[0]*BFgradJ[0]+BFgradI[1]*BFgradJ[1])*std::abs(detJ)*std::abs(detJ2);
            }
          }
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
        const double uPiece[2] = {-1.,1.};
        for(unsigned int piece = 0; piece < 2; piece++)
        {
          gmm::dense_matrix<double> J2 = jacobianHeaviside(dim, piece);
          double detJ2 = gmm::lu_det(J2);
          if(detJ2 != 0)
          {
            for(unsigned int k = 0; k < nbPoint; k++)
            {
              const double u = J2(0,0)*lx[nbPoint-1][k]+uPiece[piece];
          
              double BFEi = BFE_order1(enrichment, dim, i, u);
              double BFEj = BFE_order1(enrichment, dim, j, u);
          
              value += lp[nbPoint-1][k]*BFEi*BFEj*std::abs(detJ)*std::abs(detJ2);
            }
          }
        }
      }
    }
      break;
    case 2:
    {
      std::vector<double*> tx;
      tx.push_back(tx1);
      tx.push_back(nullptr);
      tx.push_back(tx3);
      tx.push_back(tx4);
      tx.push_back(nullptr);
      tx.push_back(tx6);
      tx.push_back(tx7);
      tx.push_back(nullptr);
      tx.push_back(nullptr);
      tx.push_back(nullptr);
      tx.push_back(nullptr);
      tx.push_back(tx12);
      tx.push_back(tx13);
      tx.push_back(nullptr);
      tx.push_back(nullptr);
      tx.push_back(tx16);
      
      std::vector<double*> ty;
      ty.push_back(ty1);
      ty.push_back(nullptr);
      ty.push_back(ty3);
      ty.push_back(ty4);
      ty.push_back(nullptr);
      ty.push_back(ty6);
      ty.push_back(ty7);
      ty.push_back(nullptr);
      ty.push_back(nullptr);
      ty.push_back(nullptr);
      ty.push_back(nullptr);
      ty.push_back(ty12);
      ty.push_back(ty13);
      ty.push_back(nullptr);
      ty.push_back(nullptr);
      ty.push_back(ty16);
      
      std::vector<double*> tp;
      tp.push_back(tp1);
      tp.push_back(nullptr);
      tp.push_back(tp3);
      tp.push_back(tp4);
      tp.push_back(nullptr);
      tp.push_back(tp6);
      tp.push_back(tp7);
      tp.push_back(nullptr);
      tp.push_back(nullptr);
      tp.push_back(nullptr);
      tp.push_back(nullptr);
      tp.push_back(tp12);
      tp.push_back(tp13);
      tp.push_back(nullptr);
      tp.push_back(nullptr);
      tp.push_back(tp16);
      
      if(enrichment == None)
      {
        for(unsigned int k = 0; k < nbPoint; k++)
        {
          double BFEi = BF_order1(dim, i, tx[nbPoint-1][k], ty[nbPoint-1][k]);
          double BFEj = BF_order1(dim, j, tx[nbPoint-1][k], ty[nbPoint-1][k]);
          
          value += tp[nbPoint-1][k]*BFEi*BFEj*std::abs(detJ);
        }
      }
      else
      {
        const double uPiece[3] = {0.,1.,0.};
        const double vPiece[3] = {0.,0.,1.};
        for(unsigned int piece = 0; piece < 3; piece++)
        {
          gmm::dense_matrix<double> J2 = jacobianHeaviside(dim, piece);
          double detJ2 = gmm::lu_det(J2);
          if(detJ2 != 0)
          {
            for(unsigned int k = 0; k < nbPoint; k++)
            {
              const double u = J2(0,0)*tx[nbPoint-1][k] + J2(0,1)*ty[nbPoint-1][k] + uPiece[piece];
              const double v = J2(1,0)*tx[nbPoint-1][k] + J2(1,1)*ty[nbPoint-1][k] + vPiece[piece];
              
              double BFEi = BFE_order1(enrichment, dim, i, u, v);
              double BFEj = BFE_order1(enrichment, dim, j, u, v);
              
              value += tp[nbPoint-1][k]*BFEi*BFEj*std::abs(detJ)*std::abs(detJ2);
            }
          }
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
      const double uBnd = (lst[0] + lst[1])/(lst[0] - lst[1]);
      value += pp1[0]*BF_order1(dim-1, i, 0.)*BFE_order1(enrichment, dim, j, uBnd);
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
      
      const double a = lst[1] - lst[0];
      const double b = lst[2] - lst[0];
      const double c = lst[0];
      
      double u0Bnd = 0.;
      double v0Bnd = -c/b;
      
      double u1Bnd = -b/(a-b);
      double v1Bnd = 1. - u0Bnd;
      
      double u2Bnd = -c/a;
      double v2Bnd = 0.;
      
      //Check if there are correct
      if(!isfinite(v0Bnd) || v0Bnd < 0 || v0Bnd > 1)
      {
        u0Bnd = 0.;
        v0Bnd = 0.;
      }
      if(!isfinite(v1Bnd) || v1Bnd < 0 || v1Bnd > 1 || u1Bnd < 0 || u1Bnd > 1)
      {
        u1Bnd = 0.;
        v1Bnd = 0.;
      }
      if(!isfinite(u2Bnd) || u2Bnd < 0 || u2Bnd > 1)
      {
        u2Bnd = 0.;
        v2Bnd = 0.;
      }
      
      const double uMed = (u0Bnd+u1Bnd+u2Bnd)/2.;
      const double vMed = (v0Bnd+v1Bnd+v2Bnd)/2.;
      
      for(unsigned int k = 0; k < nbPoint; k++)
      {
        value += lp[nbPoint-1][k]*BF_order1(dim-1, i, lx[nbPoint-1][k])*BFE_order1(enrichment, dim, j, uMed + a*lx[nbPoint-1][k], vMed * b*lx[nbPoint-1][k])*std::abs(detJ);
      }
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

void getTotalNumbersOfNodes(GModel* m, Param param, Physical physical, int *nbTotNodes, std::map<int, std::vector<int> > *enrichedNodes)
{
  const int nbNodes = m->getNumMeshVertices();
  const double bnd = param.x_bnd;
  int totNodes = nbNodes;
  
  for(unsigned int i = 0; i < physical.elmOmega.size(); i++)
  {
    if(physical.elmOmega[i]->dim() == 1)
    {
      for(std::vector<MLine*>::iterator it = static_cast<GEdge*>(physical.elmOmega[i])->lines.begin(); it != static_cast<GEdge*>(physical.elmOmega[i])->lines.end(); ++it)
      {
        MLine* line = *it;
        const int tagElm = line->getNum();
        
        double x[2];
        x[0] = line->getVertex(0)->x();
        x[1] = line->getVertex(1)->x();
        
        if((bnd-x[0])*(bnd-x[1]) < 0)
        {
          std::vector<int> nodes;
          nodes.push_back(totNodes);
          nodes.push_back(totNodes+1);
          
          enrichedNodes->insert( std::pair<int, std::vector<int> >(tagElm,nodes) );
          totNodes += 2;
        }
      }
    }
    else if(physical.elmOmega[i]->dim() == 2)
    {
      for(std::vector<MTriangle*>::iterator it = static_cast<GFace*>(physical.elmOmega[i])->triangles.begin(); it != static_cast<GFace*>(physical.elmOmega[i])->triangles.end(); ++it)
      {
        MTriangle* triangle = *it;
        const int tagElm = triangle->getNum();
        
        double x[3];
        x[0] = triangle->getVertex(0)->x();
        x[1] = triangle->getVertex(1)->x();
        x[2] = triangle->getVertex(2)->x();
        double y[3];
        y[0] = triangle->getVertex(0)->y();
        y[1] = triangle->getVertex(1)->y();
        y[2] = triangle->getVertex(2)->y();
        double r[3];
        if(param.problem2D == Square)
        {
          r[0] = x[0];
          r[1] = x[1];
          r[2] = x[2];
        }
        else if(param.problem2D == Cylinder)
        {
          r[0] = std::sqrt(x[0]*x[0] + y[0]*y[0]);
          r[1] = std::sqrt(x[1]*x[1] + y[1]*y[1]);
          r[2] = std::sqrt(x[2]*x[2] + y[2]*y[2]);
        }
        
        gmm::dense_matrix<double> J = jacobianVol(triangle, 2);
        
        if((bnd-r[0])*(bnd-r[1]) < 0 || (bnd-r[1])*(bnd-r[2]) < 0 || (bnd-r[2])*(bnd-r[0]) < 0)
        {
          std::vector<int> nodes;
          nodes.push_back(totNodes);
          nodes.push_back(totNodes+1);
          nodes.push_back(totNodes+2);
          
          enrichedNodes->insert( std::pair<int, std::vector<int> >(tagElm,nodes) );
          totNodes += 3;
        }
      }
    }
  }
  
  *nbTotNodes = totNodes;
}

void getNumbersOfLagrangian(GModel* m, Param param, Physical physical, int nbNodes, std::map<int, std::vector<int> > &enrichedNodes, std::map<int, std::vector<int> > *lagrangianNodes, int *nbLag)
{
  if(m->getDim() == 1)
  {
    for(std::map<int, std::vector<int> >::iterator it = enrichedNodes.begin(); it != enrichedNodes.end(); ++it)
    {
      std::vector<int> nodes;
      nodes.push_back(nbNodes);
      
      lagrangianNodes->insert( std::pair<int, std::vector<int> >(it->first,nodes) );
      nbNodes += 1;
      (*nbLag) += 1;
    }
  }
  else if(m->getDim() == 2)
  {
    std::map<int, std::vector<int> > neighbor = getNeighbor(m, param, enrichedNodes);
    /*
    for(std::map<int, std::vector<int> >::iterator it = enrichedNodes.begin(); it != enrichedNodes.end(); ++it)
    {
      std::cout << it->first << std::endl;
      std::vector<int> node = neighbor[it->first];
      for(unsigned int i = 0; i < node.size(); i++)
      {
        std::cout << "\t" << node[i] << std::endl;
      }
    }*/
    
    for(std::map<int, std::vector<int> >::iterator it = enrichedNodes.begin(); it != enrichedNodes.end(); ++it)
    {
      std::vector<int> nodes;
      lagrangianNodes->insert( std::pair<int, std::vector<int> >(it->first,nodes) );
    }
    
    for(std::map<int, std::vector<int> >::iterator it = enrichedNodes.begin(); it != enrichedNodes.end(); ++it)
    {
      std::vector<int> neigh = neighbor[it->first];
      
      for(unsigned int i = 0; i < neigh.size(); i++)
      {
        int tagNeighbor = neighbor[it->first][i];
        if(tagNeighbor != -1)
        {
          if(tagNeighbor == -2)
          {
            lagrangianNodes->at(it->first).push_back(nbNodes);
            neighbor[it->first][i] = -1;
          }
          else
          {
            lagrangianNodes->at(it->first).push_back(nbNodes);
          
            neighbor[it->first][i] = -1;
          
            std::vector<int> other = neighbor[tagNeighbor];
            for(unsigned int j = 0; j < other.size(); j++)
            {
              if(other[j] == it->first)
              {
                neighbor[tagNeighbor][j] = -1;
                lagrangianNodes->at(tagNeighbor).push_back(nbNodes);
              }
            }
          }
          nbNodes ++;
          (*nbLag) ++;
        }
      }
    }
  }
  /*
  for(std::map<int, std::vector<int> >::iterator it = enrichedNodes.begin(); it != enrichedNodes.end(); ++it)
  {
    std::cout << it->first << std::endl;
    std::vector<int> node = lagrangianNodes->at(it->first);
    for(unsigned int i = 0; i < node.size(); i++)
    {
      std::cout << "\t" << node[i] << std::endl;
    }
  }*/
}

std::map<int, std::vector<int> > getNeighbor(GModel* m, Param param, std::map<int, std::vector<int> > &enrichedNodes)
{
  std::map<int, std::vector<int> > neighbor;
  
  for(std::map<int, std::vector<int> >::iterator it = enrichedNodes.begin(); it != enrichedNodes.end(); ++it)
  {
    MTriangle *triangle = static_cast<MTriangle*>(m->getMeshElementByTag(it->first));
    
    std::vector<int> vec;
    neighbor.insert( std::pair<int, std::vector<int> >(triangle->getNum(), vec) );
    
    MVertex* v1[3];
    v1[0] = triangle->getVertex(0);
    v1[1] = triangle->getVertex(1);
    v1[2] = triangle->getVertex(2);
    
    for(std::map<int, std::vector<int> >::iterator it2 = enrichedNodes.begin(); it2 != enrichedNodes.end(); ++it2)
    {
      MTriangle *triangle2 = static_cast<MTriangle*>(m->getMeshElementByTag(it2->first));
      
      if(triangle2 != triangle)
      {
        MVertex* v2[3];
        v2[0] = triangle2->getVertex(0);
        v2[1] = triangle2->getVertex(1);
        v2[2] = triangle2->getVertex(2);
        
        int commonVertices = 0;
      
        for(unsigned int i = 0; i < 3; i++)
        {
          for(unsigned int j = 0; j < 3; j++)
          {
            if(v1[i]->getNum() == v2[j]->getNum())
            {
              commonVertices++;
            }
          }
        }
        
        if(commonVertices == 2)
        {
          neighbor[triangle->getNum()].push_back(triangle2->getNum());
        }
      }
    }
    
    int commonVertices = 0;
    if(param.problem2D == Square)
    {
      for(unsigned int i = 0; i < 3; i++)
      {
        if(v1[i]->y() == 0. || v1[i]->y() == 1.)
        {
          commonVertices ++;
        }
      }
    }
    
    if(commonVertices == 2)
    {
      neighbor[triangle->getNum()].push_back(-2);
    }
  }
  
  return neighbor;
}


