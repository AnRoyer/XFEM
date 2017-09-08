#include <iostream>
#include <cmath>
#include "gmm_iter_solvers.h"

#include "fem.h"
#include "computationTools.h"

#include "MLine.h"
#include "MTriangle.h"

std::vector< std::complex<double> > FEM::solve(GModel* m, Param param, Physical physical)
{
  int nbNodes = m->getNumMeshVertices();
  gmm::csr_matrix< std::complex<double> > K(nbNodes, nbNodes);
  gmm::row_matrix< gmm::wsvector< std::complex<double> > > Ktmp(nbNodes, nbNodes);
  std::vector< std::complex<double> > u(nbNodes);
  
  //K matrix
  computeK(Ktmp, physical.elmOmega1, param.k_1, param.rho_1, param.c_1);
  computeK(Ktmp, physical.elmOmega2, param.k_2, param.rho_2, param.c_2);
  
  //q vector
  std::vector< std::complex<double> > q(nbNodes, std::complex<double>(0,0));
  
  //Boundary conditions
  sommerfeldCondition(Ktmp, physical.elmInf, m->getDim(), 1./param.rho_2*param.k_2);
  for(unsigned int i = 0; i < physical.elmDir.size(); i++)
  {
    for(unsigned int j = 0; j < physical.elmDir[i]->mesh_vertices.size() ; j++)
    {
      const int tag = physical.elmDir[i]->mesh_vertices[j]->getNum()-1;
      for(unsigned int k = 0; k < nbNodes; k++)
      {
        Ktmp(tag,k) = 0.;
      }
      Ktmp(tag, tag) = 1.;
      q[tag] = param.wave;
      //q[tag] = std::complex<double>(cos(param.k_2*physical.elmDir[i]->mesh_vertices[j]->x()), -sin(param.k_2*physical.elmDir[i]->mesh_vertices[j]->x()));
    }
  }
  
  //solver
  gmm::copy(Ktmp,K);
  gmm::lu_solve(K, u, q);
  
  return u;
}

void FEM::computeK(gmm::row_matrix< gmm::wsvector< std::complex<double> > > &Ktmp, std::vector<GEntity*> elms, double k, double rho, double c)
{
  for(unsigned int i = 0; i < elms.size(); i++)
  {
    if(elms[i]->dim() == 1)
    {
      for(std::vector<MLine*>::iterator it = static_cast<GEdge*>(elms[i])->lines.begin(); it != static_cast<GEdge*>(elms[i])->lines.end(); ++it)
      {
        MLine* line = *it;
        int nbN0 = line->getVertex(0)->getNum()-1;
        int nbN1 = line->getVertex(1)->getNum()-1;
      
        double x[2];
        x[0] = line->getVertex(0)->x();
        x[1] = line->getVertex(1)->x();
      
        gmm::dense_matrix<double> Ke(2,2);
        gmm::dense_matrix<double> Me(2,2);
      
        gmm::dense_matrix<double> J = jacobianVol(line, 1);
      
        for(unsigned int i = 0; i < 2; i++)
        {
          for(unsigned int j = 0; j < 2; j++)
          {
            Ke(i,j) = 1./rho*integraleK(1, i, j, 3, J);
            Me(i,j) = - 1./rho*k*k*integraleM(1, i, j, 3, J);
          }
        }

        Ktmp(nbN0, nbN0) += Ke(0,0) + Me(0,0);
        Ktmp(nbN0, nbN1) += Ke(0,1) + Me(0,1);
        Ktmp(nbN1, nbN0) += Ke(1,0) + Me(1,0);
        Ktmp(nbN1, nbN1) += Ke(1,1) + Me(1,1);
      }
    }
    else if(elms[i]->dim() == 2)
    {
      for(std::vector<MTriangle*>::iterator it = static_cast<GFace*>(elms[i])->triangles.begin(); it != static_cast<GFace*>(elms[i])->triangles.end(); ++it)
      {
        MTriangle* triangle = *it;
        int nbN0 = triangle->getVertex(0)->getNum()-1;
        int nbN1 = triangle->getVertex(1)->getNum()-1;
        int nbN2 = triangle->getVertex(2)->getNum()-1;
        
        double x[3];
        x[0] = triangle->getVertex(0)->x();
        x[1] = triangle->getVertex(1)->x();
        x[2] = triangle->getVertex(2)->x();
        
        double y[3];
        y[0] = triangle->getVertex(0)->y();
        y[1] = triangle->getVertex(1)->y();
        y[2] = triangle->getVertex(2)->y();
        
        gmm::dense_matrix<double> Ke(3,3);
        gmm::dense_matrix<double> Me(3,3);
        
        gmm::dense_matrix<double> J = jacobianVol(triangle, 2);
        
        for(unsigned int i = 0; i < 3; i++)
        {
          for(unsigned int j = 0; j < 3; j++)
          {
            Ke(i,j) = 1./rho*integraleK(2, i, j, 7, J);
            Me(i,j) = - 1./rho*k*k*integraleM(2, i, j, 7, J);
          }
        }
        
        Ktmp(nbN0, nbN0) += Ke(0,0) + Me(0,0);
        Ktmp(nbN0, nbN1) += Ke(0,1) + Me(0,1);
        Ktmp(nbN0, nbN2) += Ke(0,2) + Me(0,2);
        
        Ktmp(nbN1, nbN0) += Ke(1,0) + Me(1,0);
        Ktmp(nbN1, nbN1) += Ke(1,1) + Me(1,1);
        Ktmp(nbN1, nbN2) += Ke(1,2) + Me(1,2);
        
        Ktmp(nbN2, nbN0) += Ke(2,0) + Me(2,0);
        Ktmp(nbN2, nbN1) += Ke(2,1) + Me(2,1);
        Ktmp(nbN2, nbN2) += Ke(2,2) + Me(2,2);
      }
    }
  }
}

