#include <iostream>
#include <complex>

#include "xfem.h"
#include "computationTools.h"

#include <MLine.h>

std::vector< std::complex<double> > XFEM::solve(GModel* m, Param param, Physical physical)
{
  int nbNodes = m->getNumMeshVertices();
  gmm::csr_matrix< std::complex<double> > K(nbNodes+2, nbNodes+2);
  gmm::row_matrix< gmm::wsvector< std::complex<double> > > Ktmp(nbNodes+2, nbNodes+2);
  std::vector< std::complex<double> > u(nbNodes+2);
  
  //K matrix
  computeK(Ktmp, physical.elmOmega, param.k_1, param.rho_1, param.c_1, param.k_2, param.rho_2, param.c_2, param.x_bnd, nbNodes);
  
  //q vector
  std::vector< std::complex<double> > q(nbNodes+2, std::complex<double>(0,0));
  
  //Boundary conditions
  sommerfeldCondition(Ktmp, physical.elmInf, m->getDim(), 1./param.rho_2*param.k_2);
  for(unsigned int i = 0; i < physical.elmDir.size(); i++)
  {
    for(unsigned int j = 0; j < physical.elmDir[i]->mesh_vertices.size() ; j++)
    {
      const int tag = physical.elmDir[i]->mesh_vertices[j]->getNum()-1;
      for(unsigned int k = 0; k < nbNodes+2; k++)
      {
        Ktmp(tag,k) = 0.;
      }
      Ktmp(tag, tag) = 1.;
      q[tag] = param.wave;
      //q[tag] = std::complex<double>(cos(param.k_2*physical.elmDir[i]->mesh_vertices[j]->x()), -sin(param.k_2*physical.elmDir[i]->mesh_vertices[j]->x()));
    }
  }
  /*std::complex<double> I(0.,0.1);
  q[nbNodes] = I;
  q[nbNodes+1] = I;*/
  
  //solver
  gmm::copy(Ktmp,K);
  //std::cout << Ktmp << std::endl;
  gmm::lu_solve(K, u, q);
  
  for(unsigned int i = 0; i < 2; i++)
  {
    std::cout << "enrichment = " << u[nbNodes+i] << std::endl;
  }
  
  std::vector< std::complex<double> > uRed(nbNodes);
  for(unsigned int i = 0; i < nbNodes; i++)
  {
    uRed[i] = u[i];
  }
  
  return uRed;
}

void XFEM::computeK(gmm::row_matrix< gmm::wsvector< std::complex<double> > > &Ktmp, std::vector<GEntity*> elms, double k1, double rho1, double c1, double k2, double rho2, double c2, double bnd, int nbNodes)
{
  //Ktmp(nbNodes, nbNodes) += 1;
  //Ktmp(nbNodes+1, nbNodes+1) += 1;
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
        
        gmm::dense_matrix<double> J = jacobianVol(line, 1);
        
        if((bnd-x[0])*(bnd-x[1]) < 0)
        {
          gmm::dense_matrix<double> Ke(4,4);
          gmm::dense_matrix<double> Me(4,4);
          
          //lst[0] must be smaller than lst[1]
          setLst(0, (x[0]-bnd)/J(0,0) );
          setLst(1, (x[1]-bnd)/J(0,0) );
          
          for(unsigned int i = 0; i < 4; i++)//i = 2, i = 3 -> enrichment
          {
            for(unsigned int j = 0; j < 4; j++)//j = 2, j = 3 -> enrichment
            {
              if(x[i%2] < bnd)
              {
                Ke(i,j) = 1./rho1*integraleK(1, i, j, 20, J, Hat);
                Me(i,j) = - 1./rho1*k1*k1*integraleM(1, i, j, 20, J, Hat);
              }
              else if(x[i%2] > bnd)
              {
                Ke(i,j) = 1./rho2*integraleK(1, i, j, 20, J, Hat);
                Me(i,j) = - 1./rho2*k2*k2*integraleM(1, i, j, 20, J, Hat);
              }
              
            }
          }
          
          Ktmp(nbN0, nbN0) += Ke(0,0) + Me(0,0);
          Ktmp(nbN0, nbN1) += Ke(0,1) + Me(0,1);
          Ktmp(nbN0, nbNodes) += Ke(0,2) + Me(0,2);
          Ktmp(nbN0, nbNodes+1) += Ke(0,3) + Me(0,3);
          
          Ktmp(nbN1, nbN0) += Ke(1,0) + Me(1,0);
          Ktmp(nbN1, nbN1) += Ke(1,1) + Me(1,1);
          Ktmp(nbN1, nbNodes) += Ke(1,2) + Me(1,2);
          Ktmp(nbN1, nbNodes+1) += Ke(1,3) + Me(1,3);
          
          Ktmp(nbNodes, nbN0) += Ke(2,0) + Me(2,0);
          Ktmp(nbNodes, nbN1) += Ke(2,1) + Me(2,1);
          Ktmp(nbNodes, nbNodes) += Ke(2,2) + Me(2,2);
          Ktmp(nbNodes, nbNodes+1) += Ke(2,3) + Me(2,3);
          
          Ktmp(nbNodes+1, nbN0) += Ke(3,0) + Me(3,0);
          Ktmp(nbNodes+1, nbN1) += Ke(3,1) + Me(3,1);
          Ktmp(nbNodes+1, nbNodes) += Ke(3,2) + Me(3,2);
          Ktmp(nbNodes+1, nbNodes+1) += Ke(3,3) + Me(3,3);
        }
        else
        {
          gmm::dense_matrix<double> Ke(2,2);
          gmm::dense_matrix<double> Me(2,2);
          
          for(unsigned int i = 0; i < 2; i++)
          {
            for(unsigned int j = 0; j < 2; j++)
            {
              if(bnd >= x[0] && bnd >= x[1])
              {
                Ke(i,j) = 1./rho1*integraleK(1, i, j, 3, J);
                Me(i,j) = - 1./rho1*k1*k1*integraleM(1, i, j, 3, J);
              }
              else if(bnd <= x[0] && bnd <= x[1])
              {
                Ke(i,j) = 1./rho2*integraleK(1, i, j, 3, J);
                Me(i,j) = - 1./rho2*k2*k2*integraleM(1, i, j, 3, J);
              }
            }
          }
        
          Ktmp(nbN0, nbN0) += Ke(0,0) + Me(0,0);
          Ktmp(nbN0, nbN1) += Ke(0,1) + Me(0,1);
          Ktmp(nbN1, nbN0) += Ke(1,0) + Me(1,0);
          Ktmp(nbN1, nbN1) += Ke(1,1) + Me(1,1);
        }
      }
    }
  }
}
