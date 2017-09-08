#include <iostream>
#include <complex>
#include <fstream>
#include <cmath>

#include "lagrange.h"
#include "computationTools.h"

#include <MLine.h>
#include <MTriangle.h>

std::vector< std::complex<double> > LAGRANGE::solve(GModel* m, Param param, Physical physical)
{
  int nbNodes = m->getNumMeshVertices();
  int nbTotNodes = 0;
  int nbLag = 0;
  int dim = m->getDim();
  std::map<int, std::vector<int> > enrichedNodes;
  std::map<int, std::vector<int> > lagrangianNodes;
  
  getTotalNumbersOfNodes(m, param, physical, &nbTotNodes, &enrichedNodes);
  getNumbersOfLagrangian(m, param, physical, nbTotNodes, enrichedNodes, &lagrangianNodes, &nbLag);
  
  std::cout << nbTotNodes << " dof and " << nbLag << " lagrangian." << std::endl;
  
  /*for(std::map<int, std::vector<int> >::iterator it = enrichedNodes.begin(); it != enrichedNodes.end(); ++it)
  {
    std::cout << "Element " << it->first << " has enriched nodes:" << std::endl;
    std::cout << "\t" << it->second[0] << std::endl;
    std::cout << "\t" << it->second[1] << std::endl;
    std::cout << "\t" << it->second[2] << std::endl;
  }*/
  
  gmm::csr_matrix< std::complex<double> > K(nbTotNodes+nbLag, nbTotNodes+nbLag);
  gmm::row_matrix< gmm::wsvector< std::complex<double> > > Ktmp(nbTotNodes+nbLag, nbTotNodes+nbLag);
  std::vector< std::complex<double> > u(nbTotNodes+nbLag);
  
  //K matrix
  computeK(Ktmp, physical.elmOmega, param, enrichedNodes, lagrangianNodes);
  
  
  /*
  gmm::dense_matrix< std::complex<double> > M(nbTotNodes+nbLag, nbTotNodes+nbLag);
  gmm::copy(Ktmp,M);
  //std::cout << M << std::endl;
  
  std::ofstream Km("Km.txt");
  Km << "K = [";
  for(unsigned int i = 0; i < nbTotNodes+nbLag; i++)
  {
    for(unsigned int j = 0; j < nbTotNodes+nbLag; j++)
    {
      Km << M(i,j).real();
      if(j != nbTotNodes+nbLag-1) Km << ", ";
    }
    if(i != nbTotNodes+nbLag-1) Km << "; ";
  }
  Km << "]; ";
  Km.close();
  */
  
  //q vector
  std::vector< std::complex<double> > q(nbTotNodes+nbLag, std::complex<double>(0,0));
  
  //Boundary conditions
  sommerfeldCondition(Ktmp, physical.elmInf, m->getDim(), 1./param.rho_2*param.k_2);
  for(unsigned int i = 0; i < physical.elmDir.size(); i++)
  {
    for(unsigned int j = 0; j < physical.elmDir[i]->mesh_vertices.size() ; j++)
    {
      const int tag = physical.elmDir[i]->mesh_vertices[j]->getNum()-1;
      for(unsigned int k = 0; k < nbTotNodes+nbLag; k++)
      {
        Ktmp(tag,k) = 0.;
        //Ktmp(nbNodes,k) = 0.;
      }
      Ktmp(tag, tag) = 1.;
      //Ktmp(nbNodes,nbNodes) = 1.;
      //Ktmp(nbNodes,nbNodes+1) = -1.;
      q[tag] = param.wave;
    }
  }
  
  //solver
  gmm::copy(Ktmp,K);
  gmm::lu_solve(K, u, q);
  
  std::vector< std::complex<double> > uRed(nbNodes);
  for(unsigned int i = 0; i < nbNodes; i++)
  {
    uRed[i] = u[i];
  }
  
  return uRed;
}

void LAGRANGE::computeK(gmm::row_matrix< gmm::wsvector< std::complex<double> > > &Ktmp, std::vector<GEntity*> elms, Param param, std::map<int, std::vector<int> > &enrichedNodes, std::map<int, std::vector<int> > &lagrangianNodes)
{
  const double k1 = param.k_1;
  const double rho1 = param.rho_1;
  const double c1 = param.c_1;
  const double k2 = param.k_2;
  const double rho2 = param.rho_2;
  const double c2 = param.c_2;
  const double bnd = param.x_bnd;
  
  for(unsigned int i = 0; i < elms.size(); i++)
  {
    if(elms[i]->dim() == 1)
    {
      for(std::vector<MLine*>::iterator it = static_cast<GEdge*>(elms[i])->lines.begin(); it != static_cast<GEdge*>(elms[i])->lines.end(); ++it)
      {
        MLine* line = *it;
        const int nbN0 = line->getVertex(0)->getNum()-1;
        const int nbN1 = line->getVertex(1)->getNum()-1;
        
        double x[2];
        x[0] = line->getVertex(0)->x();
        x[1] = line->getVertex(1)->x();
        
        gmm::dense_matrix<double> J = jacobianVol(line, 1);
        
        if((bnd-x[0])*(bnd-x[1]) < 0)
        {
          const int tagElm = line->getNum();
                    
          const int nbN0E = enrichedNodes[tagElm][0];
          const int nbN1E = enrichedNodes[tagElm][1];
          
          const int nbN0Lag = lagrangianNodes[tagElm][0];
          
          gmm::dense_matrix<double> Ke(4,4);
          gmm::dense_matrix<double> Me(4,4);
          
          const double levelSet[2] = {x[0]-bnd, x[1]-bnd};
          setLst(0, levelSet[0]);
          setLst(1, levelSet[1]);
          
          for(unsigned int i = 0; i < 4; i++)
          {
            for(unsigned int j = 0; j < 4; j++)
            {
              if(levelSet[i%2] < 0)
              {
                if(i/2 == 0)
                {
                  Ke(i,j) = 1./rho1*integraleK(1, i, j, 3, J, Heaviside);
                  Me(i,j) = - 1./rho1*k1*k1*integraleM(1, i, j, 3, J, Heaviside);
                }
                else
                {
                  Ke(i,j) = 1./rho2*integraleK(1, i, j, 3, J, Heaviside);
                  Me(i,j) = - 1./rho2*k2*k2*integraleM(1, i, j, 3, J, Heaviside);
                }
              }
              else if(levelSet[i%2] > 0)
              {
                if(i/2 == 0)
                {
                  Ke(i,j) = 1./rho2*integraleK(1, i, j, 3, J, Heaviside);
                  Me(i,j) = - 1./rho2*k2*k2*integraleM(1, i, j, 3, J, Heaviside);
                }
                else
                {
                  Ke(i,j) = 1./rho1*integraleK(1, i, j, 3, J, Heaviside);
                  Me(i,j) = - 1./rho1*k1*k1*integraleM(1, i, j, 3, J, Heaviside);
                }
              }
            }
          }
          
          Ktmp(nbN0, nbN0)   += Ke(0,0) + Me(0,0);
          Ktmp(nbN0, nbN1)   += Ke(0,1) + Me(0,1);
          Ktmp(nbN0, nbN0E)  += Ke(0,2) + Me(0,2);
          Ktmp(nbN0, nbN1E)  += Ke(0,3) + Me(0,3);
          
          Ktmp(nbN1, nbN0)   += Ke(1,0) + Me(1,0);
          Ktmp(nbN1, nbN1)   += Ke(1,1) + Me(1,1);
          Ktmp(nbN1, nbN0E)  += Ke(1,2) + Me(1,2);
          Ktmp(nbN1, nbN1E)  += Ke(1,3) + Me(1,3);
          
          Ktmp(nbN0E, nbN0)  += Ke(2,0) + Me(2,0);
          Ktmp(nbN0E, nbN1)  += Ke(2,1) + Me(2,1);
          Ktmp(nbN0E, nbN0E) += Ke(2,2) + Me(2,2);
          Ktmp(nbN0E, nbN1E) += Ke(2,3) + Me(2,3);
          
          Ktmp(nbN1E, nbN0)  += Ke(3,0) + Me(3,0);
          Ktmp(nbN1E, nbN1)  += Ke(3,1) + Me(3,1);
          Ktmp(nbN1E, nbN0E) += Ke(3,2) + Me(3,2);
          Ktmp(nbN1E, nbN1E) += Ke(3,3) + Me(3,3);
          
          // u_1 = u_2
          Ktmp(nbN0, nbN0Lag)  += (levelSet[0] < 0 ? 1. : -1.)*integraleLD(1, 0, 0, 1, J, Heaviside);
          Ktmp(nbN1, nbN0Lag)  += (levelSet[1] < 0 ? 1. : -1.)*integraleLD(1, 0, 1, 1, J, Heaviside);
          Ktmp(nbN0E, nbN0Lag) += (levelSet[0] > 0 ? 1. : -1.)*integraleLD(1, 0, 2, 1, J, Heaviside);
          Ktmp(nbN1E, nbN0Lag) += (levelSet[1] > 0 ? 1. : -1.)*integraleLD(1, 0, 3, 1, J, Heaviside);
          
          Ktmp(nbN0Lag, nbN0)  += Ktmp(nbN0, nbN0Lag);
          Ktmp(nbN0Lag, nbN1)  += Ktmp(nbN1, nbN0Lag);
          Ktmp(nbN0Lag, nbN0E) += Ktmp(nbN0E, nbN0Lag);
          Ktmp(nbN0Lag, nbN1E) += Ktmp(nbN1E, nbN0Lag);
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
    else if(elms[i]->dim() == 2)
    {
      for(std::vector<MTriangle*>::iterator it = static_cast<GFace*>(elms[i])->triangles.begin(); it != static_cast<GFace*>(elms[i])->triangles.end(); ++it)
      {
        MTriangle* triangle = *it;
        const int nbN0 = triangle->getVertex(0)->getNum()-1;
        const int nbN1 = triangle->getVertex(1)->getNum()-1;
        const int nbN2 = triangle->getVertex(2)->getNum()-1;
        
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
          const int tagElm = triangle->getNum();
          
          const int nbN0E = enrichedNodes[tagElm][0];
          const int nbN1E = enrichedNodes[tagElm][1];
          const int nbN2E = enrichedNodes[tagElm][2];
          
          const int nbN0Lag = lagrangianNodes[tagElm][0];
          const int nbN1Lag = lagrangianNodes[tagElm][1];
          
          gmm::dense_matrix<double> Ke(6,6);
          gmm::dense_matrix<double> Me(6,6);
          
          std::vector<double> v0{x[0]-bnd, y[0]};
          std::vector<double> v1{x[1]-bnd, y[1]};
          std::vector<double> v2{x[2]-bnd, y[2]};

          const double levelSet[3] = {x[0]-bnd, x[1]-bnd, x[2]-bnd};
          const int order[6] = {0,1,2,(levelSet[0]*levelSet[1] < 0 ? 1 : 2),(levelSet[1]*levelSet[2] < 0 ? 2 : 0),(levelSet[2]*levelSet[0] < 0 ? 0 : 1)};
          setLst(0, levelSet[0]);
          setLst(1, levelSet[1]);
          setLst(2, levelSet[2]);
          
          for(unsigned int i = 0; i < 6; i++)
          {
            for(unsigned int j = 0; j < 6; j++)
            {
              if(levelSet[i%3] < 0)
              {
                if(i/3 == 0)
                {
                  Ke(i,j) = 1./rho1*integraleK(2, i, j, 3, J, Heaviside);
                  Me(i,j) = - 1./rho1*k1*k1*integraleM(2, i, j, 3, J, Heaviside);
                }
                else
                {
                  Ke(i,j) = 1./rho2*integraleK(2, i, j, 3, J, Heaviside);
                  Me(i,j) = - 1./rho2*k2*k2*integraleM(2, i, j, 3, J, Heaviside);
                }
              }
              else if(levelSet[i%3] > 0)
              {
                if(i/3 == 0)
                {
                  Ke(i,j) = 1./rho2*integraleK(2, i, j, 3, J, Heaviside);
                  Me(i,j) = - 1./rho2*k2*k2*integraleM(2, i, j, 3, J, Heaviside);
                }
                else
                {
                  Ke(i,j) = 1./rho1*integraleK(2, i, j, 3, J, Heaviside);
                  Me(i,j) = - 1./rho1*k1*k1*integraleM(2, i, j, 3, J, Heaviside);
                }
              }
            }
          }
          
          Ktmp(nbN0, nbN0)  += Ke(0,0) + Me(0,0);
          Ktmp(nbN0, nbN1)  += Ke(0,1) + Me(0,1);
          Ktmp(nbN0, nbN2)  += Ke(0,2) + Me(0,2);
          Ktmp(nbN0, nbN0E) += Ke(0,3) + Me(0,3);
          Ktmp(nbN0, nbN1E) += Ke(0,4) + Me(0,4);
          Ktmp(nbN0, nbN2E) += Ke(0,5) + Me(0,5);
          
          Ktmp(nbN1, nbN0)  += Ke(1,0) + Me(1,0);
          Ktmp(nbN1, nbN1)  += Ke(1,1) + Me(1,1);
          Ktmp(nbN1, nbN2)  += Ke(1,2) + Me(1,2);
          Ktmp(nbN1, nbN0E) += Ke(1,3) + Me(1,3);
          Ktmp(nbN1, nbN1E) += Ke(1,4) + Me(1,4);
          Ktmp(nbN1, nbN2E) += Ke(1,5) + Me(1,5);
          
          Ktmp(nbN2, nbN0)  += Ke(2,0) + Me(2,0);
          Ktmp(nbN2, nbN1)  += Ke(2,1) + Me(2,1);
          Ktmp(nbN2, nbN2)  += Ke(2,2) + Me(2,2);
          Ktmp(nbN2, nbN0E) += Ke(2,3) + Me(2,3);
          Ktmp(nbN2, nbN1E) += Ke(2,4) + Me(2,4);
          Ktmp(nbN2, nbN2E) += Ke(2,5) + Me(2,5);
          
          Ktmp(nbN0E, nbN0)  += Ke(3,0) + Me(3,0);
          Ktmp(nbN0E, nbN1)  += Ke(3,1) + Me(3,1);
          Ktmp(nbN0E, nbN2)  += Ke(3,2) + Me(3,2);
          Ktmp(nbN0E, nbN0E) += Ke(3,3) + Me(3,3);
          Ktmp(nbN0E, nbN1E) += Ke(3,4) + Me(3,4);
          Ktmp(nbN0E, nbN2E) += Ke(3,5) + Me(3,5);
          
          Ktmp(nbN1E, nbN0)  += Ke(4,0) + Me(4,0);
          Ktmp(nbN1E, nbN1)  += Ke(4,1) + Me(4,1);
          Ktmp(nbN1E, nbN2)  += Ke(4,2) + Me(4,2);
          Ktmp(nbN1E, nbN0E) += Ke(4,3) + Me(4,3);
          Ktmp(nbN1E, nbN1E) += Ke(4,4) + Me(4,4);
          Ktmp(nbN1E, nbN2E) += Ke(4,5) + Me(4,5);
          
          Ktmp(nbN2E, nbN0)  += Ke(5,0) + Me(5,0);
          Ktmp(nbN2E, nbN1)  += Ke(5,1) + Me(5,1);
          Ktmp(nbN2E, nbN2)  += Ke(5,2) + Me(5,2);
          Ktmp(nbN2E, nbN0E) += Ke(5,3) + Me(5,3);
          Ktmp(nbN2E, nbN1E) += Ke(5,4) + Me(5,4);
          Ktmp(nbN2E, nbN2E) += Ke(5,5) + Me(5,5);
          
          
          // u_1 = u_2
          Ktmp(nbN0, nbN0Lag)  += (levelSet[0] < 0 ? 1. : -1.)*integraleLD(2, 0, 0, 3, J, Heaviside);
          Ktmp(nbN1, nbN0Lag)  += (levelSet[1] < 0 ? 1. : -1.)*integraleLD(2, 0, 1, 3, J, Heaviside);
          Ktmp(nbN2, nbN0Lag)  += (levelSet[2] < 0 ? 1. : -1.)*integraleLD(2, 0, 2, 3, J, Heaviside);
          Ktmp(nbN0E, nbN0Lag) += (levelSet[0] > 0 ? 1. : -1.)*integraleLD(2, 0, 3, 3, J, Heaviside);
          Ktmp(nbN1E, nbN0Lag) += (levelSet[1] > 0 ? 1. : -1.)*integraleLD(2, 0, 4, 3, J, Heaviside);
          Ktmp(nbN2E, nbN0Lag) += (levelSet[2] > 0 ? 1. : -1.)*integraleLD(2, 0, 5, 3, J, Heaviside);
          
          Ktmp(nbN0Lag, nbN0)  += Ktmp(nbN0, nbN0Lag);
          Ktmp(nbN0Lag, nbN1)  += Ktmp(nbN1, nbN0Lag);
          Ktmp(nbN0Lag, nbN2)  += Ktmp(nbN2, nbN0Lag);
          Ktmp(nbN0Lag, nbN0E) += Ktmp(nbN0E, nbN0Lag);
          Ktmp(nbN0Lag, nbN1E) += Ktmp(nbN1E, nbN0Lag);
          Ktmp(nbN0Lag, nbN2E) += Ktmp(nbN2E, nbN0Lag);
          
          Ktmp(nbN0, nbN1Lag)  += (levelSet[order[0]] < 0 ? 1. : -1.)*integraleLD(2, 1, 0, 3, J, Heaviside);
          Ktmp(nbN1, nbN1Lag)  += (levelSet[order[1]] < 0 ? 1. : -1.)*integraleLD(2, 1, 1, 3, J, Heaviside);
          Ktmp(nbN2, nbN1Lag)  += (levelSet[order[2]] < 0 ? 1. : -1.)*integraleLD(2, 1, 2, 3, J, Heaviside);
          Ktmp(nbN0E, nbN1Lag) += (levelSet[order[0]] > 0 ? 1. : -1.)*integraleLD(2, 1, 3, 3, J, Heaviside);
          Ktmp(nbN1E, nbN1Lag) += (levelSet[order[1]] > 0 ? 1. : -1.)*integraleLD(2, 1, 4, 3, J, Heaviside);
          Ktmp(nbN2E, nbN1Lag) += (levelSet[order[2]] > 0 ? 1. : -1.)*integraleLD(2, 1, 5, 3, J, Heaviside);
          
          Ktmp(nbN1Lag, nbN0)  += Ktmp(nbN0, nbN1Lag);
          Ktmp(nbN1Lag, nbN1)  += Ktmp(nbN1, nbN1Lag);
          Ktmp(nbN1Lag, nbN2)  += Ktmp(nbN2, nbN1Lag);
          Ktmp(nbN1Lag, nbN0E) += Ktmp(nbN0E, nbN1Lag);
          Ktmp(nbN1Lag, nbN1E) += Ktmp(nbN1E, nbN1Lag);
          Ktmp(nbN1Lag, nbN2E) += Ktmp(nbN2E, nbN1Lag);
          
          //Stabilized factor
          gmm::dense_matrix<double> Jl = jacobianLagrange(1, J);
          
          std::cout << integraleK(1, 0, 0, 3, Jl) << std::endl;
          Ktmp(nbN0Lag, nbN0Lag) += -9e-7*integraleK(1, 0, 0, 3, Jl);
          Ktmp(nbN0Lag, nbN1Lag) += -9e-7*integraleK(1, 0, 1, 3, Jl);
          Ktmp(nbN1Lag, nbN0Lag) += -9e-7*integraleK(1, 1, 0, 3, Jl);
          Ktmp(nbN1Lag, nbN1Lag) += -9e-7*integraleK(1, 1, 1, 3, Jl);
        }
        else
        {
          gmm::dense_matrix<double> Ke(3,3);
          gmm::dense_matrix<double> Me(3,3);
          
          for(unsigned int i = 0; i < 3; i++)
          {
            for(unsigned int j = 0; j < 3; j++)
            {
              if(bnd >= r[0] && bnd >= r[1] && bnd >= r[2])
              {
                Ke(i,j) = 1./rho1*integraleK(2, i, j, 3, J);
                Me(i,j) = - 1./rho1*k1*k1*integraleM(2, i, j, 3, J);
              }
              else if(bnd <= r[0] && bnd <= r[1] && bnd <= r[2])
              {
                Ke(i,j) = 1./rho2*integraleK(2, i, j, 3, J);
                Me(i,j) = - 1./rho2*k2*k2*integraleM(2, i, j, 3, J);
              }
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
}
