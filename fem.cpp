#include <iostream>

#include "fem.h"
#include "computationTools.h"

/*
 *  Helmholtz 1D :
 *
 *  Int_V{grad(u) grad(v)}dV - Int_V{k^2 u v}dV - Int_S{v grad(u)}dS
 *
 *
 */


std::vector< std::complex<double> > solveFEM(GModel* m, int nbNodes, Param param, Physical physical)
{
    gmm::csr_matrix< std::complex<double> > K(nbNodes, nbNodes);
    gmm::row_matrix< gmm::wsvector< std::complex<double> > > Ktmp(nbNodes, nbNodes);
    std::vector< std::complex<double> > u(nbNodes);
    
    computeK(Ktmp, m->firstEdge(), m->lastEdge(), nbNodes, param);
    
    //q vector
    std::vector< std::complex<double> > q(nbNodes, std::complex<double>(0,0));
    
    //Boundaries
    for(unsigned int i = 0; i < physical.tagDir.size(); i++)
    {
        std::cout << "Dirichlet boundary conditions." << std::endl;
        for(unsigned int j = 0; j < nbNodes; j++)
        {
            Ktmp(physical.tagDir[i]-1,j) = 0;
        }
        Ktmp(physical.tagDir[i]-1, physical.tagDir[i]-1) = 1;
        
        q[physical.tagDir[i]-1] = param.wave;
    }
    
    for(unsigned int i = 0; i < physical.tagInf.size(); i++)
    {
        std::cout << "Sommerfeld radiation condition." << std::endl;
        Ktmp(physical.tagInf[i]-1, physical.tagInf[i]-1) += std::complex<double>(0.0, -param.k_1);//Sommerfeld radiation condition
    }
    
    gmm::copy(Ktmp,K);
    
    /*for(unsigned int i = 0; i < nbNodes; i++)
    {
         for(unsigned int j = 0; j < nbNodes; j++)
        {
            std::cout << K(i,j) << " ";
        }
        std::cout << std::endl;
    }*/
    
    //solver
    gmm::lu_solve(K, u, q);
    
    return u;
}

void computeK(gmm::row_matrix< gmm::wsvector< std::complex<double> > > &Ktmp, GModel::eiter eBegin, GModel::eiter eEnd, int nbNodes, Param param)
{
    for(GModel::eiter it = eBegin; it != eEnd; ++it)
    {
        for(std::vector<MLine*>::iterator itL = (*it)->lines.begin(); itL != (*it)->lines.end(); ++itL)
        {
            MLine* line = *itL;
            int nbN0 = line->getVertex(0)->getNum()-1;
            int nbN1 = line->getVertex(1)->getNum()-1;
            
            double x0 = line->getVertex(0)->x();
            double x1 = line->getVertex(1)->x();
            
            gmm::dense_matrix<double> Ke(2,2);
            gmm::dense_matrix<double> Me(2,2);
            
            double J = jacobian1D(line);
            
            for(unsigned int i = 0; i < 2; i++)
            {
                for(unsigned int j = 0; j < 2; j++)
                {
                    Ke(i,j) = integraleK(1, i, j, 3)/J;
                    if(x0 < param.x_bnd && x1 > param.x_bnd)
                    {
                        Me(i,j) = - (param.k_0+param.k_1)/2*(param.k_0+param.k_1)/2*integraleM(1, i, j, 3)*J;
                    }
                    else if(x0 < param.x_bnd && x1 <= param.x_bnd)
                    {
                        Me(i,j) = - param.k_0*param.k_0*integraleM(1, i, j, 3)*J;
                    }
                    else if(x0 >= param.x_bnd && x1 > param.x_bnd)
                    {
                        Me(i,j) = - param.k_1*param.k_1*integraleM(1, i, j, 3)*J;
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





