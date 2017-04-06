#include <iostream>

#include "xfem.h"
#include "computationTools.h"

#include <MLine.h>

/*
 *  Helmholtz 1D :
 *
 *  Int_V{grad(u) grad(v)}dV - Int_V{k^2 u v}dV - Int_S{v grad(u)}dS
 *
 *
 */


std::vector< std::complex<double> > XFEM::solve(GModel* m, int nbNodes, Param param, Physical physical)
{
    //The last componant (nbNodes) correspond to the dof that is added by XFEM
    gmm::csr_matrix< std::complex<double> > K(nbNodes+2, nbNodes+2);
    gmm::row_matrix< gmm::wsvector< std::complex<double> > > Ktmp(nbNodes+2, nbNodes+2);
    std::vector< std::complex<double> > u(nbNodes+2);
    
    computeK(Ktmp, m->firstEdge(), m->lastEdge(), nbNodes, param);
    
    //q vector
    std::vector< std::complex<double> > q(nbNodes+2, std::complex<double>(0,0));
    
    //Boundaries
    for(unsigned int i = 0; i < physical.tagDir.size(); i++)
    {
        for(unsigned int j = 0; j < nbNodes; j++)
        {
            Ktmp(physical.tagDir[i]-1,j) = 0;
        }
        Ktmp(physical.tagDir[i]-1, physical.tagDir[i]-1) = 1;
        
        q[physical.tagDir[i]-1] = param.wave;
    }
    
    sommerfeldCondition(Ktmp, physical.elmInf, m->getDim(), param.k_1);
    
    gmm::copy(Ktmp,K);
    
    gmm::dense_matrix< std::complex<double> > Kprint(nbNodes+2, nbNodes+2);
    gmm::copy(Ktmp,Kprint);
    std::cout << Kprint << std::endl;
    
    //solver
    gmm::lu_solve(K, u, q);
    
    std::vector< std::complex<double> > uRed(nbNodes);
    for(unsigned int i = 0; i < nbNodes; i++)
    {
        uRed[i] = u[i];
    }
    
    
    return uRed;
}

void XFEM::computeK(gmm::row_matrix< gmm::wsvector< std::complex<double> > > &Ktmp, GModel::eiter eBegin, GModel::eiter eEnd, int nbNodes, Param param)
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
            
            gmm::dense_matrix<double> J = jacobianVol(line, 1);
            
            if(x0 < param.x_bnd && x1 > param.x_bnd)
            {
                //Real XFEM case
                gmm::dense_matrix<double> Ke(4,4);
                gmm::dense_matrix<double> Me(4,4);
                
                setLsu(0, (param.x_bnd-x0)/J(0,0));
                setLsu(1, (param.x_bnd-x1)/J(0,0));
                
                for(unsigned int i = 0; i < 4; i++)//i = 2, i = 3 -> enrichment
                {
                    for(unsigned int j = 0; j < 4; j++)//j = 2, j = 3 -> enrichment
                    {
                        Ke(i,j) = integraleK(1, i, j, 3, J);
                        if(i == 0 || i == 2)
                        {
                            Me(i,j) = - param.k_0*param.k_0*integraleM(1, i, j, 3, J);
                        }
                        else if(i == 1 || i == 3)
                        {
                            Me(i,j) = - param.k_1*param.k_1*integraleM(1, i, j, 3, J);
                        }
                    }
                }
                
                Ktmp(nbN0, nbN0) += Ke(0,0) + Me(0,0);
                Ktmp(nbN0, nbN1) += Ke(0,1) + Me(0,1);
                Ktmp(nbN1, nbN0) += Ke(1,0) + Me(1,0);
                Ktmp(nbN1, nbN1) += Ke(1,1) + Me(1,1);
                
                Ktmp(nbNodes, nbN0) += Ke(2,0) + Me(2,0);
                Ktmp(nbNodes, nbN1) += Ke(2,1) + Me(2,1);
                Ktmp(nbN0, nbNodes) += Ke(0,2) + Me(0,2);
                Ktmp(nbN1, nbNodes) += Ke(1,2) + Me(1,2);
                
                Ktmp(nbNodes+1, nbN0) += Ke(3,0) + Me(3,0);
                Ktmp(nbNodes+1, nbN1) += Ke(3,1) + Me(3,1);
                Ktmp(nbN0, nbNodes+1) += Ke(0,3) + Me(0,3);
                Ktmp(nbN1, nbNodes+1) += Ke(1,3) + Me(1,3);
                
                Ktmp(nbNodes+1, nbNodes) += Ke(3,2) + Me(3,2);
                Ktmp(nbNodes+1, nbNodes+1) += Ke(3,3) + Me(3,3);
                Ktmp(nbNodes, nbNodes+1) += Ke(2,3) + Me(2,3);
                Ktmp(nbNodes, nbNodes) += Ke(2,2) + Me(2,2);
            }
            else
            {
                //Classical FEM
                gmm::dense_matrix<double> Ke(2,2);
                gmm::dense_matrix<double> Me(2,2);
                
                for(unsigned int i = 0; i < 2; i++)
                {
                    for(unsigned int j = 0; j < 2; j++)
                    {
                        Ke(i,j) = integraleK(1, i, j, 3, J);
                        if(x0 < param.x_bnd && x1 <= param.x_bnd)
                        {
                            Me(i,j) = - param.k_0*param.k_0*integraleM(1, i, j, 3, J);
                        }
                        else if(x0 >= param.x_bnd && x1 > param.x_bnd)
                        {
                            Me(i,j) = - param.k_1*param.k_1*integraleM(1, i, j, 3, J);
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

