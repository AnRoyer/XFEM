#include <iostream>

#include "fem.h"
#include "computationTools.h"

#include "MLine.h"
#include "MTriangle.h"

/*
 *  Helmholtz 1D :
 *
 *  Int_V{grad(u) grad(v)}dV - Int_V{k^2 u v}dV - Int_S{v grad(u) n}dS
 *= Int_V{grad(u) grad(v)}dV - Int_V{k^2 u v}dV - Int_S{i k u v}dS
 *
 */


std::vector< std::complex<double> > FEM::solve(GModel* m, int nbNodes, Param param, Physical physical)
{
    gmm::csr_matrix< std::complex<double> > K(nbNodes, nbNodes);
    gmm::row_matrix< gmm::wsvector< std::complex<double> > > Ktmp(nbNodes, nbNodes);
    std::vector< std::complex<double> > u(nbNodes);
    
    if(m->getDim() == 1)
    {
        computeK(Ktmp, m->firstEdge(), m->lastEdge(), nbNodes, param);
    }
    else if(m->getDim() == 2)
    {
        computeK(Ktmp, m->firstFace(), m->lastFace(), nbNodes, param);
    }
    
    gmm::dense_matrix< std::complex<double> > Kprint(nbNodes, nbNodes);
    gmm::copy(Ktmp,Kprint);
    std::cout << Kprint << std::endl;
    
    //q vector
    std::vector< std::complex<double> > q(nbNodes, std::complex<double>(0,0));
    
    //Boundaries
    for(unsigned int i = 0; i < physical.tagDir.size(); i++)
    {
        for(unsigned int j = 0; j < nbNodes; j++)
        {
            Ktmp(physical.tagDir[i]-1,j) = 0.;
        }
        Ktmp(physical.tagDir[i]-1, physical.tagDir[i]-1) = 1.;
        
        q[physical.tagDir[i]-1] = param.wave;
    }
    
    sommerfeldCondition(Ktmp, physical.elmInf, m->getDim(), param.k_1);
    
    gmm::copy(Ktmp,K);
    
    /*gmm::dense_matrix< std::complex<double> > Kprint(nbNodes, nbNodes);
    gmm::copy(Ktmp,Kprint);
    std::cout << Kprint << std::endl;*/

    
    //solver
    gmm::lu_solve(K, u, q);
    
    return u;
}

void FEM::computeK(gmm::row_matrix< gmm::wsvector< std::complex<double> > > &Ktmp, GModel::eiter eBegin, GModel::eiter eEnd, int nbNodes, Param param)
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
            
            gmm::dense_matrix<double> J = jacobianVol(line, 1);
            
            for(unsigned int i = 0; i < 2; i++)
            {
                for(unsigned int j = 0; j < 2; j++)
                {
                    Ke(i,j) = integraleK(1, i, j, 3, J);
                    if(x0 < param.x_bnd && x1 > param.x_bnd)
                    {
                        if(i == 0)
                        {
                            Me(i,j) = - param.k_0*param.k_0*integraleM(1, i, j, 3, J);
                        }
                        else if(i == 1)
                        {
                            Me(i,j) = - param.k_1*param.k_1*integraleM(1, i, j, 3, J);
                        }
                    }
                    else if(x0 < param.x_bnd && x1 <= param.x_bnd)
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

void FEM::computeK(gmm::row_matrix< gmm::wsvector< std::complex<double> > > &Ktmp, GModel::fiter fBegin, GModel::fiter fEnd, int nbNodes, Param param)
{
    for(GModel::fiter it = fBegin; it != fEnd; ++it)
    {
        for(std::vector<MTriangle*>::iterator itT = (*it)->triangles.begin(); itT != (*it)->triangles.end(); ++itT)
        {
            MTriangle* triangle = *itT;
            int nbN0 = triangle->getVertex(0)->getNum()-1;
            int nbN1 = triangle->getVertex(1)->getNum()-1;
            int nbN2 = triangle->getVertex(2)->getNum()-1;
            
            double x[3];
            x[0] = triangle->getVertex(0)->x();
            x[1] = triangle->getVertex(1)->x();
            x[2] = triangle->getVertex(2)->x();
            
            gmm::dense_matrix<double> Ke(3,3);
            gmm::dense_matrix<double> Me(3,3);
            
            gmm::dense_matrix<double> J = jacobianVol(triangle, 2);
            
            for(unsigned int i = 0; i < 3; i++)
            {
                for(unsigned int j = 0; j < 3; j++)
                {
                    Ke(i,j) = integraleK(2, i, j, 7, J);
                    if(x[i] <= param.x_bnd)
                    {
                        Me(i,j) = - param.k_0*param.k_0*integraleM(2, i, j, 7, J);
                    }
                    else
                    {
                        Me(i,j) = - param.k_1*param.k_1*integraleM(2, i, j, 7, J);
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





