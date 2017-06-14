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
    gmm::csr_matrix< std::complex<double> > K(nbNodes+2, nbNodes+2);
    gmm::row_matrix< gmm::wsvector< std::complex<double> > > Ktmp(nbNodes+2, nbNodes+2);
    std::vector< std::complex<double> > u(nbNodes+2);
    
    //K matrix
    computeK(Ktmp, physical.elmOmega, param.k_1, param.rho_1, param.c_1, param.k_2, param.rho_2, param.c_2, param.x_bnd, nbNodes);
    
    //q vector
    std::vector< std::complex<double> > q(nbNodes+2, std::complex<double>(0,0));
    
    //Boundary conditions
    sommerfeldCondition(Ktmp, physical.elmInf, m->getDim(), param.rho_2*param.c_2*param.c_2*param.k_2);
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
    
    //std::cout << Ktmp << std::endl;
    
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

void XFEM::computeK(gmm::row_matrix< gmm::wsvector< std::complex<double> > > &Ktmp, std::vector<GEntity*> elms, double k1, double rho1, double c1, double k2, double rho2, double c2, double bnd, int nbNodes)
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
                
                gmm::dense_matrix<double> J = jacobianVol(line, 1);
                
                if((bnd-x[0])*(bnd-x[1]) < 0)
                {
                    gmm::dense_matrix<double> Ke(4,4);
                    gmm::dense_matrix<double> Me(4,4);
                    
                    setLsu(0, (bnd-x[0])/J(0,0));
                    setLsu(1, (bnd-x[1])/J(0,0));
                    
                    for(unsigned int i = 0; i < 4; i++)//i = 2, i = 3 -> enrichment
                    {
                        for(unsigned int j = 0; j < 4; j++)//j = 2, j = 3 -> enrichment
                        {
                            //Better to use an even number of intergration points to avoid the point x=0 where the deriative is not defined
                            if(i == 0 || i == 2)
                            {
                                Ke(i,j) = rho1*c1*c1*integraleK(1, i, j, 20, J);
                                Me(i,j) = - rho1*c1*c1*k1*k1*integraleM(1, i, j, 20, J);
                            }
                            else if(i == 1 || i == 3)
                            {
                                Ke(i,j) = rho2*c2*c2*integraleK(1, i, j, 20, J);
                                Me(i,j) = - rho2*c2*c2*k2*k2*integraleM(1, i, j, 20, J);
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
                    gmm::dense_matrix<double> Ke(2,2);
                    gmm::dense_matrix<double> Me(2,2);
                    
                    for(unsigned int i = 0; i < 2; i++)
                    {
                        for(unsigned int j = 0; j < 2; j++)
                        {
                            if(bnd <= x[0] && bnd <= x[1])
                            {
                                Ke(i,j) = rho1*c1*c1*integraleK(1, i, j, 3, J);
                                Me(i,j) = - rho1*c1*c1*k1*k1*integraleM(1, i, j, 3, J);
                            }
                            else if(bnd >= x[0] && bnd >= x[1])
                            {
                                Ke(i,j) = rho2*c2*c2*integraleK(1, i, j, 3, J);
                                Me(i,j) = - rho2*c2*c2*k2*k2*integraleM(1, i, j, 3, J);
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
