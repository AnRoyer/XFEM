#include <complex>

#include "computationTools.h"
#include "gauss.h"

#include "GVertex.h"
#include "GEdge.h"

#include "MPoint.h"
#include "MLine.h"
#include "MTriangle.h"

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

double BF_order1(unsigned int dim, unsigned int num, double u, double v)
{
    double value = 0.;
    switch (dim) {
        case 0:
            /*
             *  phi_0 = 1;
             */
            switch (num) {
                case 0:
                    value = 1.;
                    break;
                default:
                    break;
            }
            break;
        case 1:
            /*
             *  phi_0 = (1-u)/2;
             *  phi_1 = (1+u)/2;
             *
             *  (0)                 (1)
             *  -|---------|---------|-->u
             *   -1        0         1
             */
            switch (num) {
                case 0:
                    value = (1-u)/2;
                    break;
                case 1:
                    value = (1+u)/2;
                    break;
                default:
                    break;
            }
            break;
        case 2:
            /*
             *  phi_0 = 1-u-v;
             *  phi_1 = u;
             *  phi_2 = v;
             *
             * v ^
             * 1 |- (2)
             *   | \__
             *   |    \__
             *   |       \__
             *   |          \__
             *   |             \__
             *   |(0)             \ (1)
             *  -|-----------------|-->u
             *   0                 1
             */
            switch (num) {
                case 0:
                    value = 1-u-v;
                    break;
                case 1:
                    value = u;
                    break;
                case 2:
                    value = v;
                    break;
                default:
                    break;
            }
            break;
        default:
            break;
    }
    
    return value;
}

std::vector<double> BFgrad_order1(unsigned int dim, unsigned int num, double u, double v)
{
    std::vector<double> value(2, 0.);
    switch (dim) {
        case 0:
            /*
             *  grad(phi_0) = 0;
             */
            switch (num) {
                case 0:
                    value.push_back(0.);
                    break;
                default:
                    break;
            }
            break;
        case 1:
            /*
             *  grad(phi_0) = -1/2;
             *  grad(phi_1) = 1/2;
             *
             *  -|---------|---------|-->u
             *   -1        0         1
             */
            switch (num) {
                case 0:
                    value[0] = -1./2;
                    break;
                case 1:
                    value[0] = 1./2;
                    break;
                default:
                    break;
            }
            break;
        case 2:
            /*
             *  phi_0 = 1-u-v;
             *  phi_1 = u;
             *  phi_2 = v;
             *
             * v ^
             * 1 |- (2)
             *   | \__
             *   |    \__
             *   |       \__
             *   |          \__
             *   |             \__
             *   |(0)             \ (1)
             *  -|-----------------|-->u
             *   0                 1
             */
            switch (num) {
                case 0:
                    value[0] = -1.;
                    value[1] = -1.;
                    break;
                case 1:
                    value[0] = 1.;
                    value[1] = 0.;
                    break;
                case 2:
                    value[0] = 0.;
                    value[1] = 1.;
                    break;
                default:
                    break;
            }
            break;
        default:
            break;
    }
    
    return value;
}

std::vector<double> BFEgrad_order1(unsigned int dim, unsigned int num, double u, double v)
{
    std::vector<double> value(2,0.);
    std::vector<double> BFgrad = BFgrad_order1(dim, num, u, v);
    std::vector<double> Fgrad = Fgrad_enrichment(dim, u, v);
    
    for(unsigned int i = 0; i < 2; i++)
    {
        value[i] = F_enrichment(dim, u, v)*BFgrad[i] + BF_order1(dim, num, u, v)*Fgrad[i];
    }
    
    return value;
}

double lsu[3];
void setLsu(int num, double value)
{
    lsu[num] = value;
}

//Enrichment function : sum_i(|ls_i|N_i(x)) - |sum_i(ls_i*N_i(x)|
double F_enrichment(unsigned int dim, double u, double v)
{
    double value = 0.;
    
    double term1 = 0.;
    double term2 = 0.;
    
    for(unsigned int i = 0; i < dim+1; i++)
    {
        term1 += std::abs(lsu[i])*BF_order1(dim, i, u);
        term2 += lsu[i]*BF_order1(dim, i, u);
    }
    
    value = term1 - std::abs(term2);
    
    return value;
}

std::vector<double> Fgrad_enrichment(unsigned int dim, double u, double v)
{
    std::vector<double> value(2, 0.);
    
    std::vector<double> term1(2, 0.);
    std::vector<double> term2(2, 0.);
    
    for(unsigned int i = 0; i < dim+1; i++)
    {
        std::vector<double> BFgrad = BFgrad_order1(dim, i, u);
        
        term1[0] += std::abs(lsu[i])*BFgrad[0];
        term1[1] += std::abs(lsu[i])*BFgrad[1];
        term2[0] += lsu[i]*BFgrad[0];
        term2[1] += lsu[i]*BFgrad[1];
    }
    
    for(unsigned int i = 0; i < dim; i++)
    {
        value[i] = term1[i] - std::abs(term2[i]);
    }
    
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

double integraleK(unsigned int dim, unsigned int i, unsigned int j, unsigned int nbPoint, gmm::dense_matrix<double> J)
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
            switch (nbPoint) {
                case 1:
                    for(unsigned int k = 0; k < 1; k++)
                    {
                        if(i >= 2 || j >= 2)
                        {
                            if(i >= 2 && j < 2)
                            {
                                std::vector<double>BFEgradI = prodMatVec(J, BFEgrad_order1(dim, i-2, lx1[k]));
                                std::vector<double>BFgradJ = prodMatVec(J, BFgrad_order1(dim, j, lx1[k]));
                                
                                value += lp1[k]*(BFEgradI[0]*BFgradJ[0]+BFEgradI[1]*BFgradJ[1])*std::abs(detJ);
                            }
                            else if(j >= 2 && i < 2)
                            {
                                std::vector<double>BFgradI = prodMatVec(J, BFgrad_order1(dim, i, lx1[k]));
                                std::vector<double>BFEgradJ = prodMatVec(J, BFEgrad_order1(dim, j-2, lx1[k]));
                                
                                value += lp1[k]*(BFgradI[0]*BFEgradJ[0]+BFgradI[1]*BFEgradJ[1])*std::abs(detJ);
                            }
                            else
                            {
                                std::vector<double>BFEgradI = prodMatVec(J, BFEgrad_order1(dim, i-2, lx1[k]));
                                std::vector<double>BFEgradJ = prodMatVec(J, BFEgrad_order1(dim, j-2, lx1[k]));
                                
                                value += lp1[k]*(BFEgradI[0]*BFEgradJ[0]+BFEgradI[1]*BFEgradJ[1])*std::abs(detJ);
                            }
                        }
                        else
                        {
                            std::vector<double>BFgradI = prodMatVec(J, BFgrad_order1(dim, i, lx1[k]));
                            std::vector<double>BFgradJ = prodMatVec(J, BFgrad_order1(dim, j, lx1[k]));
                            
                            value += lp1[k]*(BFgradI[0]*BFgradJ[0]+BFgradI[1]*BFgradJ[1])*std::abs(detJ);
                        }
                    }
                    break;
                case 2:
                    for(unsigned int k = 0; k < 2; k++)
                    {
                        if(i >= 2 || j >= 2)
                        {
                            if(i >= 2 && j < 2)
                            {
                                std::vector<double>BFEgradI = prodMatVec(J, BFEgrad_order1(dim, i-2, lx2[k]));
                                std::vector<double>BFgradJ = prodMatVec(J, BFgrad_order1(dim, j, lx2[k]));
                                
                                value += lp2[k]*(BFEgradI[0]*BFgradJ[0]+BFEgradI[1]*BFgradJ[1])*std::abs(detJ);
                            }
                            else if(j >= 2 && i < 2)
                            {
                                std::vector<double>BFgradI = prodMatVec(J, BFgrad_order1(dim, i, lx2[k]));
                                std::vector<double>BFEgradJ = prodMatVec(J, BFEgrad_order1(dim, j-2, lx2[k]));
                                
                                value += lp2[k]*(BFgradI[0]*BFEgradJ[0]+BFgradI[1]*BFEgradJ[1])*std::abs(detJ);
                            }
                            else
                            {
                                std::vector<double>BFEgradI = prodMatVec(J, BFEgrad_order1(dim, i-2, lx2[k]));
                                std::vector<double>BFEgradJ = prodMatVec(J, BFEgrad_order1(dim, j-2, lx2[k]));
                                
                                value += lp2[k]*(BFEgradI[0]*BFEgradJ[0]+BFEgradI[1]*BFEgradJ[1])*std::abs(detJ);
                            }
                        }
                        else
                        {
                            std::vector<double>BFgradI = prodMatVec(J, BFgrad_order1(dim, i, lx2[k]));
                            std::vector<double>BFgradJ = prodMatVec(J, BFgrad_order1(dim, j, lx2[k]));
                            
                            value += lp2[k]*(BFgradI[0]*BFgradJ[0]+BFgradI[1]*BFgradJ[1])*std::abs(detJ);
                        }
                    }
                    break;
                case 3:
                    for(unsigned int k = 0; k < 3; k++)
                    {
                        if(i >= 2 || j >= 2)
                        {
                            if(i >= 2 && j < 2)
                            {
                                std::vector<double>BFEgradI = prodMatVec(J, BFEgrad_order1(dim, i-2, lx3[k]));
                                std::vector<double>BFgradJ = prodMatVec(J, BFgrad_order1(dim, j, lx3[k]));
                                
                                value += lp3[k]*(BFEgradI[0]*BFgradJ[0]+BFEgradI[1]*BFgradJ[1])*std::abs(detJ);
                            }
                            else if(j >= 2 && i < 2)
                            {
                                std::vector<double>BFgradI = prodMatVec(J, BFgrad_order1(dim, i, lx3[k]));
                                std::vector<double>BFEgradJ = prodMatVec(J, BFEgrad_order1(dim, j-2, lx3[k]));
                                
                                value += lp3[k]*(BFgradI[0]*BFEgradJ[0]+BFgradI[1]*BFEgradJ[1])*std::abs(detJ);
                            }
                            else
                            {
                                std::vector<double>BFEgradI = prodMatVec(J, BFEgrad_order1(dim, i-2, lx3[k]));
                                std::vector<double>BFEgradJ = prodMatVec(J, BFEgrad_order1(dim, j-2, lx3[k]));
                                
                                value += lp3[k]*(BFEgradI[0]*BFEgradJ[0]+BFEgradI[1]*BFEgradJ[1])*std::abs(detJ);
                            }
                        }
                        else
                        {
                            std::vector<double>BFgradI = prodMatVec(J, BFgrad_order1(dim, i, lx3[k]));
                            std::vector<double>BFgradJ = prodMatVec(J, BFgrad_order1(dim, j, lx3[k]));
                            
                            value += lp3[k]*(BFgradI[0]*BFgradJ[0]+BFgradI[1]*BFgradJ[1])*std::abs(detJ);
                        }
                    }
                    break;
                case 4:
                    for(unsigned int k = 0; k < 4; k++)
                    {
                        if(i >= 2 || j >= 2)
                        {
                            if(i >= 2 && j < 2)
                            {
                                std::vector<double>BFEgradI = prodMatVec(J, BFEgrad_order1(dim, i-2, lx4[k]));
                                std::vector<double>BFgradJ = prodMatVec(J, BFgrad_order1(dim, j, lx4[k]));
                                
                                value += lp2[k]*(BFEgradI[0]*BFgradJ[0]+BFEgradI[1]*BFgradJ[1])*std::abs(detJ);
                            }
                            else if(j >= 2 && i < 2)
                            {
                                std::vector<double>BFgradI = prodMatVec(J, BFgrad_order1(dim, i, lx4[k]));
                                std::vector<double>BFEgradJ = prodMatVec(J, BFEgrad_order1(dim, j-2, lx4[k]));
                                
                                value += lp2[k]*(BFgradI[0]*BFEgradJ[0]+BFgradI[1]*BFEgradJ[1])*std::abs(detJ);
                            }
                            else
                            {
                                std::vector<double>BFEgradI = prodMatVec(J, BFEgrad_order1(dim, i-2, lx4[k]));
                                std::vector<double>BFEgradJ = prodMatVec(J, BFEgrad_order1(dim, j-2, lx4[k]));
                                
                                value += lp2[k]*(BFEgradI[0]*BFEgradJ[0]+BFEgradI[1]*BFEgradJ[1])*std::abs(detJ);
                            }
                        }
                        else
                        {
                            std::vector<double>BFgradI = prodMatVec(J, BFgrad_order1(dim, i, lx4[k]));
                            std::vector<double>BFgradJ = prodMatVec(J, BFgrad_order1(dim, j, lx4[k]));
                            
                            value += lp2[k]*(BFgradI[0]*BFgradJ[0]+BFgradI[1]*BFgradJ[1])*std::abs(detJ);
                        }
                    }
                    break;
                default:
                    break;
            }
            break;
        case 2:
            switch (nbPoint) {
                case 1:
                    for(unsigned int k = 0; k < 1; k++)
                    {
                        std::vector<double>BFgradI = prodMatVec(J, BFgrad_order1(dim, i, tx1[k], ty1[k]));
                        std::vector<double>BFgradJ = prodMatVec(J, BFgrad_order1(dim, j, tx1[k], ty1[k]));
                        
                        value += tp1[k]*(BFgradI[0]*BFgradJ[0]+BFgradI[1]*BFgradJ[1])*std::abs(detJ);
                    }
                    break;
                case 3:
                    for(unsigned int k = 0; k < 3; k++)
                    {
                        std::vector<double>BFgradI = prodMatVec(J, BFgrad_order1(dim, i, tx3[k], ty3[k]));
                        std::vector<double>BFgradJ = prodMatVec(J, BFgrad_order1(dim, j, tx3[k], ty3[k]));
                        
                        value += tp3[k]*(BFgradI[0]*BFgradJ[0]+BFgradI[1]*BFgradJ[1])*std::abs(detJ);
                    }
                    break;
                case 4:
                    for(unsigned int k = 0; k < 4; k++)
                    {
                        std::vector<double>BFgradI = prodMatVec(J, BFgrad_order1(dim, i, tx4[k], ty4[k]));
                        std::vector<double>BFgradJ = prodMatVec(J, BFgrad_order1(dim, j, tx4[k], ty4[k]));
                        
                        value += tp4[k]*(BFgradI[0]*BFgradJ[0]+BFgradI[1]*BFgradJ[1])*std::abs(detJ);
                    }
                    break;
                case 6:
                    for(unsigned int k = 0; k < 6; k++)
                    {
                        std::vector<double>BFgradI = prodMatVec(J, BFgrad_order1(dim, i, tx6[k], ty6[k]));
                        std::vector<double>BFgradJ = prodMatVec(J, BFgrad_order1(dim, j, tx6[k], ty6[k]));
                        
                        value += tp6[k]*(BFgradI[0]*BFgradJ[0]+BFgradI[1]*BFgradJ[1])*std::abs(detJ);
                    }
                    break;
                case 7:
                    for(unsigned int k = 0; k < 7; k++)
                    {
                        std::vector<double>BFgradI = prodMatVec(J, BFgrad_order1(dim, i, tx7[k], ty7[k]));
                        std::vector<double>BFgradJ = prodMatVec(J, BFgrad_order1(dim, j, tx7[k], ty7[k]));
                        
                        value += tp7[k]*(BFgradI[0]*BFgradJ[0]+BFgradI[1]*BFgradJ[1])*std::abs(detJ);
                    }
                    break;
                default:
                    break;
            }
            break;
        default:
            break;
    }
    
    return value;
}

double integraleM(unsigned int dim, unsigned int i, unsigned int j, unsigned int nbPoint, gmm::dense_matrix<double> J)
{
    double detJ = gmm::lu_det(J);
    
    double value = 0.;
    switch (dim) {
        case 0:
            value += pp1[0]*BF_order1(dim, i, px1[0])*BF_order1(dim, j, px1[0])*std::abs(detJ);
            break;
        case 1:
            switch (nbPoint) {
                case 1:
                    for(unsigned int k = 0; k < 1; k++)
                    {
                        if(i >= 2 || j >= 2)
                        {
                            if(i >= 2 && j < 2)
                            {
                                value += lp1[k]*BF_order1(dim, i-2, lx1[k])*F_enrichment(dim, lx1[k])*BF_order1(dim, j, lx1[k])*std::abs(detJ);
                            }
                            else if(j >= 2 && i < 2)
                            {
                                value += lp1[k]*BF_order1(dim, i, lx1[k])*BF_order1(dim, j-2, lx1[k])*F_enrichment(dim, lx1[k])*std::abs(detJ);
                            }
                            else
                            {
                                value += lp1[k]*BF_order1(dim, i-2, lx1[k])*F_enrichment(dim, lx1[k])*BF_order1(dim, j-2, lx1[k])*F_enrichment(dim, lx1[k])*std::abs(detJ);
                            }
                        }
                        else
                        {
                            value += lp1[k]*BF_order1(dim, i, lx1[k])*BF_order1(dim, j, lx1[k])*std::abs(detJ);
                        }
                    }
                    break;
                case 2:
                    for(unsigned int k = 0; k < 2; k++)
                    {
                        if(i >= 2 || j >= 2)
                        {
                            if(i >= 2 && j < 2)
                            {
                                value += lp2[k]*BF_order1(dim, i-2, lx2[k])*F_enrichment(dim, lx2[k])*BF_order1(dim, j, lx2[k])*std::abs(detJ);
                            }
                            else if(j >= 2 && i < 2)
                            {
                                value += lp2[k]*BF_order1(dim, i, lx2[k])*BF_order1(dim, j-2, lx2[k])*F_enrichment(dim, lx2[k])*std::abs(detJ);
                            }
                            else
                            {
                                value += lp2[k]*BF_order1(dim, i-2, lx2[k])*F_enrichment(dim, lx2[k])*BF_order1(dim, j-2, lx2[k])*F_enrichment(dim, lx2[k])*std::abs(detJ);
                            }
                        }
                        else
                        {
                            value += lp2[k]*BF_order1(dim, i, lx2[k])*BF_order1(dim, j, lx2[k])*std::abs(detJ);
                        }
                    }
                    break;
                case 3:
                    for(unsigned int k = 0; k < 3; k++)
                    {
                        if(i >= 2 || j >= 2)
                        {
                            if(i >= 2 && j < 2)
                            {
                                value += lp3[k]*BF_order1(dim, i-2, lx3[k])*F_enrichment(dim, lx3[k])*BF_order1(dim, j, lx3[k])*std::abs(detJ);
                            }
                            else if(j >= 2 && i < 2)
                            {
                                value += lp3[k]*BF_order1(dim, i, lx3[k])*BF_order1(dim, j-2, lx3[k])*F_enrichment(dim, lx3[k])*std::abs(detJ);
                            }
                            else
                            {
                                value += lp3[k]*BF_order1(dim, i-2, lx3[k])*F_enrichment(dim, lx3[k])*BF_order1(dim, j-2, lx3[k])*F_enrichment(dim, lx3[k])*std::abs(detJ);
                            }
                        }
                        else
                        {
                            value += lp3[k]*BF_order1(dim, i, lx3[k])*BF_order1(dim, j, lx3[k])*std::abs(detJ);
                        }
                    }
                    break;
                case 4:
                    for(unsigned int k = 0; k < 4; k++)
                    {
                        if(i >= 2 || j >= 2)
                        {
                            if(i >= 2 && j < 2)
                            {
                                value += lp4[k]*BF_order1(dim, i-2, lx4[k])*F_enrichment(dim, lx4[k])*BF_order1(dim, j, lx4[k])*std::abs(detJ);
                            }
                            else if(j >= 2 && i < 2)
                            {
                                value += lp4[k]*BF_order1(dim, i, lx4[k])*BF_order1(dim, j-2, lx4[k])*F_enrichment(dim, lx4[k])*std::abs(detJ);
                            }
                            else
                            {
                                value += lp4[k]*BF_order1(dim, i-2, lx4[k])*F_enrichment(dim, lx4[k])*BF_order1(dim, j-2, lx4[k])*F_enrichment(dim, lx4[k])*std::abs(detJ);
                            }
                        }
                        else
                        {
                            value += lp4[k]*BF_order1(dim, i, lx4[k])*BF_order1(dim, j, lx4[k])*std::abs(detJ);
                        }
                    }
                    break;
                default:
                    break;
            }
            break;
        case 2:
            switch (nbPoint) {
                case 1:
                    for(unsigned int k = 0; k < 1; k++)
                    {
                        value += tp1[k]*BF_order1(dim, i, tx1[k], ty1[k])*BF_order1(dim, j, tx1[k], ty1[k])*std::abs(detJ);
                    }
                    break;
                case 3:
                    for(unsigned int k = 0; k < 3; k++)
                    {
                        value += tp3[k]*BF_order1(dim, i, tx3[k], ty3[k])*BF_order1(dim, j, tx3[k], ty3[k])*std::abs(detJ);
                    }
                    break;
                case 4:
                    for(unsigned int k = 0; k < 4; k++)
                    {
                        value += tp4[k]*BF_order1(dim, i, tx4[k], ty4[k])*BF_order1(dim, j, tx4[k], ty4[k])*std::abs(detJ);
                    }
                    break;
                case 6:
                    for(unsigned int k = 0; k < 6; k++)
                    {
                        value += tp6[k]*BF_order1(dim, i, tx6[k], ty6[k])*BF_order1(dim, j, tx6[k], ty6[k])*std::abs(detJ);
                    }
                    break;
                case 7:
                    for(unsigned int k = 0; k < 7; k++)
                    {
                        value += tp7[k]*BF_order1(dim, i, tx7[k], ty7[k])*BF_order1(dim, j, tx7[k], ty7[k])*std::abs(detJ);
                    }
                    break;
                default:
                    break;
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
                
                Ktmp(nbN0, nbN0) += std::complex<double>(0., -k*integraleM(dim-1, 0, 0, 1, J));
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
                        Ke(i,j) = -k*integraleM(dim-1, i, j, 3, J);
                    }
                }
                
                /*std::cout << "Nodes : " << nbN0+1 << " ; " << nbN1+1 << std::endl;
                std::cout << Ke << std::endl;*/
                
                Ktmp(nbN0, nbN0) += std::complex<double>(0., Ke(0,0));
                Ktmp(nbN0, nbN1) += std::complex<double>(0., Ke(0,1));
                Ktmp(nbN1, nbN0) += std::complex<double>(0., Ke(1,0));
                Ktmp(nbN1, nbN1) += std::complex<double>(0., Ke(1,1));
            }
        }
    }
}
