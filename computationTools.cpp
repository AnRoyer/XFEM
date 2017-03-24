#include "computationTools.h"
#include "gauss.h"

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
             *  phi_0 = (1-x)/2;
             *  phi_1 = (1+x)/2;
             *
             *  -|---------|---------|-->x
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
        default:
            break;
    }
    
    return value;
}

double BFgrad_order1(unsigned int dim, unsigned int num, double u, double v)
{
    double value = 0.;
    switch (dim) {
        case 0:
            /*
             *  grad(phi_0) = 0;
             */
            switch (num) {
                case 0:
                    value = 0.;
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
             *  -|---------|---------|-->x
             *   -1        0         1
             */
            switch (num) {
                case 0:
                    value = -1./2;
                    break;
                case 1:
                    value = 1./2;
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

double jacobian1D(MElement* e)
{
    double x0 = e->getVertex(0)->x();
    double x1 = e->getVertex(1)->x();
    
    double J = (x1-x0)/2;
    
    return J;
}

double integraleK(unsigned int dim, unsigned int i, unsigned int j, unsigned int nbPoint)
{
    double value = 0.;
    switch (dim) {
        case 0:
            break;
        case 1:
            switch (nbPoint) {
                case 1:
                    for(unsigned int k = 0; k < 1; k++)
                    {
                        value += lp1[k]*BFgrad_order1(dim, i, lx1[k])*BFgrad_order1(dim, j, lx1[k]);
                    }
                    break;
                case 2:
                    for(unsigned int k = 0; k < 2; k++)
                    {
                        value += lp2[k]*BFgrad_order1(dim, i, lx2[k])*BFgrad_order1(dim, j, lx2[k]);
                    }
                    break;
                case 3:
                    for(unsigned int k = 0; k < 3; k++)
                    {
                        value += lp3[k]*BFgrad_order1(dim, i, lx3[k])*BFgrad_order1(dim, j, lx3[k]);
                    }
                    break;
                case 4:
                    for(unsigned int k = 0; k < 4; k++)
                    {
                        value += lp4[k]*BFgrad_order1(dim, i, lx4[k])*BFgrad_order1(dim, j, lx4[k]);
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

double integraleM(unsigned int dim, unsigned int i, unsigned int j, unsigned int nbPoint)
{
    double value = 0.;
    switch (dim) {
        case 0:
            break;
        case 1:
            switch (nbPoint) {
                case 1:
                    for(unsigned int k = 0; k < 1; k++)
                    {
                        value += lp1[k]*BF_order1(dim, i, lx1[k])*BF_order1(dim, j, lx1[k]);
                    }
                    break;
                case 2:
                    for(unsigned int k = 0; k < 2; k++)
                    {
                        value += lp2[k]*BF_order1(dim, i, lx2[k])*BF_order1(dim, j, lx2[k]);
                    }
                    break;
                case 3:
                    for(unsigned int k = 0; k < 3; k++)
                    {
                        value += lp3[k]*BF_order1(dim, i, lx3[k])*BF_order1(dim, j, lx3[k]);
                    }
                    break;
                case 4:
                    for(unsigned int k = 0; k < 4; k++)
                    {
                        value += lp4[k]*BF_order1(dim, i, lx4[k])*BF_order1(dim, j, lx4[k]);
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
