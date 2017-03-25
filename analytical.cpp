#include <cmath>

#include "analytical.h"

using namespace std;

std::vector< std::complex<double> > ANALYTICAL::solve(GModel* m, int nbNodes, Param param, Physical physical)
{
    vector< complex<double> > u(nbNodes);
    
    const double k1 = param.k_0;
    const double k2 = param.k_1;
    const double a = param.x_bnd;
    const double L = m->getMeshVertexByTag(physical.tagInf[0])->x();
    complex<double> U = param.wave;
    
    const double Ar = -(k2*k2*sin(2*a*k1)*U.real() - k1*k1*sin(2*a*k1)*U.real() + 2*k1*k2*U.imag())/(k1*k1*cos(2*a*k1) - k2*k2*cos(2*a*k1) + k1*k1 + k2*k2);
    const double Ai = (k1*k1*U.imag()*sin(2*a*k1) - k2*k2*U.imag()*sin(2*a*k1) + 2*k1*k2*U.real())/(k1*k1*cos(2*a*k1) - k2*k2*cos(2*a*k1) + k1*k1 + k2*k2);
    
    complex<double> A(Ar, Ai);
    
    complex<double> B(U);

    const double Cr = -(k1*(k1*cos(a*k1)*cos(a*k2)*U.imag() - k1*cos(a*k1)*sin(a*k2)*U.real() + k2*cos(a*k2)*sin(a*k1)*U.real() + k2*U.imag()*sin(a*k1)*sin(a*k2)))/(k1*k1*cos(a*k1)*cos(a*k1) - k2*k2*cos(a*k1)*cos(a*k1) + k2*k2);
    const double Ci = (k1*(k1*cos(a*k1)*cos(a*k2)*U.real() + k1*cos(a*k1)*U.imag()*sin(a*k2) - k2*cos(a*k2)*U.imag()*sin(a*k1) + k2*sin(a*k1)*sin(a*k2)*U.real()))/(k1*k1*cos(a*k1)*cos(a*k1) - k2*k2*cos(a*k1)*cos(a*k1) + k2*k2);
    
    complex<double> C(Cr, Ci);
    
    const double Dr =  (k1*(k1*cos(a*k1)*cos(a*k2)*U.real() + k1*cos(a*k1)*U.imag()*sin(a*k2) - k2*cos(a*k2)*U.imag()*sin(a*k1) + k2*sin(a*k1)*sin(a*k2)*U.real()))/(k1*k1*cos(a*k1)*cos(a*k1) - k2*k2*cos(a*k1)*cos(a*k1) + k2*k2);
    const double Di = (k1*(k1*cos(a*k1)*cos(a*k2)*U.imag() - k1*cos(a*k1)*sin(a*k2)*U.real() + k2*cos(a*k2)*sin(a*k1)*U.real() + k2*U.imag()*sin(a*k1)*sin(a*k2)))/(k1*k1*cos(a*k1)*cos(a*k1) - k2*k2*cos(a*k1)*cos(a*k1) + k2*k2);
    
    complex<double> D(Dr, Di);
    
    for(unsigned int i = 0; i < nbNodes; i++)
    {
        double x = m->getMeshVertexByTag(i+1)->x();
        if(x <= a)
        {
            u[i] = A*sin(k1*x) + B*cos(k1*x);
        }
        else
        {
            u[i] = C*sin(k2*x) + D*cos(k2*x);
        }
    }
    
    return u;
}
