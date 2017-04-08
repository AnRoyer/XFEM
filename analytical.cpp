#include <cmath>
#include <complex>

#include "analytical.h"
#include "MPoint.h"

using namespace std;

std::vector< std::complex<double> > ANALYTICAL::solve(GModel* m, int nbNodes, Param param, Physical physical)
{
    vector< complex<double> > u(nbNodes);
    /*
    const double k1 = param.k_0;
    const double k2 = param.k_1;
    const double a = param.x_bnd;
    const double L = m->getMeshVertexByTag(static_cast<GVertex*>(physical.elmInf[0])->points[0]->getVertex(0)->getNum()-1)->x();
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
    */
    
    const double k1 = param.k_0;
    const double k2 = param.k_1;
    const double Z1 = 1./param.k_0;
    const double Z2 = 1./param.k_1;
    const double a = param.x_bnd;
    complex<double> U = param.wave;
    
    complex<double> e_ik1 = exp(complex<double>(0., -k1*a));
    complex<double> e_2ik1 = exp(complex<double>(0., -2*k1*a));
    complex<double> e_ik2 = exp(complex<double>(0., -k2*a));
    
    
    complex<double> A = U/((Z2-Z1)/(Z1+Z2) * e_2ik1 + 1.);
    
    complex<double> B = U - A;
    
    complex<double> C = 0.;
    
    complex<double> D = 2*Z2/(Z1+Z2) * A * e_ik1 / e_ik2;
    
    for(unsigned int i = 0; i < nbNodes; i++)
    {
        double x = m->getMeshVertexByTag(i+1)->x();
        if(x <= a)
        {
            u[i] = A*exp(complex<double>(0., -k1*x)) + B*exp(complex<double>(0., k1*x));
        }
        else
        {
            u[i] = C*exp(complex<double>(0., k2*x)) + D*exp(complex<double>(0., -k2*x));
        }
    }
    
    //Let's take the conjugate to have a wave that propagate in the same direction.
    for(unsigned int i = 0; i < nbNodes; i++)
    {
        u[i] = conj(u[i]);
    }
    
    return u;
}
