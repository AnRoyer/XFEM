#include <cmath>
#include <complex>

#include "analytical.h"
#include "MPoint.h"
#include "MLine.h"
#include "bessel.h"

using namespace std;

std::vector< std::complex<double> > ANALYTICAL::solve(GModel* m, Param param, Physical physical, bool xfem)
{
  int nbNodes = m->getNumMeshVertices();
  vector< complex<double> > u(nbNodes);
  
  const double k1 = param.k_1;
  const double k2 = param.k_2;
  double a = 0.;
  if(xfem == false)
  {
    if(m->getDim() == 1)
    {
      a = static_cast<GVertex*>(physical.elmInter[0])->points[0]->getVertex(0)->x();
    }
    else if(m->getDim() == 2)
    {
      a = static_cast<GEdge*>(physical.elmInter[0])->lines[0]->getVertex(0)->x();
    }
  }
  else
  {
    a = param.x_bnd;
  }
  complex<double> U = param.wave;
  
  std::cout << "k1 = " << k1 << std::endl;
  std::cout << "k2 = " << k2 << std::endl;
  
  if(m->getDim() == 1)
  {
    complex<double> e_ik1 = exp(complex<double>(0., -k1*a));
    complex<double> e_2ik1 = exp(complex<double>(0., -2*k1*a));
    complex<double> e_ik2 = exp(complex<double>(0., -k2*a));
  
    //Z_ac = rho * c
    const double Z2 = param.rho_1*param.c_1;
    const double Z1 = param.rho_2*param.c_2;
  
  
    complex<double> A = U/((Z1-Z2)/(Z1+Z2) * e_2ik1 + 1.);
  
    complex<double> B = U - A;
  
    complex<double> C = 0.;
  
    complex<double> D = 2*Z1/(Z1+Z2) * A * e_ik1 / e_ik2;
  
    for(unsigned int i = 0; i < nbNodes; i++)
    {
      const double x = m->getMeshVertexByTag(i+1)->x();
      if(x <= a)
      {
        u[i] = A*exp(complex<double>(0., -k1*x)) + B*exp(complex<double>(0., k1*x));
      }
      else
      {
        u[i] = C*exp(complex<double>(0., k2*x)) + D*exp(complex<double>(0., -k2*x));
      }
    }
  }
  else if(m->getDim() == 2)
  {
    if(param.problem2D == Square)
    {
      complex<double> e_ik1 = exp(complex<double>(0., -k1*a));
      complex<double> e_2ik1 = exp(complex<double>(0., -2*k1*a));
      complex<double> e_ik2 = exp(complex<double>(0., -k2*a));
      
      //Z_ac = rho * c
      const double Z2 = param.rho_1*param.c_1;
      const double Z1 = param.rho_2*param.c_2;
      
      
      complex<double> A = U/((Z1-Z2)/(Z1+Z2) * e_2ik1 + 1.);
      
      complex<double> B = U - A;
      
      complex<double> C = 0.;
      
      complex<double> D = 2*Z1/(Z1+Z2) * A * e_ik1 / e_ik2;
      
      for(unsigned int i = 0; i < nbNodes; i++)
      {
        const double x = m->getMeshVertexByTag(i+1)->x();
        if(x <= a)
        {
          u[i] = A*exp(complex<double>(0., -k1*x)) + B*exp(complex<double>(0., k1*x));
        }
        else
        {
          u[i] = C*exp(complex<double>(0., k2*x)) + D*exp(complex<double>(0., -k2*x));
        }
      }
    }
  }
  
  return u;
}
