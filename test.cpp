#include <iostream>
#include "test.h"
#include "computationTools.h"
#include "struct.h"

void testH(int dim)
{
  if(dim == 1)
  {
    for(unsigned int i = 0; i < 4; i++)
    {
      std::cout << "Basic function " << i << ":" << std::endl;
      setLst(0, -1.);
      setLst(1, 1.);
      
      std::cout << "\t{";
      for(unsigned int j = 0; j < 10; j++)
      {
        double bfe = BFE_order1(Heaviside, dim, i, 2.*j/9.-1., 0);
        std::cout << bfe << ", ";
      }
      std::cout << "}\n\t";
      for(unsigned int j = 0; j < 10; j++)
      {
        double bfe = BFE_order1(Heaviside, dim, i, 2.*j/9.-1., 0);
        if(bfe == 0)
        {
          std::cout << "0";
        }
        else if(bfe > 0)
        {
          std::cout << "+";
        }
      }
      
      std::cout << "\n\n" << std::endl;
    }
  }
  if(dim == 2)
  {
    for(unsigned int i = 0; i < 6; i++)
    {
      std::cout << "Basic function " << i << ":" << std::endl;
      setLst(0, 1.);
      setLst(1, 1.);
      setLst(2, -1.);
      
      std::cout << "\t";
      for(unsigned int j = 0; j < 5; j++)
      {
        for(unsigned int k = 0; k < 2*j+1; k++)
        {
          double bfe = BFE_order1(Heaviside, dim, i, (double)k/10., -(double)j/4.+1.);
          if(bfe == 0)
          {
            std::cout << "0";
          }
          else if(bfe > 0)
          {
            std::cout << "+";
          }
        }
        std::cout << std::endl;
        std::cout << "\t";
      }
      
      std::cout << "\n\n" << std::endl;
    }
  }
}
