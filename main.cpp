#include <iostream>
#include <string>
#include <fstream>
#include <complex>
#include <cmath>

#include "Gmsh.h"
#include "GModel.h"
#include "MElement.h"
#include "MPoint.h"
#include "MLine.h"
#include "MTriangle.h"
#include "MVertex.h"

#include "analytical.h"
#include "fem.h"
#include "xfem.h"
#include "lagrange.h"
#include "error.h"

#include "computationTools.h"

#include "test.h"

#define GAMMADIR 1
#define GAMMAINF 2
#define OMEGA1 3
#define OMEGA2 4
#define INTERFACE 5
#define OMEGA 6

void writePOS(GModel* m, std::vector< std::complex<double> > u, std::string name);
void writeL2Error(std::vector< std::complex<double> > u, std::vector< std::complex<double> > e);
Param readParam(int argc, char **argv);
Physical checkPhysical(GModel *m);

int main(int argc, char **argv)
{
  std::cout << "#################################################" << std::endl;
  std::cout << "#                                               #" << std::endl;
  std::cout << "#  #       #  ########  ########  #      #      #" << std::endl;
  std::cout << "#   #     #   #         #         ##    ##      #" << std::endl;
  std::cout << "#    #   #    #         #         # #  # #      #" << std::endl;
  std::cout << "#     # #     #         #         #  ##  #      #" << std::endl;
  std::cout << "#      #      #####     #####     #      #      #" << std::endl;
  std::cout << "#     # #     #         #         #      #      #" << std::endl;
  std::cout << "#    #   #    #         #         #      #      #" << std::endl;
  std::cout << "#   #     #   #         #         #      #      #" << std::endl;
  std::cout << "#  #       #  #         ########  #      #      #" << std::endl;
  std::cout << "#                                               #" << std::endl;
  std::cout << "#################################################" << std::endl << std::endl;
  
  if(argc < 1)
  {
    std::cout << "Missing a mesh file!" << std::endl;
    return 0;
  }
  
  GmshInitialize(1, argv);
  GModel *m = new GModel();
  std::cout << "Reading msh file... " << std::flush;
  m->readMSH(argv[1]);
  std::cout << "Done!" << std::endl;
  
  Param param = readParam(argc, argv);
  Physical physical = checkPhysical(m);
    
  int dim = m->getDim();
  int nbNodes = m->getNumMeshVertices();
  
  //*************************************************************************
  
  if(dim == 1)
  {
    std::cout << "1D analysis..." << std::endl;
    
    if(param.method == Fem)
    {
      std::cout << "### ANALYTICAL ###" << std::endl;
      std::vector< std::complex<double> > uANALYTICAL = ANALYTICAL::solve(m, param, physical, false);
      writePOS(m, uANALYTICAL, "uANALYTICAL");
      std::cout << std::endl;
      
      std::cout << "### FEM ###" << std::endl;
      std::cout << "-> u" << std::endl;
      std::vector< std::complex<double> > uFEM = FEM::solve(m, param, physical);
      writePOS(m, uFEM, "uFEM");
      std::cout << "-> e" << std::endl;
      std::vector< std::complex<double> > eFEM = error(uANALYTICAL, uFEM);
      writePOS(m, eFEM, "eFEM");
      std::cout << std::endl;
      
      writeL2Error(uFEM, uANALYTICAL);
    }
    else if(param.method == Xfem)
    {
      std::cout << "### ANALYTICAL ###" << std::endl;
      std::vector< std::complex<double> > uANALYTICAL = ANALYTICAL::solve(m, param, physical, true);
      writePOS(m, uANALYTICAL, "uANALYTICAL");
      std::cout << std::endl;
      
      std::cout << "### XFEM ###" << std::endl;
      std::cout << "-> u" << std::endl;
      std::vector< std::complex<double> > uXFEM = XFEM::solve(m, param, physical);
      writePOS(m, uXFEM, "uXFEM");
      std::cout << "-> e" << std::endl;
      std::vector< std::complex<double> > eXFEM = error(uANALYTICAL, uXFEM);
      writePOS(m, eXFEM, "eXFEM");
      std::cout << std::endl;
      
      writeL2Error(uXFEM, uANALYTICAL);
    }
    else if(param.method == Lagrange)
    {
      std::cout << "### ANALYTICAL ###" << std::endl;
      std::vector< std::complex<double> > uANALYTICAL = ANALYTICAL::solve(m, param, physical, true);
      writePOS(m, uANALYTICAL, "uANALYTICAL");
      std::cout << std::endl;
      
      std::cout << "### LAGRANGE ###" << std::endl;
      std::cout << "-> u" << std::endl;
      std::vector< std::complex<double> > uLAGRANGE = LAGRANGE::solve(m, param, physical);
      writePOS(m, uLAGRANGE, "uLAGRANGE");
      std::cout << "-> e" << std::endl;
      std::vector< std::complex<double> > eLAGRANGE = error(uANALYTICAL, uLAGRANGE);
      writePOS(m, eLAGRANGE, "eLAGRANGE");
      std::cout << std::endl;
      
      writeL2Error(uLAGRANGE, uANALYTICAL);
    }
  }
  else if(dim == 2)
  {
    std::cout << "2D analysis..." << std::endl;
    if(param.method == Fem)
    {
      std::cout << "### ANALYTICAL ###" << std::endl;
      std::vector< std::complex<double> > uANALYTICAL = ANALYTICAL::solve(m, param, physical, false);
      writePOS(m, uANALYTICAL, "uANALYTICAL");
      std::cout << std::endl;
      
      std::cout << "### FEM ###" << std::endl;
      std::cout << "-> u" << std::endl;
      std::vector< std::complex<double> > uFEM = FEM::solve(m, param, physical);
      writePOS(m, uFEM, "uFEM");
      std::cout << "-> e" << std::endl;
      std::vector< std::complex<double> > eFEM = error(uANALYTICAL, uFEM);
      writePOS(m, eFEM, "eFEM");
      std::cout << std::endl;
      
      writeL2Error(uFEM, uANALYTICAL);
    }
    if(param.method == Xfem)
    {
      std::cout << "Not implemented yet !" << std::endl;
    }
    else if(param.method == Lagrange)
    {
      std::cout << "### ANALYTICAL ###" << std::endl;
      std::vector< std::complex<double> > uANALYTICAL = ANALYTICAL::solve(m, param, physical, true);
      writePOS(m, uANALYTICAL, "uANALYTICAL");
      std::cout << std::endl;
      
      std::cout << "### LAGRANGE ###" << std::endl;
      std::cout << "-> u" << std::endl;
      std::vector< std::complex<double> > uLAGRANGE = LAGRANGE::solve(m, param, physical);
      writePOS(m, uLAGRANGE, "uLAGRANGE");
      std::cout << "-> e" << std::endl;
      std::vector< std::complex<double> > eLAGRANGE = error(uANALYTICAL, uLAGRANGE);
      writePOS(m, eLAGRANGE, "eLAGRANGE");
      std::cout << std::endl;
      
      writeL2Error(uLAGRANGE, uANALYTICAL);
    }
  }
  
  
  delete m;
  GmshFinalize();
  return 1;
}

void writePOS(GModel* m, std::vector< std::complex<double> > u, std::string name)
{
  std::ofstream pos(name + ".pos", std::ofstream::trunc);
  
  pos << "View \"" << name << "\" {" << std::endl;
  pos << "TIME{0,0};" << std::endl;
  
  //Loop over vertices
  for(GModel::viter it = m->firstVertex(); it != m->lastVertex(); ++it)
  {
    GVertex *v = *it;
    
    for(unsigned int i = 0; i < v->points.size(); i++)
    {
      pos << "SP(" << v->points[i]->getVertex(0)->x() << ","  << v->points[i]->getVertex(0)->y() << ","  << v->points[i]->getVertex(0)->z() << "){" << u[v->points[i]->getVertex(0)->getNum()-1].real() << "," << u[v->points[i]->getVertex(0)->getNum()-1].imag() << "};" << std::endl;
    }
  }
  
  if(m->getDim() == 1)
  {
    //Loop over edges
    for(GModel::eiter it = m->firstEdge(); it != m->lastEdge(); ++it)
    {
      GEdge *e = *it;
    
      for(unsigned int i = 0; i < e->lines.size(); i++)
      {
        pos << "SL(" << e->lines[i]->getVertex(0)->x() << ","  << e->lines[i]->getVertex(0)->y() << ","  << e->lines[i]->getVertex(0)->z() << "," << e->lines[i]->getVertex(1)->x() << ","  << e->lines[i]->getVertex(1)->y() << ","  << e->lines[i]->getVertex(1)->z() << "){" << u[e->lines[i]->getVertex(0)->getNum()-1].real() << "," << u[e->lines[i]->getVertex(1)->getNum()-1].real() << "," << u[e->lines[i]->getVertex(0)->getNum()-1].imag() << "," << u[e->lines[i]->getVertex(1)->getNum()-1].imag() << "};" << std::endl;
      }
    }
  }
  else if(m->getDim() == 2)
  {
    //Loop over faces
    for(GModel::fiter it = m->firstFace(); it != m->lastFace(); ++it)
    {
      GFace *f = *it;
    
      for(unsigned int i = 0; i < f->triangles.size(); i++)
      {
        pos << "ST(";
        for(unsigned int j = 0; j < f->triangles[i]->getNumVertices(); j++)
        {
          pos << f->triangles[i]->getVertex(j)->x() << "," << f->triangles[i]->getVertex(j)->y() << "," << f->triangles[i]->getVertex(j)->z();
          if(j != f->triangles[i]->getNumVertices()-1)
          {
            pos << ",";
          }
        }
        pos << "){";
        for(unsigned int j = 0; j < f->triangles[i]->getNumVertices(); j++)
        {
          pos << u[f->triangles[i]->getVertex(j)->getNum()-1].real() << ",";
        }
        for(unsigned int j = 0; j < f->triangles[i]->getNumVertices(); j++)
        {
          pos << u[f->triangles[i]->getVertex(j)->getNum()-1].imag();
          if(j != f->triangles[i]->getNumVertices()-1)
          {
            pos << ",";
          }
        }
        pos << "};" << std::endl;
      }
    }
  }
  
  pos << "};" << std::endl;
  
  pos.close();
}

void writeL2Error(std::vector< std::complex<double> > u, std::vector< std::complex<double> > e)
{
  std::ofstream file("error.txt");
  
  double error = 0., num = 0., denum = 0.;
  
  for(unsigned int i = 0; i < u.size(); i++)
  {
    num += (std::norm(u[i]) - std::norm(e[i]))*(std::norm(u[i]) - std::norm(e[i]));
    denum += std::norm(e[i])*std::norm(e[i]);
  }
  
  error = num/denum;
  
  file << error;
  
  file.close();
}

Param readParam(int argc, char **argv)
{
  Param param;
  
  param.c_1 = 1;
  param.c_2 = 1;
  param.rho_1 = 1;
  param.rho_2 = 1;
  param.w = 50;
  param.x_bnd = 0.5;
  param.method = Fem;
  
  double phi = 0., A = 1.;
  
  for(unsigned int i = 1; i < argc; i++)
  {
    if(strcmp(argv[i], "-c1") == 0)
    {
      param.c_1 = atof(argv[i+1]);
      i++;
    }
    else if(strcmp(argv[i], "-c2") == 0)
    {
      param.c_2 = atof(argv[i+1]);
      i++;
    }
    else if(strcmp(argv[i], "-rho1") == 0)
    {
      param.rho_1 = atof(argv[i+1]);
      i++;
    }
    else if(strcmp(argv[i], "-rho2") == 0)
    {
      param.rho_2 = atof(argv[i+1]);
      i++;
    }
    else if(strcmp(argv[i], "-w") == 0)
    {
      param.w = atof(argv[i+1]);
      i++;
    }
    else if(strcmp(argv[i], "-x") == 0)
    {
      param.x_bnd = atof(argv[i+1]);
      i++;
    }
    else if(strcmp(argv[i], "-a") == 0)
    {
      A = atof(argv[i+1]);
      i++;
    }
    else if(strcmp(argv[i], "-phi") == 0)
    {
      phi = atof(argv[i+1]);
      i++;
    }
    else if(strcmp(argv[i], "-xfem") == 0)
    {
      param.method = Xfem;
      i++;
    }
    else if(strcmp(argv[i], "-fem") == 0)
    {
      param.method = Fem;
      i++;
    }
    else if(strcmp(argv[i], "-lag") == 0)
    {
      param.method = Lagrange;
      i++;
    }
    else if(strcmp(argv[i], "-prob") == 0)
    {
      if(atoi(argv[i+1]) == 0)
      {
        param.problem2D = Square;
      }
      else if(atoi(argv[i+1]) == 1)
      {
        param.problem2D = Cylinder;
      }
      else
      {
        param.problem2D = Square;
      }
    }
    else if(strcmp(argv[i], "-help") == 0)
    {
      std::cout << "-c1\t\tWave velocity imposed before x_bnd." << std::endl;
      std::cout << "-c2\t\tWave velocity imposed after x_bnd." << std::endl;
      std::cout << "-rho1\t\tDensity before x_bnd." << std::endl;
      std::cout << "-rho2\t\tDensity after x_bnd." << std::endl;
      std::cout << "-w\t\tThe angular frequency." << std::endl;
      std::cout << "-x\t\tLimit between two area with different wave number." << std::endl;
      std::cout << "-a\t\tAmplitude of the incidence wave." << std::endl;
      std::cout << "-phi\t\tPhase of the incidence wave." << std::endl;
      
      exit(1);
    }
  }
  param.wave = A*std::complex<double>(cos(phi),sin(phi));
  param.k_1 = param.w/param.c_1;
  param.k_2 = param.w/param.c_2;
  
  return param;
}

Physical checkPhysical(GModel *m)
{
  Physical phy;
  
  //Loop over faces
  for(GModel::fiter it = m->firstFace(); it != m->lastFace(); ++it)
  {
    GFace *f = *it;
    
    std::vector<int> physicals = f->physicals;
    
    for(unsigned int i = 0; i < physicals.size(); i++)
    {
      if(physicals[i] == OMEGA1)
      {
        phy.elmOmega1.push_back(f);
      }
      else if(physicals[i] == OMEGA2)
      {
        phy.elmOmega2.push_back(f);
      }
      else if(physicals[i] == OMEGA)
      {
        phy.elmOmega.push_back(f);
      }
    }
  }
  
  //Loop over edges
  for(GModel::eiter it = m->firstEdge(); it != m->lastEdge(); ++it)
  {
    GEdge *e = *it;
    
    std::vector<int> physicals = e->physicals;
    
    for(unsigned int i = 0; i < physicals.size(); i++)
    {
      if(physicals[i] == GAMMADIR)
      {
        phy.elmDir.push_back(e);
      }
      else if(physicals[i] == GAMMAINF)
      {
        phy.elmInf.push_back(e);
      }
      else if(physicals[i] == INTERFACE)
      {
        phy.elmInter.push_back(e);
      }
      else if(physicals[i] == OMEGA1)
      {
        phy.elmOmega1.push_back(e);
      }
      else if(physicals[i] == OMEGA2)
      {
        phy.elmOmega2.push_back(e);
      }
      else if(physicals[i] == OMEGA)
      {
        phy.elmOmega.push_back(e);
      }
    }
  }
  
  //Loop over vertices
  for(GModel::viter it = m->firstVertex(); it != m->lastVertex(); ++it)
  {
    GVertex *v = *it;
    
    std::vector<int> physicals = v->physicals;
    
    for(unsigned int i = 0; i < physicals.size(); i++)
    {
      if(physicals[i] == GAMMADIR)
      {
        phy.elmDir.push_back(v);
      }
      else if(physicals[i] == GAMMAINF)
      {
        phy.elmInf.push_back(v);
      }
      else if(physicals[i] == INTERFACE)
      {
        phy.elmInter.push_back(v);
      }
    }
    
  }
  
  return phy;
}



