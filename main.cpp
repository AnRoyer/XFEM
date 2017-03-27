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
#include "error.h"

#define GAMMADIR 1
#define GAMMAINF 2
#define OMEGA 3

void writePOS(GModel* m, std::vector< std::complex<double> > u, std::string name);
Param readParam(int argc, char **argv);
Physical checkPhysical(GModel *m);

int main(int argc, char **argv)
{
    std::cout << "#################################################" << std::endl;
    std::cout << "#                                               #" << std::endl;
    std::cout << "#    #       #  ########  ########  #      #    #" << std::endl;
    std::cout << "#     #     #   #         #         ##    ##    #" << std::endl;
    std::cout << "#      #   #    #         #         # #  # #    #" << std::endl;
    std::cout << "#       # #     #         #         #  ##  #    #" << std::endl;
    std::cout << "#        #      #####     #####     #      #    #" << std::endl;
    std::cout << "#       # #     #         #         #      #    #" << std::endl;
    std::cout << "#      #   #    #         #         #      #    #" << std::endl;
    std::cout << "#     #     #   #         #         #      #    #" << std::endl;
    std::cout << "#    #       #  #         ########  #      #    #" << std::endl;
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
    
    if(dim == 1)
    {
        std::cout << "1D analysis..." << std::endl;
        
        std::cout << "### ANALYTICAL ###" << std::endl;
        
        std::vector< std::complex<double> > uANALYTICAL = ANALYTICAL::solve(m, nbNodes, param, physical);
        writePOS(m, uANALYTICAL, "uANALYTICAL");
        
        std::cout << std::endl;
        std::cout << "### FEM ###" << std::endl;
        std::cout << "-> u" << std::endl;
        std::vector< std::complex<double> > uFEM = FEM::solve(m, nbNodes, param, physical);
        writePOS(m, uFEM, "uFEM");
        std::cout << "-> e" << std::endl;
        std::vector< std::complex<double> > eFEM = error(uANALYTICAL, uFEM);
        writePOS(m, eFEM, "eFEM");
        
        std::cout << std::endl;
        std::cout << "### XFEM ###" << std::endl;
        std::cout << "-> u" << std::endl;
        std::vector< std::complex<double> > uXFEM = XFEM::solve(m, nbNodes, param, physical);
        writePOS(m, uXFEM, "uXFEM");
        std::cout << "-> e" << std::endl;
        std::vector< std::complex<double> > eXFEM = error(uANALYTICAL, uXFEM);
        writePOS(m, eXFEM, "eXFEM");
        
        std::cout << std::endl;
    }
    else if(dim == 2)
    {
        std::cout << "2D analysis..." << std::endl;
        
        std::cout << "### FEM ###" << std::endl;
        std::cout << "-> u" << std::endl;
        std::vector< std::complex<double> > uFEM = FEM::solve(m, nbNodes, param, physical);
        writePOS(m, uFEM, "uFEM");
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

Param readParam(int argc, char **argv)
{
    Param param;
    
    param.k_0 = 10;
    param.k_1 = 10;
    param.x_bnd = 0.5;
    
    double phi = 0., A = 1.;
    
    for(unsigned int i = 1; i < argc; i++)
    {
        if(strcmp(argv[i], "-k0") == 0)
        {
            param.k_0 = atof(argv[i+1]);
            i++;
        }
        else if(strcmp(argv[i], "-k1") == 0)
        {
            param.k_1 = atof(argv[i+1]);
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
        else if(strcmp(argv[i], "-help") == 0)
        {
            std::cout << "-k0\t\tWave number imposed before x_bnd." << std::endl;
            std::cout << "-k1\t\tWave number imposed after x_bnd." << std::endl;
            std::cout << "-x\t\tLimit between two area with different wave number." << std::endl;
            std::cout << "-a\t\tAmplitude of the incidence wave." << std::endl;
            std::cout << "-a\t\tPhase of the incidence wave." << std::endl;
            
            exit(1);
        }
    }
    param.wave = A*std::complex<double>(cos(phi),sin(phi));
    
    return param;
}

Physical checkPhysical(GModel *m)
{
    Physical phy;
    /*
    //Loop over faces
    for(GModel::fiter it = gModel->firstFace(); it != gModel->lastFace(); ++it)
    {
        GFace *f = *it;
        
        fillNodesToElements(nodesToElements, f->triangles.begin(), f->triangles.end());
        fillNodesToElements(nodesToElements, f->quadrangles.begin(), f->quadrangles.end());
        fillNodesToElements(nodesToElements, f->polygons.begin(), f->polygons.end());
    }
    */
    
    //Loop over edges
    for(GModel::eiter it = m->firstEdge(); it != m->lastEdge(); ++it)
    {
        GEdge *e = *it;
        
        std::vector<int> physicals = e->physicals;
        
        for(unsigned int i = 0; i < physicals.size(); i++)
        {
            if(physicals[i] == GAMMADIR)
            {
                std::vector<MVertex*> vertices = e->mesh_vertices;
                for(unsigned int j = 0; j < vertices.size(); j++)
                {
                    phy.tagDir.push_back(vertices[j]->getNum());
                }
            }
            else if(physicals[i] == GAMMAINF)
            {
                phy.elmInf.push_back(e);
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
                std::vector<MVertex*> vertices = v->mesh_vertices;
                for(unsigned int j = 0; j < vertices.size(); j++)
                {
                    phy.tagDir.push_back(vertices[j]->getNum());
                }
            }
            else if(physicals[i] == GAMMAINF)
            {
                phy.elmInf.push_back(v);
            }
        }
        
    }
    
    return phy;
}



