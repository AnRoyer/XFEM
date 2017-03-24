#include <iostream>
#include <string>
#include <fstream>
#include <complex>
#include <cmath>

#include "Gmsh.h"
#include "GModel.h"
#include "MElement.h"
#include "MVertex.h"

#include "fem.h"

#define GAMMADIR 1
#define GAMMAINF 2
#define OMEGA 3

void writePOS(std::vector< std::complex<double> > u);
Param readParam(int argc, char **argv);
Physical checkPhysical(GModel *m);

int main(int argc, char **argv)
{
    std::cout << "###############################################" << std::endl;
    std::cout << "#                                             #" << std::endl;
    std::cout << "#   #       #  ########  ########  #      #   #" << std::endl;
    std::cout << "#    #     #   #         #         ##    ##   #" << std::endl;
    std::cout << "#     #   #    #         #         # #  # #   #" << std::endl;
    std::cout << "#      # #     #         #         #  ##  #   #" << std::endl;
    std::cout << "#       #      #####     #####     #      #   #" << std::endl;
    std::cout << "#      # #     #         #         #      #   #" << std::endl;
    std::cout << "#     #   #    #         #         #      #   #" << std::endl;
    std::cout << "#    #     #   #         #         #      #   #" << std::endl;
    std::cout << "#   #       #  #         ########  #      #   #" << std::endl;
    std::cout << "#                                             #" << std::endl;
    std::cout << "###############################################" << std::endl << std::endl;
    
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
    
    std::cout << physical.tagDir[0] << " " << physical.tagInf[0] << " " << param.wave << std::endl;
    
    int dim = m->getDim();
    int nbNodes = m->getNumMeshVertices();
    
    if(dim == 1)
    {
        std::cout << "1D analysis..." << std::endl;
        
        std::vector< std::complex<double> > u = solveFEM(m, nbNodes, param, physical);
        
        writePOS(u);
    }
    else if(dim == 2)
    {
        std::cout << "2D analysis..." << std::endl;
    }

    
    delete m;
    GmshFinalize();
    return 1;
}

void writePOS(std::vector< std::complex<double> > u)
{
    std::ofstream posR("u_real.pos", std::ofstream::trunc);
    
    posR << "$MeshFormat" << std::endl;
    posR << "2.2 0 8" << std::endl;
    posR << "$EndMeshFormat" << std::endl;
    posR << "$NodeData" << std::endl;
    posR << "1" << std::endl;
    posR << "\"u_real\"" << std::endl;
    posR << "1" << std::endl;
    posR << "1.0" << std::endl;
    posR << "3" << std::endl;
    posR << "0" << std::endl;
    posR << "1" << std::endl;
    posR << u.size() << std::endl;
    
    for(unsigned int i = 0; i < u.size(); i++)
    {
        posR << i+1 << " " << u[i].real() << std::endl;
    }
    posR << "$EndNodeData" << std::endl;
    
    posR.close();
    
    
    std::ofstream posI("u_imag.pos", std::ofstream::trunc);
    
    posI << "$MeshFormat" << std::endl;
    posI << "2.2 0 8" << std::endl;
    posI << "$EndMeshFormat" << std::endl;
    posI << "$NodeData" << std::endl;
    posI << "1" << std::endl;
    posI << "\"u_imag\"" << std::endl;
    posI << "1" << std::endl;
    posI << "1.0" << std::endl;
    posI << "3" << std::endl;
    posI << "0" << std::endl;
    posI << "1" << std::endl;
    posI << u.size() << std::endl;
    
    for(unsigned int i = 0; i < u.size(); i++)
    {
        posI << i+1 << " " << u[i].imag() << std::endl;
    }
    posI << "$EndNodeData" << std::endl;
    
    posI.close();
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
                std::vector<MVertex*> vertices = e->mesh_vertices;
                for(unsigned int j = 0; j < vertices.size(); j++)
                {
                    phy.tagInf.push_back(vertices[j]->getNum());
                }
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
                std::vector<MVertex*> vertices = v->mesh_vertices;
                for(unsigned int j = 0; j < vertices.size(); j++)
                {
                    phy.tagInf.push_back(vertices[j]->getNum());
                }
            }
        }
        
    }
    
    return phy;
}



