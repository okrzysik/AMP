// This file contains classes for generating meshes that are used for different tests
#ifndef included_AMP_Unit_test_Mesh_Generators_h
#define included_AMP_Unit_test_Mesh_Generators_h

#include "ampmesh/Mesh.h"
#include "ampmesh/libmesh/libMesh.h"
#include "utils/MemoryDatabase.h"

namespace AMP {
namespace unit_test {


// Base class for Mesh Generators
class MeshGenerator
{
public:
    // Routine to build the mesh
    virtual void build_mesh() { AMP_ERROR("ERROR"); }
    // Routine to get the pointer to the mesh
    virtual AMP::Mesh::Mesh::shared_ptr getMesh() {
        if ( mesh.get()==NULL ) 
            this->build_mesh();
        return mesh;
    }
protected:
    AMP::Mesh::Mesh::shared_ptr  mesh;
};


// Class to create a cube in Libmesh
template <int SIZE>
class  LibMeshCubeGenerator : public MeshGenerator
{
public:

    virtual void build_mesh() {
        // Create the parameter object
        boost::shared_ptr<AMP::MemoryDatabase> database(new AMP::MemoryDatabase("Mesh"));
        database->putInteger("dim",3);
        database->putString("MeshName","cube_mesh");
        database->putString("Generator","cube");
        database->putIntegerArray("size",std::vector<int>(3,SIZE));
        database->putDoubleArray("xmin",std::vector<double>(3,-1.0));
        database->putDoubleArray("xmax",std::vector<double>(3,1.0));
        boost::shared_ptr<AMP::Mesh::MeshParameters> params(new AMP::Mesh::MeshParameters(database));
        params->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));
        // Create a libMesh mesh
        mesh = boost::shared_ptr<AMP::Mesh::libMesh> (new AMP::Mesh::libMesh(params));    
    }
};


// Class to read in a default exodus file
class  ExodusReaderGenerator : public MeshGenerator
{
public:
    virtual void build_mesh() {
        // Create the parameter object
        boost::shared_ptr<AMP::MemoryDatabase> database(new AMP::MemoryDatabase("Mesh"));
        database->putInteger("dim",3);
        database->putString("MeshName","exodus reader mesh");
        database->putString("FileName","clad_1x_1pellet.e");
        boost::shared_ptr<AMP::Mesh::MeshParameters> params(new AMP::Mesh::MeshParameters(database));
        params->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));
        // Create a libMesh mesh
        mesh = boost::shared_ptr<AMP::Mesh::libMesh>(new AMP::Mesh::libMesh(params));    
    }
};


// MulitMesh generator
class   MultiMeshGenerator : public MeshGenerator
{
public:
    virtual void build_mesh() {
        int N_meshes = 4;
        // Create the multimesh database
        boost::shared_ptr<AMP::MemoryDatabase> meshDatabase(new AMP::MemoryDatabase("Mesh"));
        meshDatabase->putString("MeshName","PelletMeshes");
        meshDatabase->putString("MeshType","Multimesh");
        meshDatabase->putString("MeshDatabasePrefix","Mesh_");
        meshDatabase->putString("MeshArrayDatabasePrefix","MeshArray_");
        // Create the mesh array database
        boost::shared_ptr<Database> meshArrayDatabase = meshDatabase->putDatabase("MeshArray_1");
        meshArrayDatabase->putInteger("N",N_meshes);
        meshArrayDatabase->putString("iterator","%i");
        std::vector<int> indexArray(N_meshes);
        for (int i=0; i<N_meshes; i++)
            indexArray[i] = i+1;
        meshArrayDatabase->putIntegerArray("indicies",indexArray);
        meshArrayDatabase->putString("MeshName","pellet_%i");
        meshArrayDatabase->putString("FileName","pellet_lo_res.e");
        meshArrayDatabase->putString("MeshType","libMesh");
        meshArrayDatabase->putInteger("dim",3);
        meshArrayDatabase->putDouble("x_offset",0.0);
        meshArrayDatabase->putDouble("y_offset",0.0);
        std::vector<double> offsetArray(N_meshes);
        for (int i=0; i<N_meshes; i++)
            offsetArray[i] = ((double) i)*0.0105;
        meshArrayDatabase->putDoubleArray("z_offset",offsetArray);
        meshArrayDatabase->putInteger("NumberOfElements",80);
        // Create the parameter object
        boost::shared_ptr<AMP::Mesh::MeshParameters> params(new AMP::Mesh::MeshParameters(meshDatabase));
        params->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));
        // Create the mesh
        mesh = AMP::Mesh::Mesh::buildMesh(params);
    }
};


// Surface subset generator
template <class GENERATOR> 
class   SurfaceSubsetGenerator : public MeshGenerator
{
public:
    virtual void build_mesh() {
        boost::shared_ptr<MeshGenerator> generator( new GENERATOR );
        generator->build_mesh();
        AMP::Mesh::Mesh::shared_ptr mesh1 = generator->getMesh();
        AMP::Mesh::GeomType type = mesh1->getGeomType();
        AMP::Mesh::GeomType type2 = (AMP::Mesh::GeomType) ((int) type - 1);
        AMP::Mesh::MeshIterator iterator = mesh1->getSurfaceIterator(type2,1);
        mesh = mesh1->Subset(iterator);
    }
};


// ThreeElement generator
/*class   ThreeElementLGenerator : public MeshGenerator
{
public:

    static std::vector<unsigned int> getBndDofIndices() {
        std::vector<unsigned int> bndDofIndices(4);
        bndDofIndices[0] = 0;
        bndDofIndices[1] = 3;
        bndDofIndices[2] = 4;
        bndDofIndices[3] = 7;
        return bndDofIndices;
    }

    static std::vector<std::vector<unsigned int> > getElemNodeMap() {
        std::vector<std::vector<unsigned int> > elemNodeMap(3);

        elemNodeMap[0].resize(8);
        elemNodeMap[1].resize(8);
        elemNodeMap[2].resize(8);

        elemNodeMap[0][0] = 0;
        elemNodeMap[0][1] = 1;
        elemNodeMap[0][2] = 2;
        elemNodeMap[0][3] = 3;
        elemNodeMap[0][4] = 4;
        elemNodeMap[0][5] = 5;
        elemNodeMap[0][6] = 6;
        elemNodeMap[0][7] = 7;

        elemNodeMap[1][0] = 1;
        elemNodeMap[1][1] = 8;
        elemNodeMap[1][2] = 9;
        elemNodeMap[1][3] = 2;
        elemNodeMap[1][4] = 5;
        elemNodeMap[1][5] = 10;
        elemNodeMap[1][6] = 11;
        elemNodeMap[1][7] = 6;

        elemNodeMap[2][0] = 2;
        elemNodeMap[2][1] = 9;
        elemNodeMap[2][2] = 12;
        elemNodeMap[2][3] = 13;
        elemNodeMap[2][4] = 6;
        elemNodeMap[2][5] = 11;
        elemNodeMap[2][6] = 14;
        elemNodeMap[2][7] = 15;

        return elemNodeMap;
    }

    virtual void build_mesh() {
        const unsigned int mesh_dim = 3;
        const unsigned int num_elem = 3;
        const unsigned int num_nodes = 16;

        boost::shared_ptr< ::Mesh > local_mesh(new ::Mesh(mesh_dim));
        local_mesh->reserve_elem(num_elem);
        local_mesh->reserve_nodes(num_nodes);

        local_mesh->add_point(::Point(0.0, 0.0, 0.0), 0);
        local_mesh->add_point(::Point(0.5, 0.0, 0.0), 1);
        local_mesh->add_point(::Point(0.5, 0.5, 0.0), 2);
        local_mesh->add_point(::Point(0.0, 0.5, 0.0), 3);
        local_mesh->add_point(::Point(0.0, 0.0, 0.5), 4);
        local_mesh->add_point(::Point(0.5, 0.0, 0.5), 5);
        local_mesh->add_point(::Point(0.5, 0.5, 0.5), 6);
        local_mesh->add_point(::Point(0.0, 0.5, 0.5), 7);
        local_mesh->add_point(::Point(1.0, 0.0, 0.0), 8);
        local_mesh->add_point(::Point(1.0, 0.5, 0.0), 9);
        local_mesh->add_point(::Point(1.0, 0.0, 0.5), 10);
        local_mesh->add_point(::Point(1.0, 0.5, 0.5), 11);
        local_mesh->add_point(::Point(1.0, 1.0, 0.0), 12);
        local_mesh->add_point(::Point(0.5, 1.0, 0.0), 13);
        local_mesh->add_point(::Point(1.0, 1.0, 0.5), 14);
        local_mesh->add_point(::Point(0.5, 1.0, 0.5), 15);

        std::vector<std::vector<unsigned int> > elemNodeMap = getElemNodeMap();

        for(size_t i = 0; i < elemNodeMap.size(); i++) {
            ::Elem* elem = local_mesh->add_elem(new ::Hex8);
            for(int j = 0; j < 8; j++) {
                elem->set_node(j) = local_mesh->node_ptr(elemNodeMap[i][j]);
            }
        }

        const short int boundaryId = 1;
        std::vector<unsigned int> bndDofIndices = getBndDofIndices(); 

        for(size_t i = 0; i < bndDofIndices.size(); i++) {
            local_mesh->boundary_info->add_node(local_mesh->node_ptr(bndDofIndices[i]), boundaryId);
        }

        local_mesh->prepare_for_use(true);

        mesh = AMP::Mesh::MeshAdapter::shared_ptr ( new AMP::Mesh::MeshAdapter (local_mesh) );
    }

};*/

 
}
}

#endif
