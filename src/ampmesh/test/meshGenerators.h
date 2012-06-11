// This file contains classes for generating meshes that are used for different tests
#ifndef included_AMP_Unit_test_Mesh_Generators_h
#define included_AMP_Unit_test_Mesh_Generators_h

#include "ampmesh/Mesh.h"
#include "ampmesh/structured/BoxMesh.h"
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


// Class to create a cube 
template <int SIZE>
class  AMPCubeGenerator : public MeshGenerator
{
public:
    virtual void build_mesh() {
        // Set the dimensions of the mesh
        std::vector<int> size(3,SIZE);
        std::vector<double> range(6,0.0);
        range[1] = 1.0;
        range[3] = 1.0;
        range[5] = 1.0;
        // Create a generic MeshParameters object
        boost::shared_ptr<AMP::MemoryDatabase> database(new AMP::MemoryDatabase("Mesh"));
        database->putInteger("dim",3);
        database->putString("MeshName","mesh1");
        database->putString("Generator","cube");
        database->putIntegerArray("Size",size);
        database->putDoubleArray("Range",range);
        boost::shared_ptr<AMP::Mesh::MeshParameters> params(new AMP::Mesh::MeshParameters(database));
        params->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));
        // Create an AMP mesh
        mesh = boost::shared_ptr<AMP::Mesh::BoxMesh>(new AMP::Mesh::BoxMesh(params));      
    }
};


// Class to create a cylinder 
class  AMPCylinderGenerator : public MeshGenerator
{
public:
    virtual void build_mesh() {
        // Set the dimensions of the mesh
        std::vector<int> size(2);
        size[0] = 10;
        size[1] = 10;
        std::vector<double> range(3);
        range[0] = 1.0;
        range[1] = 0.0;
        range[2] = 1.0;
        // Create a generic MeshParameters object
        boost::shared_ptr<AMP::MemoryDatabase> database(new AMP::MemoryDatabase("Mesh"));
        database->putInteger("dim",3);
        database->putString("MeshName","mesh1");
        database->putString("Generator","cylinder");
        database->putIntegerArray("Size",size);
        database->putDoubleArray("Range",range);
        boost::shared_ptr<AMP::Mesh::MeshParameters> params(new AMP::Mesh::MeshParameters(database));
        params->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));
        // Create an AMP mesh
        mesh = boost::shared_ptr<AMP::Mesh::BoxMesh>(new AMP::Mesh::BoxMesh(params));      
    }
};


// Class to create a tube 
class  AMPTubeGenerator : public MeshGenerator
{
public:
    virtual void build_mesh() {
        // Set the dimensions of the mesh
        std::vector<int> size(3);
        size[0] = 3;
        size[1] = 12;
        size[2] = 10;
        std::vector<double> range(4);
        range[0] = 0.7;
        range[1] = 1.0;
        range[2] = 0.0;
        range[3] = 1.0;
        // Create a generic MeshParameters object
        boost::shared_ptr<AMP::MemoryDatabase> database(new AMP::MemoryDatabase("Mesh"));
        database->putInteger("dim",3);
        database->putString("MeshName","mesh1");
        database->putString("Generator","tube");
        database->putIntegerArray("Size",size);
        database->putDoubleArray("Range",range);
        boost::shared_ptr<AMP::Mesh::MeshParameters> params(new AMP::Mesh::MeshParameters(database));
        params->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));
        // Create an AMP mesh
        mesh = boost::shared_ptr<AMP::Mesh::BoxMesh>(new AMP::Mesh::BoxMesh(params));      
    }
};


// MulitMesh generator
class AMPMultiMeshGenerator : public MeshGenerator
{
public:
    virtual void build_mesh() {
        // Create the multimesh database
        boost::shared_ptr<AMP::MemoryDatabase> meshDatabase(new AMP::MemoryDatabase("Mesh"));
        meshDatabase->putString("MeshName","SinglePin");
        meshDatabase->putString("MeshType","Multimesh");
        meshDatabase->putString("MeshDatabasePrefix","Mesh_");
        meshDatabase->putString("MeshArrayDatabasePrefix","MeshArray_");
        // Create the mesh array database (PelletMeshes)
        boost::shared_ptr<Database> pelletMeshDatabase = meshDatabase->putDatabase("Mesh_1");
        createPelletMeshDatabase( pelletMeshDatabase );
        // Create the mesh database (clad)
        boost::shared_ptr<Database> cladMeshDatabase = meshDatabase->putDatabase("Mesh_2");
        createCladMeshDatabase( cladMeshDatabase );
        // Create the parameter object
        boost::shared_ptr<AMP::Mesh::MeshParameters> params(new AMP::Mesh::MeshParameters(meshDatabase));
        params->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));
        // Create the mesh
        mesh = AMP::Mesh::Mesh::buildMesh(params);
    }
private:
    void createPelletMeshDatabase( boost::shared_ptr<Database> db ) {
        int N_pellet = 2;
        // Create the multimesh database
        boost::shared_ptr<AMP::MemoryDatabase> meshDatabase(new AMP::MemoryDatabase("Mesh"));
        db->putString("MeshName","PelletMeshes");
        db->putString("MeshType","Multimesh");
        db->putString("MeshDatabasePrefix","Mesh_");
        db->putString("MeshArrayDatabasePrefix","MeshArray_");
        // Create the mesh array database (PelletMeshes)
        boost::shared_ptr<Database> meshArrayDatabase = db->putDatabase("MeshArray_1");
        meshArrayDatabase->putInteger("N",N_pellet);
        meshArrayDatabase->putString("iterator","%i");
        std::vector<int> indexArray(N_pellet);
        for (int i=0; i<N_pellet; i++)
            indexArray[i] = i+1;
        meshArrayDatabase->putIntegerArray("indicies",indexArray);
        meshArrayDatabase->putString("MeshName","pellet_%i");
        std::vector<int> size(2);
        size[0] = 5;
        size[1] = 8;
        std::vector<double> range(3);
        range[0] = 0.004025;
        range[1] = 0;
        range[2] = 0.0105;
        meshArrayDatabase->putString("MeshType","AMP");
        meshArrayDatabase->putString("Generator","cylinder");
        meshArrayDatabase->putIntegerArray("Size",size);
        meshArrayDatabase->putDoubleArray("Range",range);
        meshArrayDatabase->putInteger("dim",3);
        meshArrayDatabase->putDouble("x_offset",0.0);
        meshArrayDatabase->putDouble("y_offset",0.0);
        std::vector<double> offsetArray(N_pellet);
        for (int i=0; i<N_pellet; i++)
            offsetArray[i] = ((double) i)*0.0105;
        meshArrayDatabase->putDoubleArray("z_offset",offsetArray);
    }

    void createCladMeshDatabase( boost::shared_ptr<Database> db ) {
        int N_meshes = 3;
        std::vector<int> size(3);
        std::vector<double> range(4);
        size[0] = 3;
        size[1] = 36;
        size[2] = 32;
        range[0] = 0.00411;
        range[1] = 0.00475;
        range[2] = 0;
        range[3] = 0.042;
        // Create the multimesh database
        boost::shared_ptr<AMP::MemoryDatabase> meshDatabase(new AMP::MemoryDatabase("Mesh"));
        db->putString("MeshName","clad");
        db->putString("MeshType","AMP");
        db->putString("Generator","tube");
        db->putIntegerArray("Size",size);
        db->putDoubleArray("Range",range);
        db->putInteger("dim",3);
    }

};


// Surface subset generator
template <class GENERATOR,int GCW> 
class   SurfaceSubsetGenerator : public MeshGenerator
{
public:
    virtual void build_mesh() {
        boost::shared_ptr<MeshGenerator> generator( new GENERATOR );
        generator->build_mesh();
        AMP::Mesh::Mesh::shared_ptr mesh1 = generator->getMesh();
        AMP::Mesh::GeomType type = mesh1->getGeomType();
        AMP::Mesh::GeomType type2 = (AMP::Mesh::GeomType) ((int) type - 1);
        AMP::Mesh::MeshIterator iterator = mesh1->getSurfaceIterator(type2,GCW);
        mesh = mesh1->Subset(iterator);
    }
};


 
}
}


// Include libmesh generators
#ifdef USE_LIBMESH
    #include "libmeshGenerators.h"
#endif


#endif
