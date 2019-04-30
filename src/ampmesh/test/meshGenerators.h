// This file contains classes for generating meshes that are used for different tests
#ifndef included_AMP_Unit_test_Mesh_Generators_h
#define included_AMP_Unit_test_Mesh_Generators_h

#include "AMP/ampmesh/Mesh.h"
#include "AMP/ampmesh/structured/BoxMesh.h"
#include "AMP/utils/MemoryDatabase.h"

namespace AMP {
namespace unit_test {


// Base class for Mesh Generators
class MeshGenerator
{
public:
    // Routine to build the mesh
    virtual void build_mesh() { AMP_ERROR( "ERROR" ); }
    // Routine to get the pointer to the mesh
    virtual AMP::Mesh::Mesh::shared_ptr getMesh()
    {
        if ( mesh.get() == nullptr )
            this->build_mesh();
        return mesh;
    }
    virtual ~MeshGenerator(){};

protected:
    AMP::Mesh::Mesh::shared_ptr mesh;
};


// Class to create a cube
template<int SIZE_X, int SIZE_Y, int SIZE_Z>
class AMPCubeGenerator3 : public MeshGenerator
{
public:
    virtual void build_mesh() override
    {
        // Create a generic MeshParameters object
        auto database = AMP::make_shared<AMP::MemoryDatabase>( "Mesh" );
        database->putInteger( "dim", 3 );
        database->putString( "MeshName", "AMP::cube" );
        database->putString( "Generator", "cube" );
        database->putIntegerArray( "Size", { SIZE_X, SIZE_Y, SIZE_Z } );
        database->putDoubleArray( "Range", { 0.0, 1.0, 0.0, 1.0, 0.0, 1.0 } );
        auto params = AMP::make_shared<AMP::Mesh::MeshParameters>( database );
        params->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
        // Create an AMP mesh
        mesh = AMP::Mesh::BoxMesh::generate( params );
    }
};
template<int SIZE>
class AMPCubeGenerator : public MeshGenerator
{
public:
    virtual void build_mesh() override
    {
        AMPCubeGenerator3<SIZE, SIZE, SIZE> gen;
        gen.build_mesh();
        mesh = gen.getMesh();
    }
    static std::string name()
    {
        char tmp[128];
        sprintf( tmp, "AMPCubeGenerator<%i>", SIZE );
        return std::string( tmp );
    }
};


// Class to create a cylinder
class AMPCylinderGenerator : public MeshGenerator
{
public:
    virtual void build_mesh() override
    {
        // Create a generic MeshParameters object
        auto database = AMP::make_shared<AMP::MemoryDatabase>( "Mesh" );
        database->putInteger( "dim", 3 );
        database->putString( "MeshName", "AMP::cylinder" );
        database->putString( "Generator", "cylinder" );
        database->putIntegerArray( "Size", { 10, 10 } );
        database->putDoubleArray( "Range", { 1.0, 0.0, 1.0 } );
        auto params = AMP::make_shared<AMP::Mesh::MeshParameters>( database );
        params->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
        // Create an AMP mesh
        mesh = AMP::Mesh::BoxMesh::generate( params );
    }
    static std::string name() { return "AMPCylinderGenerator"; }
};


// Class to create a tube
class AMPTubeGenerator : public MeshGenerator
{
public:
    virtual void build_mesh() override
    {
        // Create a generic MeshParameters object
        auto database = AMP::make_shared<AMP::MemoryDatabase>( "Mesh" );
        database->putInteger( "dim", 3 );
        database->putString( "MeshName", "AMP::tube" );
        database->putString( "Generator", "tube" );
        database->putIntegerArray( "Size", { 3, 12, 10 } );
        database->putDoubleArray( "Range", { 0.7, 1.0, 0.0, 1.0 } );
        auto params = AMP::make_shared<AMP::Mesh::MeshParameters>( database );
        params->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
        // Create an AMP mesh
        mesh = AMP::Mesh::BoxMesh::generate( params );
    }
    static std::string name() { return "AMPTubeGenerator"; }
};


// MulitMesh generator
class AMPMultiMeshGenerator : public MeshGenerator
{
public:
    virtual void build_mesh() override
    {
        // Create the multimesh database
        auto meshDatabase = AMP::make_shared<AMP::MemoryDatabase>( "Mesh" );
        meshDatabase->putString( "MeshName", "SinglePin" );
        meshDatabase->putString( "MeshType", "Multimesh" );
        meshDatabase->putString( "MeshDatabasePrefix", "Mesh_" );
        meshDatabase->putString( "MeshArrayDatabasePrefix", "MeshArray_" );
        // Create the mesh array database (PelletMeshes)
        auto pelletMeshDatabase = meshDatabase->putDatabase( "Mesh_1" );
        createPelletMeshDatabase( pelletMeshDatabase );
        // Create the mesh database (clad)
        auto cladMeshDatabase = meshDatabase->putDatabase( "Mesh_2" );
        createCladMeshDatabase( cladMeshDatabase );
        // Create the parameter object
        auto params = AMP::make_shared<AMP::Mesh::MeshParameters>( meshDatabase );
        params->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
        // Create the mesh
        mesh = AMP::Mesh::Mesh::buildMesh( params );
    }
    static std::string name() { return "AMPMultiMeshGenerator"; }

private:
    void createPelletMeshDatabase( AMP::shared_ptr<Database> db )
    {
        int N_pellet = 3;
        // Create the multimesh database
        db->putString( "MeshName", "PelletMeshes" );
        db->putString( "MeshType", "Multimesh" );
        db->putString( "MeshDatabasePrefix", "Mesh_" );
        db->putString( "MeshArrayDatabasePrefix", "MeshArray_" );
        // Create the mesh array database (PelletMeshes)
        auto meshArrayDatabase = db->putDatabase( "MeshArray_1" );
        meshArrayDatabase->putInteger( "N", N_pellet );
        meshArrayDatabase->putString( "iterator", "%i" );
        std::vector<int> indexArray( N_pellet );
        for ( int i = 0; i < N_pellet; i++ )
            indexArray[i] = i + 1;
        meshArrayDatabase->putIntegerArray( "indicies", indexArray );
        meshArrayDatabase->putString( "MeshName", "pellet_%i" );
        meshArrayDatabase->putString( "MeshType", "AMP" );
        meshArrayDatabase->putString( "Generator", "cylinder" );
        meshArrayDatabase->putIntegerArray( "Size", { 5, 8 } );
        meshArrayDatabase->putDoubleArray( "Range", { 0.004025, 0.0, 0.0105 } );
        meshArrayDatabase->putInteger( "dim", 3 );
        meshArrayDatabase->putDouble( "x_offset", 0.0 );
        meshArrayDatabase->putDouble( "y_offset", 0.0 );
        std::vector<double> offsetArray( N_pellet );
        for ( int i = 0; i < N_pellet; i++ )
            offsetArray[i] = ( (double) i ) * 0.0105;
        meshArrayDatabase->putDoubleArray( "z_offset", offsetArray );
    }
    void createCladMeshDatabase( AMP::shared_ptr<Database> db )
    {
        // Create the multimesh database
        db->putString( "MeshName", "clad" );
        db->putString( "MeshType", "AMP" );
        db->putString( "Generator", "tube" );
        db->putIntegerArray( "Size", { 3, 36, 32 } );
        db->putDoubleArray( "Range", { 0.00411, 0.00475, 0, 0.0315 } );
        db->putInteger( "dim", 3 );
    }
};


// Surface subset generator
template<class GENERATOR, int GCW>
class SurfaceSubsetGenerator : public MeshGenerator
{
public:
    virtual void build_mesh() override
    {
        GENERATOR generator;
        generator.build_mesh();
        auto mesh1 = generator.getMesh();
        auto type  = mesh1->getGeomType();
        auto type2 = ( AMP::Mesh::GeomType )( (int) type - 1 );
        auto it    = mesh1->getSurfaceIterator( type2, GCW );
        mesh       = mesh1->Subset( it );
    }
    static std::string name() { return "SurfaceSubsetGenerator"; }
};


} // namespace unit_test
} // namespace AMP


// Include libmesh generators
#ifdef USE_EXT_LIBMESH
#include "libmeshGenerators.h"
#endif


#endif
