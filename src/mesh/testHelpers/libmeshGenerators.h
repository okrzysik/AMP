// This file contains classes for generating meshes that are based in libmesh
#ifndef included_AMP_Unit_test_Libmesh_Generators_h
#define included_AMP_Unit_test_Libmesh_Generators_h

#include "AMP/mesh/MeshParameters.h"
#include "AMP/mesh/libmesh/initializeLibMesh.h"
#include "AMP/mesh/libmesh/libmeshMesh.h"
#include "AMP/mesh/testHelpers/meshGenerators.h"


namespace libMesh::Parallel {
class Communicator;
}


namespace AMP::unit_test {


// Class to create a cube in Libmesh
template<int SIZE>
class LibMeshCubeGenerator : public MeshGenerator
{
public:
    void build_mesh() override
    {
        // Create the parameter object
        auto database = std::make_shared<AMP::Database>( "Mesh" );
        database->putScalar<int>( "dim", 3 );
        database->putScalar<std::string>( "MeshName", "cube_mesh" );
        database->putScalar<std::string>( "Generator", "cube" );
        database->putVector<int>( "size", std::vector<int>( 3, SIZE ) );
        database->putVector<double>( "xmin", std::vector<double>( 3, -1.0 ) );
        database->putVector<double>( "xmax", std::vector<double>( 3, 1.0 ) );
        auto params = std::make_shared<AMP::Mesh::MeshParameters>( database );
        params->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
        // Create a libMesh mesh
        mesh = std::make_shared<AMP::Mesh::libmeshMesh>( params );
    }

    static std::string name() { return "LibMeshCubeGenerator"; }
};


// Class to read in a default exodus file
template<int FILE = 1>
class ExodusReaderGenerator : public MeshGenerator
{
public:
    void build_mesh() override
    {
        // Create the parameter object
        auto database = std::make_shared<AMP::Database>( "Mesh" );
        database->putScalar( "dim", 3 );
        database->putScalar<std::string>( "MeshName", "exodus reader mesh" );
        if ( FILE == 1 ) {
            database->putScalar<std::string>( "FileName", "clad_1x_1pellet.e" );
        } else if ( FILE == 2 ) {
            database->putScalar<std::string>( "FileName", "multiElementMesh.e" );
        } else if ( FILE == 3 ) {
            database->putScalar<std::string>( "FileName", "pellet_1x.e" );
        } else {
            AMP_ERROR( "Bad file for generator" );
        }
        auto params = std::make_shared<AMP::Mesh::MeshParameters>( database );
        params->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
        // Create a libMesh mesh
        mesh = std::make_shared<AMP::Mesh::libmeshMesh>( params );
    }

    static std::string name() { return "ExodusReaderGenerator"; }
};


// MulitMesh generator
class MultiMeshGenerator : public MeshGenerator
{
public:
    void build_mesh() override;
    static std::string name() { return "MultiMeshGenerator"; }
};


// libMeshThreeElement generator
class libMeshThreeElementGenerator : public MeshGenerator
{
public:
    static std::string name() { return "libMeshThreeElementGenerator"; }

    static std::vector<unsigned int> getBndDofIndices();

    static std::vector<std::vector<unsigned int>> getElemNodeMap();

    void build_mesh() override;

    virtual ~libMeshThreeElementGenerator()
    {
        mesh.reset();
        libmeshInit.reset();
    }

protected:
    std::shared_ptr<libMesh::Parallel::Communicator> libMeshComm;
    std::shared_ptr<AMP::Mesh::initializeLibMesh> libmeshInit;
};


} // namespace AMP::unit_test


#endif
