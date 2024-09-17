#ifndef included_AMP_Unit_test_Libmesh_Generators_hpp
#define included_AMP_Unit_test_Libmesh_Generators_hpp

#include "AMP/mesh/testHelpers/libmeshGenerators.h"


#ifdef AMP_USE_LIBMESH


    #include "AMP/mesh/libmesh/initializeLibMesh.h"
    #include "AMP/mesh/libmesh/libmeshMesh.h"

template<int SIZE>
void AMP::unit_test::LibMeshCubeGenerator<SIZE>::build_mesh()
{
    auto database = std::make_shared<AMP::Database>( "Mesh" );
    database->putScalar<int>( "dim", 3 );
    database->putScalar<std::string>( "MeshName", "cube_mesh" );
    database->putScalar<std::string>( "Generator", "cube" );
    database->putVector<int>( "size", std::vector<int>( 3, SIZE ) );
    database->putVector<double>( "xmin", std::vector<double>( 3, -1.0 ) );
    database->putVector<double>( "xmax", std::vector<double>( 3, 1.0 ) );
    auto params = std::make_shared<AMP::Mesh::MeshParameters>( database );
    params->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    mesh = std::make_shared<AMP::Mesh::libmeshMesh>( params );
}
template<int FILE>
void AMP::unit_test::ExodusReaderGenerator<FILE>::build_mesh()
{
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
    mesh = std::make_shared<AMP::Mesh::libmeshMesh>( params );
}


#else


    #include "AMP/utils/UtilityMacros.h"

template<int SIZE>
void AMP::unit_test::LibMeshCubeGenerator<SIZE>::build_mesh()
{
    AMP_ERROR( "LibMeshCubeGenerator requires libMesh" );
}
template<int FILE>
void AMP::unit_test::ExodusReaderGenerator<FILE>::build_mesh()
{
    AMP_ERROR( "ExodusReaderGenerator requires libMesh" );
}


#endif


#endif
