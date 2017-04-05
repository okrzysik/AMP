#include <iostream>

#include "utils/AMPManager.h"
#include "utils/InputManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"

#include "ampmesh/Mesh.h"

#include "discretization/simpleDOF_Manager.h"
#include "matrices/MatrixBuilder.h"
#include "vectors/VectorBuilder.h"


void myTest( AMP::UnitTest *ut, std::string exeName )
{
#if !defined( USE_EXT_PETSC ) && !defined( USE_EXT_TRILINOS )
    if ( AMP::AMP_MPI( AMP_COMM_WORLD ).getSize() > 1 ) {
        ut->expected_failure( "No parallel matrix to test" );
        return;
    }
#endif
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;
    AMP::PIO::logOnlyNodeZero( log_file );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
    // build mesh database
    AMP::shared_ptr<AMP::Database> database( new AMP::MemoryDatabase( "Mesh" ) );
    database->putString( "MeshName", "cube" );
    database->putString( "MeshType", "AMP" );
    database->putInteger( "dim", 3 );
    database->putString( "Generator", "cube" );
    database->putDouble( "x_offset", 0.0 );
    database->putDouble( "y_offset", 0.0 );
    database->putDouble( "z_offset", 0.0 );
    std::vector<int> size( 3, 8 );
    database->putIntegerArray( "Size", size );
    std::vector<double> range( 6, 0.0 );
    range[1] = range[3] = range[5] = 1.0;
    database->putDoubleArray( "Range", range );
    // create mesh
    AMP::shared_ptr<AMP::Mesh::MeshParameters> params( new AMP::Mesh::MeshParameters( database ) );
    params->setComm( globalComm );
    AMP::shared_ptr<AMP::Mesh::Mesh> mesh = AMP::Mesh::Mesh::buildMesh( params );
    // create two different dof managers
    AMP::Discretization::DOFManager::shared_ptr firstDofManager =
        AMP::Discretization::simpleDOFManager::create( mesh, AMP::Mesh::GeomType::Vertex, 1, 1, true );
    AMP::Discretization::DOFManager::shared_ptr secondDofManager =
        AMP::Discretization::simpleDOFManager::create( mesh, AMP::Mesh::GeomType::Volume, 0, 2, true );
    size_t n = firstDofManager->numGlobalDOF();
    size_t m = secondDofManager->numGlobalDOF();
    AMP_ASSERT( n != m );
    // create the two corresponding vectors
    AMP::LinearAlgebra::Variable::shared_ptr var( new AMP::LinearAlgebra::Variable( "var" ) );
    AMP::LinearAlgebra::Vector::shared_ptr firstVec =
        AMP::LinearAlgebra::createVector( firstDofManager, var, true ); // n
    AMP::LinearAlgebra::Vector::shared_ptr secondVec =
        AMP::LinearAlgebra::createVector( secondDofManager, var, true ); // m
    firstVec->zero();
    secondVec->zero();
    // Loop through the valid matrix build types
    std::vector<int> types;
    types.push_back( 0 );
#if defined( USE_EXT_TRILINOS )
    types.push_back( 1 );
#endif
    if ( globalComm.getSize() == 1 )
        types.push_back( 2 );
    for ( auto &type : types ) {
        char tmp[100];
        sprintf( tmp, "%s: %i", exeName.c_str(), type );
        // create four matrices
        AMP::LinearAlgebra::Matrix::shared_ptr firstMat =
            AMP::LinearAlgebra::createMatrix( firstVec, secondVec, type ); // mxn
        AMP::LinearAlgebra::Matrix::shared_ptr secondMat =
            AMP::LinearAlgebra::createMatrix( secondVec, firstVec, type ); // nxm
        AMP_ASSERT( firstMat->numGlobalRows() == m && firstMat->numGlobalColumns() == n );
        firstMat->setScalar( 1.0 );
        secondMat->setScalar( 2.0 );
        // do matrix multiplications
        try {
            AMP::LinearAlgebra::Matrix::shared_ptr thirdMat =
                AMP::LinearAlgebra::Matrix::matMultiply( secondMat, firstMat );
            AMP::LinearAlgebra::Matrix::shared_ptr fourthMat =
                AMP::LinearAlgebra::Matrix::matMultiply( firstMat, secondMat );
            AMP_ASSERT( thirdMat->numGlobalRows() == n && thirdMat->numGlobalColumns() == n );
            AMP_ASSERT( fourthMat->numGlobalRows() == m && fourthMat->numGlobalColumns() == m );
            ut->passes( tmp );
        } catch ( ... ) {
            ut->failure( tmp );
        }
    }
}

int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    std::string exeName = "test_RectangularMatricesMultiplication";

    myTest( &ut, exeName );

    ut.report();

    int numFailed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return numFailed;
}
