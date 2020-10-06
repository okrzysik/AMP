#ifndef included_AMP_MeshVectorTests
#define included_AMP_MeshVectorTests

#include "AMP/ampmesh/testHelpers/meshTests.h"

#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/DOF_ManagerParameters.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"

#include "AMP/vectors/testHelpers/VectorTests.h"


namespace AMP {
namespace Mesh {


// Simple nodal vector tests
template<int DOF_PER_NODE, bool SPLIT>
void meshTests::simpleNodalVectorTests( AMP::UnitTest *utils,
                                        AMP::Mesh::Mesh::shared_ptr mesh,
                                        AMP::Discretization::DOFManager::shared_ptr DOFs,
                                        int gcw )
{

    // Create a nodal variable
    auto variable = std::make_shared<AMP::LinearAlgebra::Variable>( "test vector" );

    // Create the vectors
    auto vectora = AMP::LinearAlgebra::createVector( DOFs, variable, SPLIT );
    auto vectorb = AMP::LinearAlgebra::createVector( DOFs, variable, SPLIT );

    // Check the size of the vector
    size_t num_dofs = mesh->numGlobalElements( AMP::Mesh::GeomType::Vertex ) * DOF_PER_NODE;
    if ( vectora->getGlobalSize() == num_dofs )
        utils->passes( "global vector size" );
    else
        utils->failure( "global vector size: " + mesh->getName() );
    num_dofs = mesh->numLocalElements( AMP::Mesh::GeomType::Vertex ) * DOF_PER_NODE;
    if ( vectora->getLocalSize() == num_dofs )
        utils->passes( "local vector size" );
    else
        utils->failure( "local vector size: " + mesh->getName() );
    num_dofs = mesh->numGhostElements( AMP::Mesh::GeomType::Vertex, gcw ) * DOF_PER_NODE;
    if ( vectora->getGhostSize() == num_dofs )
        utils->passes( "ghost vector size" );
    else
        utils->failure( "ghost vector size: " + mesh->getName() );

    // Try some trival operations
    vectora->setRandomValues();
    double t1 = static_cast<double>(vectora->L2Norm());
    vectora->abs( *vectora );
    if ( fabs( double( vectora->L2Norm() ) - t1 ) < 0.0000001 )
        utils->passes( "non-trivial random vector" );
    else
        utils->failure( "non-trivial random vector: + " + mesh->getName() );
    vectorb->setToScalar( 3. );
    vectora->multiply( *vectora, *vectorb );
    if ( fabs( double( vectora->L2Norm() ) - 3. * t1 ) < 0.00000001 )
        utils->passes( "trivial usage" );
    else
        utils->failure( "trivial usage: " + mesh->getName() );
}


// VerifyGetVectorTest
template<int DOF_PER_NODE, bool SPLIT>
void meshTests::VerifyGetVectorTest( AMP::UnitTest *utils, AMP::Mesh::Mesh::shared_ptr mesh )
{

    for ( int gcw = 0; gcw <= 1; gcw++ ) {

        // Create the DOF_Manager
        auto DOFparams = std::make_shared<AMP::Discretization::DOFManagerParameters>( mesh );
        auto DOFs      = AMP::Discretization::simpleDOFManager::create(
            mesh, AMP::Mesh::GeomType::Vertex, gcw, DOF_PER_NODE, SPLIT );

        // Run some basic nodal vector tests
        simpleNodalVectorTests<DOF_PER_NODE, SPLIT>( utils, mesh, DOFs, gcw );

        // Run the vector tests
        globalMeshForMeshVectorFactory = mesh;
        auto factory                   = std::make_shared<MeshVectorFactory>(
            mesh, AMP::Mesh::GeomType::Vertex, gcw, DOF_PER_NODE, SPLIT );
        AMP::LinearAlgebra::VectorTests vectorTests( factory );
        vectorTests.testManagedVector( utils );
        vectorTests.testParallelVectors( utils );
        globalMeshForMeshVectorFactory.reset();
    }
}


} // namespace Mesh
} // namespace AMP

#endif
