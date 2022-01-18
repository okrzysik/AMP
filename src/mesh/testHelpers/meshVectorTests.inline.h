#ifndef included_AMP_MeshVectorTests
#define included_AMP_MeshVectorTests

#include "AMP/mesh/testHelpers/meshTests.h"

#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/DOF_ManagerParameters.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"

#include "AMP/vectors/testHelpers/VectorTests.h"


namespace AMP::Mesh {


// Simple nodal vector tests
template<int DOF_PER_NODE, bool SPLIT>
void meshTests::simpleNodalVectorTests( AMP::UnitTest &ut,
                                        AMP::Mesh::Mesh::shared_ptr mesh,
                                        std::shared_ptr<AMP::Discretization::DOFManager> DOFs,
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
        ut.passes( "global vector size" );
    else
        ut.failure( "global vector size: " + mesh->getName() );
    num_dofs = mesh->numLocalElements( AMP::Mesh::GeomType::Vertex ) * DOF_PER_NODE;
    if ( vectora->getLocalSize() == num_dofs )
        ut.passes( "local vector size" );
    else
        ut.failure( "local vector size: " + mesh->getName() );
    num_dofs = mesh->numGhostElements( AMP::Mesh::GeomType::Vertex, gcw ) * DOF_PER_NODE;
    if ( vectora->getGhostSize() == num_dofs )
        ut.passes( "ghost vector size" );
    else
        ut.failure( "ghost vector size: " + mesh->getName() );

    // Try some trival operations
    vectora->setRandomValues();
    double t1 = static_cast<double>( vectora->L2Norm() );
    vectora->abs( *vectora );
    if ( fabs( double( vectora->L2Norm() ) - t1 ) < 0.0000001 )
        ut.passes( "non-trivial random vector" );
    else
        ut.failure( "non-trivial random vector: + " + mesh->getName() );
    vectorb->setToScalar( 3. );
    vectora->multiply( *vectora, *vectorb );
    if ( fabs( double( vectora->L2Norm() ) - 3. * t1 ) < 0.00000001 )
        ut.passes( "trivial usage" );
    else
        ut.failure( "trivial usage: " + mesh->getName() );
}


// VerifyGetVectorTest
template<int DOF_PER_NODE, bool SPLIT>
void meshTests::VerifyGetVectorTest( AMP::UnitTest &ut, AMP::Mesh::Mesh::shared_ptr mesh )
{

    for ( int gcw = 0; gcw <= 1; gcw++ ) {

        // Create the DOF_Manager
        auto DOFparams = std::make_shared<AMP::Discretization::DOFManagerParameters>( mesh );
        auto DOFs      = AMP::Discretization::simpleDOFManager::create(
            mesh, AMP::Mesh::GeomType::Vertex, gcw, DOF_PER_NODE, SPLIT );

        // Run some basic nodal vector tests
        simpleNodalVectorTests<DOF_PER_NODE, SPLIT>( ut, mesh, DOFs, gcw );

        // Run the vector tests
        globalMeshForMeshVectorFactory = mesh;
        auto factory                   = std::make_shared<MeshVectorFactory>(
            mesh, AMP::Mesh::GeomType::Vertex, gcw, DOF_PER_NODE, SPLIT );
        AMP::LinearAlgebra::VectorTests vectorTests( factory );
        vectorTests.testManagedVector( &ut );
        vectorTests.testParallelVectors( &ut );
        globalMeshForMeshVectorFactory.reset();
    }
}


} // namespace AMP::Mesh

#endif
