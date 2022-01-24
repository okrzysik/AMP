//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   operators/test/testNekOperator.cc
 * \brief  Tests Nek pipe problem run through NekMoabOperator
 *
 * This test is intended to test our ability to use a Nek-generated Moab
 * instance.  It duplicates the behavior of testNekPipe, but maps the Nek
 * output onto an actual AMP mesh.
 */
//---------------------------------------------------------------------------//

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "AMP/IO/PIO.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"

#include "AMP/IO/Writer.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/vectors/VectorBuilder.h"

#include "AMP/operators/moab/MoabMapOperator.h"

// Nek includes
#include "nek/NekMoabOperator.h"
#include "nek/NekMoabOperatorParameters.h"

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

static void nekPipeOperator( AMP::UnitTest *ut )
{
#ifdef USE_EXT_NEK
    // Print Banner
    AMP::Utilities::printBanner();

    // Log all nodes
    AMP::logAllNodes( "output_testNekOperator" );

    // Build database
    AMP::pout << "Building Input Database" << std::endl;
    auto nekDB = std::make_shared<AMP::Database>( "Nek_DB" );
    nekDB->putScalar<std::string>( "NekProblemName", "pipe" );

    // Build operator params
    AMP::pout << "Building Nek Operator Parameters" << std::endl;
    auto nekParams = std::make_shared<AMP::Operator::NekMoabOperatorParameters>( nekDB );

    // Build operator
    AMP::pout << "Building Nek Operator" << std::endl;
    auto nekOp = std::make_shared<AMP::Operator::NekMoabOperator>( nekParams );

    // Call apply
    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    nekOp->apply( nullVec, nullVec, nullVec, 0.0, 0.0 );

    // Read AMP pellet mesh from file
    auto meshDB = nekDB->putDatabase( "Mesh" );
    meshDB->putScalar<std::string>( "FileName", "pellet_1x.e" );
    meshDB->putScalar<std::string>( "MeshName", "fuel" );
    meshDB->putScalar<std::string>( "MeshType", "libMesh" );
    meshDB->putScalar<int>( "dim", 3 );
    meshDB->putScalar<double>( "x_offset", 0.0 );
    meshDB->putScalar<double>( "y_offset", 0.0 );
    meshDB->putScalar<double>( "z_offset", 0.0 );
    meshDB->putScalar<int>( "NumberOfElements", 300 );


    // Create Mesh
    AMP::pout << "Creating AMP mesh" << std::endl;
    auto meshParams = std::make_shared<AMP::Mesh::MeshParameters>( meshDB );
    meshParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    auto mesh = AMP::Mesh::Mesh::buildMesh( meshParams );


    // Create Parameters for Map Operator
    AMP::pout << "Creating map operator" << std::endl;
    nekDB->putScalar<std::string>( "MoabMapVariable", "VPRESS" );
    nekDB->putScalar<std::string>( "InterpolateToType", "GaussPoint" );
    auto mapParams = std::make_shared<AMP::Operator::MoabMapOperatorParameters>( nekDB );
    mapParams->setMoabOperator( nekOp );
    mapParams->setMesh( mesh );

    AMP::pout << "Creating GP-Based Moab Map Operator" << std::endl;
    auto moabGPMap = std::make_shared<AMP::Operator::MoabMapOperator>( mapParams );

    // Create variable to hold pressure data
    auto allGPPressures =
        std::make_shared<AMP::LinearAlgebra::Variable>( "AllGaussPointPressures" );
    auto allNodePressures = std::make_shared<AMP::LinearAlgebra::Variable>( "AllNodalPressures" );

    // Create DOF managers
    size_t DOFsPerElement    = 8;
    size_t DOFsPerNode       = 1;
    int gaussPointGhostWidth = 1;
    int nodalGhostWidth      = 1;
    bool split               = true;
    auto gaussPointDofMap    = AMP::Discretization::simpleDOFManager::create(
        mesh, AMP::Mesh::GeomType::Volume, gaussPointGhostWidth, DOFsPerElement, split );
    auto nodalDofMap = AMP::Discretization::simpleDOFManager::create(
        mesh, AMP::Mesh::GeomType::Vertex, nodalGhostWidth, DOFsPerNode, split );

    // Have mesh manager create vector over all meshes
    auto r_gp   = AMP::LinearAlgebra::createVector( gaussPointDofMap, allGPPressures );
    auto r_node = AMP::LinearAlgebra::createVector( nodalDofMap, allNodePressures );

    AMP::pout << "Gauss Point MultiVector size: " << r_gp->getGlobalSize() << std::endl;
    AMP::pout << "Nodal MultiVector size: " << r_node->getGlobalSize() << std::endl;

    AMP::pout << "Calling apply" << std::endl;
    moabGPMap->apply( nullVec, nullVec, r_gp, 0.0, 0.0 );

    AMP::pout << "Creating Node-Based Moab Map Operator" << std::endl;
    nekDB->putScalar<std::string>( "InterpolateToType", "GeomType::Vertex" );
    auto moabNodeMap = std::make_shared<MoabMap>( mapParams );

    moabNodeMap->apply( nullVec, nullVec, r_node, 0.0, 0.0 );

    // Did we actually get data?
    typedef AMP::LinearAlgebra::Vector AMPVec;
    AMPVec::iterator myIter;
    int ctr      = 0;
    bool nonZero = false;
    for ( myIter = r_gp->begin(); myIter != r_gp->end(); ++myIter ) {
        // AMP::pout << "GP Vector Element " << ctr << " is " << *myIter << std::endl;

        if ( *myIter != 0.0 )
            nonZero = true;

        ctr++;
    }

    if ( nonZero )
        ut->passes( "Gauss point vector is not identically zero" );
    else
        ut->failure( "Gauss point vector is identically zero" );

    ctr     = 0;
    nonZero = false;
    for ( myIter = r_node->begin(); myIter != r_node->end(); ++myIter ) {
        if ( *myIter != 0.0 )
            nonZero = true;

        ctr++;
    }

    if ( nonZero )
        ut->passes( "Nodal vector is not identically zero" );
    else
        ut->failure( "Nodal vector is identically zero" );

    // How about some output?
    auto siloWriter = AMP::IO::Writer::buildWriter( "Silo" );
    siloWriter->registerMesh( mesh );
    siloWriter->registerVector( r_gp, mesh, AMP::Mesh::GeomType::Volume, "AllGaussPointPressures" );
    siloWriter->registerVector( r_node, mesh, AMP::Mesh::GeomType::Vertex, "AllNodalPressures" );
    siloWriter->writeFile( "Nek_Pressure", 0 );

    // Finalize Nek Operator
    nekOp->finalize();

#else
    ut->passes( "Nek was not used." );
#endif

    if ( ut->NumPassGlobal() == 0 )
        ut->failure( "if it doesn't pass, it must have failed." );
}


int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    nekPipeOperator( &ut );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}

//---------------------------------------------------------------------------//
//                 end of testNekOperator.cc
//---------------------------------------------------------------------------//
