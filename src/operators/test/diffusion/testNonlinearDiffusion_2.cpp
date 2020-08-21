#include "AMP/ampmesh/Mesh.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/materials/Material.h"
#include "AMP/operators/ElementPhysicsModelFactory.h"
#include "AMP/operators/ElementPhysicsModelParameters.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/diffusion/DiffusionConstants.h"
#include "AMP/operators/diffusion/DiffusionLinearElement.h"
#include "AMP/operators/diffusion/DiffusionLinearFEOperator.h"
#include "AMP/operators/diffusion/DiffusionLinearFEOperatorParameters.h"
#include "AMP/operators/diffusion/DiffusionNonlinearElement.h"
#include "AMP/operators/diffusion/DiffusionNonlinearFEOperator.h"
#include "AMP/operators/diffusion/DiffusionNonlinearFEOperatorParameters.h"
#include "AMP/operators/diffusion/DiffusionTransportModel.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/VectorBuilder.h"

#include "patchfunctions.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>


/**
 * This test is patch test for the diffusion operator.
 */
static void nonlinearTest( AMP::UnitTest *ut,
                           const std::string &exeName,
                           double function( const double, const double, const double ) )
{
    // Initialization
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::PIO::logOnlyNodeZero( log_file );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    std::cout << "testing with input file " << input_file << std::endl;
    std::cout.flush();

    // Test create
    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    // Get the Mesh database and create the mesh parameters
    auto database = input_db->getDatabase( "Mesh" );
    auto params   = std::make_shared<AMP::Mesh::MeshParameters>( database );
    params->setComm( globalComm );

    // Create the meshes from the input database
    auto meshAdapter = AMP::Mesh::Mesh::buildMesh( params );

    //----------------------------------------------------------------------------------------------------------------------------------------------//

    std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementModel;
    auto diffFEOp_db =
        std::dynamic_pointer_cast<AMP::Database>( input_db->getDatabase( "NonlinearDiffusionOp" ) );
    auto nonlinearOperator = AMP::Operator::OperatorBuilder::createOperator(
        meshAdapter, "NonlinearDiffusionOp", input_db, elementModel );
    auto diffOp =
        std::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>( nonlinearOperator );

    ut->passes( exeName + ": create" );
    std::cout.flush();

    // set up defaults for materials arguments and create transport model
    std::shared_ptr<AMP::Database> transportModel_db;
    if ( input_db->keyExists( "DiffusionTransportModel" ) )
        transportModel_db = input_db->getDatabase( "DiffusionTransportModel" );
    auto elementPhysicsModel =
        AMP::Operator::ElementPhysicsModelFactory::createElementPhysicsModel( transportModel_db );
    auto transportModel =
        std::dynamic_pointer_cast<AMP::Operator::DiffusionTransportModel>( elementPhysicsModel );

    double defTemp = transportModel_db->getWithDefault<double>( "Default_Temperature", 400.0 );
    double defConc = transportModel_db->getWithDefault<double>( "Default_Concentration", .33 );
    double defBurn = transportModel_db->getWithDefault<double>( "Default_Burnup", .5 );

    // create parameters
    auto diffOpParams =
        std::make_shared<AMP::Operator::DiffusionNonlinearFEOperatorParameters>( diffFEOp_db );

    // nullify vectors in parameters
    diffOpParams->d_FrozenTemperature.reset();
    diffOpParams->d_FrozenConcentration.reset();
    diffOpParams->d_FrozenBurnup.reset();

    // create vectors for parameters
    auto active_db = diffFEOp_db->getDatabase( "ActiveInputVariables" );
    auto tVar      = std::make_shared<AMP::LinearAlgebra::Variable>(
        active_db->getWithDefault<std::string>( "Temperature", "not_specified" ) );
    auto cVar = std::make_shared<AMP::LinearAlgebra::Variable>(
        active_db->getWithDefault<std::string>( "Concentration", "not_specified" ) );
    auto bVar = std::make_shared<AMP::LinearAlgebra::Variable>(
        active_db->getWithDefault<std::string>( "Burnup", "not_specified" ) );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // Create a DOF manager for a nodal vector
    int DOFsPerNode     = 1;
    int nodalGhostWidth = 1;
    bool split          = true;
    auto nodalDofMap    = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, nodalGhostWidth, DOFsPerNode, split );
    //----------------------------------------------------------------------------------------------------------------------------------------------//

    // create solution, rhs, and residual vectors
    auto tVec = AMP::LinearAlgebra::createVector( nodalDofMap, tVar );
    auto cVec = AMP::LinearAlgebra::createVector( nodalDofMap, cVar );
    auto bVec = AMP::LinearAlgebra::createVector( nodalDofMap, bVar );
    tVec->setToScalar( defTemp );
    cVec->setToScalar( defConc );
    bVec->setToScalar( defBurn );

    // set principal variable vector
    if ( diffOp->getPrincipalVariableId() == AMP::Operator::Diffusion::TEMPERATURE )
        diffOpParams->d_FrozenTemperature = tVec;
    if ( diffOp->getPrincipalVariableId() == AMP::Operator::Diffusion::CONCENTRATION )
        diffOpParams->d_FrozenConcentration = cVec;
    if ( diffOp->getPrincipalVariableId() == AMP::Operator::Diffusion::BURNUP )
        diffOpParams->d_FrozenBurnup = bVec;

    // set frozen vectors in parameters
    if ( diffFEOp_db->getWithDefault( "FreezeTemperature", false ) )
        diffOpParams->d_FrozenTemperature = tVec;
    if ( diffFEOp_db->getWithDefault( "FreezeConcentration", false ) )
        diffOpParams->d_FrozenConcentration = cVec;
    if ( diffFEOp_db->getWithDefault( "FreezeBurnup", false ) )
        diffOpParams->d_FrozenBurnup = bVec;

    // set transport model
    diffOpParams->d_transportModel = transportModel;

    // Reset with new parameters
    diffOp->reset( diffOpParams );

    // set  up variables for apply tests
    // auto diffOp->getInputVariable(diffOp->getPrincipalVariableId());
    auto diffSolVar = diffOp->getOutputVariable();

    auto diffRhsVar     = diffOp->getOutputVariable();
    auto diffResVar     = diffOp->getOutputVariable();
    auto nonPrincIds    = diffOp->getNonPrincipalVariableIds();
    auto numNonPrincIds = nonPrincIds.size();
    std::vector<AMP::LinearAlgebra::Variable::shared_ptr> nonPrincVars( numNonPrincIds );
    AMP::LinearAlgebra::Variable::shared_ptr inputVar = diffOp->getInputVariable();
    for ( size_t i = 0; i < numNonPrincIds; i++ ) {
        // nonPrincVars[i] = diffOp->getInputVariable(nonPrincIds[i]);
        nonPrincVars[i] = std::dynamic_pointer_cast<AMP::LinearAlgebra::MultiVariable>( inputVar )
                              ->getVariable( i );
    }

    // set up vectors for apply tests
    // std::string msgPrefix=exeName+": apply";
    auto diffSolVec = AMP::LinearAlgebra::createVector( nodalDofMap, diffSolVar );
    auto diffRhsVec = AMP::LinearAlgebra::createVector( nodalDofMap, diffRhsVar );
    auto diffResVec = AMP::LinearAlgebra::createVector( nodalDofMap, diffResVar );

    std::vector<AMP::LinearAlgebra::Vector::shared_ptr> nonPrincVecs( numNonPrincIds );
    for ( unsigned int i = 0; i < numNonPrincIds; i++ ) {
        nonPrincVecs[i] = AMP::LinearAlgebra::createVector( nodalDofMap, nonPrincVars[i] );
        if ( nonPrincIds[i] == AMP::Operator::Diffusion::TEMPERATURE )
            nonPrincVecs[i]->setToScalar( defTemp );
        if ( nonPrincIds[i] == AMP::Operator::Diffusion::CONCENTRATION )
            nonPrincVecs[i]->setToScalar( defConc );
        if ( nonPrincIds[i] == AMP::Operator::Diffusion::BURNUP )
            nonPrincVecs[i]->setToScalar( defBurn );
    }
    diffRhsVec->setToScalar( 0.0 );

    int zeroGhostWidth = 0;
    auto curNode       = meshAdapter->getIterator( AMP::Mesh::GeomType::Vertex, zeroGhostWidth );
    auto endNode       = curNode.end();

    for ( curNode = curNode.begin(); curNode != endNode; ++curNode ) {
        // double x = curNode->x();
        double x = ( curNode->coord() )[0];
        double y = ( curNode->coord() )[1];
        double z = ( curNode->coord() )[2];
        std::vector<size_t> dofs;
        nodalDofMap->getDOFs( curNode->globalID(), dofs );
        double fval = function( x, y, z );
        diffSolVec->setValueByGlobalID( dofs[0], fval );
    }
    diffSolVec->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );

    // Compute finite element operator
    diffOp->apply( diffSolVec, diffResVec );

    // Check that interior values are zero.
    double totalBnd = 0.;
    for ( size_t face = 1; face < 64; face++ ) {
        for ( curNode = curNode.begin(); curNode != endNode; ++curNode ) {
            std::vector<size_t> dofs;
            nodalDofMap->getDOFs( curNode->globalID(), dofs );
            double fval = diffResVec->getValueByGlobalID( dofs[0] );
            totalBnd += fabs( fval );
        }
    }
    double totalAll = 0.;
    for ( curNode = curNode.begin(); curNode != endNode; ++curNode ) {
        std::vector<size_t> dofs;
        nodalDofMap->getDOFs( curNode->globalID(), dofs );
        double fval = diffResVec->getValueByGlobalID( dofs[0] );
        totalAll += fabs( fval );
    }
    totalBnd = globalComm.sumReduce( totalBnd );
    totalAll = globalComm.sumReduce( totalAll );
    int rank = globalComm.getRank();
    if ( rank == 0 ) {
        std::cout << "******** All = " << totalAll << ", Bnd = " << totalBnd
                  << ", Err = " << totalAll - totalBnd << "********" << std::endl;
    }

    // write values in mathematica form
    int nranks = globalComm.getSize();
    if ( nranks == 1 ) {
        size_t nnodes        = meshAdapter->numLocalElements( AMP::Mesh::GeomType::Vertex );
        int proc             = globalComm.getRank();
        int nproc            = globalComm.getSize();
        std::string filename = "values-" + exeName;
        std::ofstream file( filename.c_str() );
        if ( proc == 0 ) {
            file << "values={"
                 << "\n";
        }
        AMP::Mesh::MeshIterator node =
            meshAdapter->getIterator( AMP::Mesh::GeomType::Vertex, zeroGhostWidth );

        int i = 0;
        for ( ; node != node.end(); ++node ) {
            double x = ( node->coord() )[0];
            double y = ( node->coord() )[1];
            double z = ( node->coord() )[2];

            int ii = i;
            i += 1;
            double rval = diffResVec->getValueByLocalID( ii );
            double fval = function( x, y, z );
            file << "{" << x << "," << y << "," << z << "," << rval << "," << fval << "}";
            if ( i < (int) nnodes - 1 )
                file << ",\n";
        }
        if ( proc < nproc - 1 ) {
            file << ",\n";
        } else {
            file << "}\n";
        }
    }

    ut->passes( "values-" + exeName );
}

int testNonlinearDiffusion_2( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;
    const int NUMFILES          = 1;
    std::string files[NUMFILES] = { "Diffusion-TUI-Thermal-1" };

    for ( auto &file : files )
        nonlinearTest( &ut, file, x_linear );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
