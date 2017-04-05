// This test checks the verification problem in SubChannelFlow.tex
#include "ampmesh/Mesh.h"
#include "operators/IdentityOperator.h"
#include "operators/OperatorBuilder.h"
#include "operators/subchannel/SubchannelConstants.h"
#include "operators/subchannel/SubchannelTwoEqLinearOperator.h"
#include "operators/subchannel/SubchannelTwoEqNonlinearOperator.h"
#include "solvers/ColumnSolver.h"
#include "solvers/trilinos/TrilinosMLSolver.h"
#include "solvers/trilinos/nox/TrilinosNOXSolver.h"
#include "solvers/trilinos/nox/TrilinosNOXSolverParameters.h"
#include "utils/AMPManager.h"
#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "utils/Utilities.h"
#include "utils/shared_ptr.h"
#include "vectors/SimpleVector.h"
#include "vectors/Variable.h"
#include "vectors/Vector.h"
#include <string>

#include "ampmesh/StructuredMeshHelper.h"
#include "discretization/simpleDOF_Manager.h"
#include "discretization/structuredFaceDOFManager.h"
#include "utils/Writer.h"
#include "vectors/VectorBuilder.h"


typedef std::vector<double> doubleVec;


// Function to get the linear heat generation rate
double getLinearHeatGeneration( double Q, double H, double z )
{
    const double pi = 3.141592653589793;
    return 0.5 * pi * Q / H * sin( pi * z / H );
}


// Function to get the enthalpy solution
// Note: this is only an approximation that assumes incompressible water and no friction
double getSolutionEnthalpy( double Q, double H, double m, double hin, double z )
{
    const double pi = 3.141592653589793;
    return hin + 0.5 * Q / m * ( 1.0 - cos( pi * z / H ) );
}


// Function to get the pressure solution
// Note: this is only an approximation for an incompressible fluid with a fixed density
double
getSolutionPressure( AMP::Database::shared_ptr db, double H, double Pout, double p, double z )
{
    if ( db->keyExists( "Inlet_Pressure" ) )
        return Pout + ( 1. - z / H ) * ( db->getDouble( "Inlet_Pressure" ) - Pout );
    else
        return Pout + ( H - z ) * 9.80665 * p;
}


void flowTest( AMP::UnitTest *ut, std::string exeName )
{
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;
    AMP::PIO::logAllNodes( log_file );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    // Read the input file
    AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );

    // Get the Mesh database and create the mesh parameters
    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    AMP::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase( "Mesh" );
    AMP::shared_ptr<AMP::Mesh::MeshParameters> meshParams(
        new AMP::Mesh::MeshParameters( mesh_db ) );
    meshParams->setComm( globalComm );

    // Create the meshes from the input database
    AMP::shared_ptr<AMP::Mesh::Mesh> subchannelMesh = AMP::Mesh::Mesh::buildMesh( meshParams );

    // get subchannel physics model
    AMP::shared_ptr<AMP::Database> subchannelPhysics_db =
        input_db->getDatabase( "SubchannelPhysicsModel" );
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModelParameters> params(
        new AMP::Operator::ElementPhysicsModelParameters( subchannelPhysics_db ) );
    AMP::shared_ptr<AMP::Operator::SubchannelPhysicsModel> subchannelPhysicsModel(
        new AMP::Operator::SubchannelPhysicsModel( params ) );

    // Create the SubchannelOperatorParameters
    AMP::shared_ptr<AMP::Database> nonlinearOperator_db =
        input_db->getDatabase( "SubchannelTwoEqNonlinearOperator" );
    AMP::shared_ptr<AMP::Operator::SubchannelOperatorParameters> subchannelOpParams(
        new AMP::Operator::SubchannelOperatorParameters( nonlinearOperator_db ) );
    subchannelOpParams->d_Mesh                   = subchannelMesh;
    subchannelOpParams->d_subchannelPhysicsModel = subchannelPhysicsModel;
    subchannelOpParams->clad_x = input_db->getDatabase( "CladProperties" )->getDoubleArray( "x" );
    subchannelOpParams->clad_y = input_db->getDatabase( "CladProperties" )->getDoubleArray( "y" );
    subchannelOpParams->clad_d = input_db->getDatabase( "CladProperties" )->getDoubleArray( "d" );

    // create nonlinear operator
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> elementModel;
    AMP::shared_ptr<AMP::Operator::SubchannelTwoEqNonlinearOperator> nonlinearOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::SubchannelTwoEqNonlinearOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                subchannelMesh, "SubchannelTwoEqNonlinearOperator", input_db, elementModel ) );

    // create linear operator
    AMP::shared_ptr<AMP::Operator::LinearOperator> linearOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::LinearOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                subchannelMesh, "SubchannelTwoEqLinearOperator", input_db, elementModel ) );

    // pass creation test
    ut->passes( exeName + ": creation" );
    std::cout.flush();

    // get input and output variables
    AMP::LinearAlgebra::Variable::shared_ptr inputVariable = nonlinearOperator->getInputVariable();
    AMP::LinearAlgebra::Variable::shared_ptr outputVariable =
        nonlinearOperator->getOutputVariable();

    // get dof manager
    int DOFsPerFace[3] = { 0, 0, 2 };
    AMP::Discretization::DOFManager::shared_ptr faceDOFManager =
        AMP::Discretization::structuredFaceDOFManager::create( subchannelMesh, DOFsPerFace, 0 );

    // create solution, rhs, and residual vectors
    AMP::LinearAlgebra::Vector::shared_ptr manufacturedVec =
        AMP::LinearAlgebra::createVector( faceDOFManager, inputVariable, true );
    AMP::LinearAlgebra::Vector::shared_ptr solVec =
        AMP::LinearAlgebra::createVector( faceDOFManager, inputVariable, true );
    AMP::LinearAlgebra::Vector::shared_ptr rhsVec =
        AMP::LinearAlgebra::createVector( faceDOFManager, outputVariable, true );
    AMP::LinearAlgebra::Vector::shared_ptr resVec =
        AMP::LinearAlgebra::createVector( faceDOFManager, outputVariable, true );

    // Get the problem parameters
    std::vector<double> box = subchannelMesh->getBoundingBox();
    AMP_ASSERT( box[4] == 0.0 );
    double H    = box[5] - box[4];
    double m    = nonlinearOperator_db->getDouble( "Inlet_Mass_Flow_Rate" );
    double Q    = nonlinearOperator_db->getDouble( "Rod_Power" );
    double Pout = nonlinearOperator_db->getDouble( "Exit_Pressure" );
    double Tin  = nonlinearOperator_db->getDouble( "Inlet_Temperature" );

    // compute inlet enthalpy
    double Pin = Pout;
    double hin = 0.0;
    double rho = 1000;
    for ( int i = 0; i < 3; i++ ) {
        std::map<std::string, AMP::shared_ptr<std::vector<double>>> enthalpyArgMap;
        enthalpyArgMap.insert( std::make_pair( "temperature", AMP::make_shared<doubleVec>( 1, Tin ) ) );
        enthalpyArgMap.insert( std::make_pair( "pressure", AMP::make_shared<doubleVec>( 1, Pin ) ) );
        std::vector<double> enthalpyResult( 1 );
        subchannelPhysicsModel->getProperty( "Enthalpy", enthalpyResult, enthalpyArgMap );
        hin = enthalpyResult[0];
        std::map<std::string, AMP::shared_ptr<std::vector<double>>> volumeArgMap_plus;
        volumeArgMap_plus.insert( std::make_pair( "enthalpy", AMP::make_shared<doubleVec>( 1, hin ) ) );
        volumeArgMap_plus.insert( std::make_pair( "pressure", AMP::make_shared<doubleVec>( 1, Pin ) ) );
        std::vector<double> volumeResult_plus( 1 );
        subchannelPhysicsModel->getProperty(
            "SpecificGeomType::Volume", volumeResult_plus, volumeArgMap_plus );
        rho = 1.0 / volumeResult_plus[0];
        Pin = getSolutionPressure( input_db, H, Pout, rho, 0 );
    }
    std::cout << "Inlet density:" << rho << std::endl;
    std::cout << "Enthalpy Solution:" << hin << std::endl;

    // Compute the manufactured solution
    AMP::Mesh::Mesh::shared_ptr xyFaceMesh = subchannelMesh->Subset(
        AMP::Mesh::StructuredMeshHelper::getXYFaceIterator( subchannelMesh, 0 ) );
    AMP::Mesh::MeshIterator face = xyFaceMesh->getIterator( AMP::Mesh::GeomType::Face, 0 );
    std::vector<size_t> dofs;
    const double h_scale = 1.0 / AMP::Operator::Subchannel::scaleEnthalpy; // Scale to change the
                                                                           // input vector back to
                                                                           // correct units
    const double P_scale = 1.0 / AMP::Operator::Subchannel::scalePressure; // Scale to change the
                                                                           // input vector back to
                                                                           // correct units
    for ( int i = 0; i < (int) face.size(); i++ ) {
        faceDOFManager->getDOFs( face->globalID(), dofs );
        std::vector<double> coord = face->centroid();
        double z                  = coord[2];
        double h                  = getSolutionEnthalpy( Q, H, m, hin, z );
        double P                  = getSolutionPressure( input_db, H, Pout, rho, z );
        manufacturedVec->setValueByGlobalID( dofs[0],
                                             AMP::Operator::Subchannel::scaleEnthalpy * h );
        manufacturedVec->setValueByGlobalID( dofs[1],
                                             AMP::Operator::Subchannel::scalePressure * P );
        ++face;
    }

    // Compute the initial guess solution
    face = xyFaceMesh->getIterator( AMP::Mesh::GeomType::Face, 0 );
    for ( int i = 0; i < (int) face.size(); i++ ) {
        faceDOFManager->getDOFs( face->globalID(), dofs );
        solVec->setValueByGlobalID( dofs[0], AMP::Operator::Subchannel::scaleEnthalpy * hin );
        solVec->setValueByGlobalID( dofs[1], AMP::Operator::Subchannel::scalePressure * Pout );
        ++face;
    }
    solVec->copyVector( manufacturedVec );

    // get nonlinear solver database
    AMP::shared_ptr<AMP::Database> nonlinearSolver_db = input_db->getDatabase( "NonlinearSolver" );

    // get linear solver database
    AMP::shared_ptr<AMP::Database> linearSolver_db =
        nonlinearSolver_db->getDatabase( "LinearSolver" );

    // put manufactured RHS into resVec
    nonlinearOperator->reset( subchannelOpParams );
    AMP::shared_ptr<AMP::Operator::SubchannelOperatorParameters> subchannelLinearParams =
        AMP::dynamic_pointer_cast<AMP::Operator::SubchannelOperatorParameters>(
            nonlinearOperator->getParameters( "Jacobian", solVec ) );
    subchannelLinearParams->d_initialize = false;
    linearOperator->reset( subchannelLinearParams );
    linearOperator->residual( rhsVec, solVec, resVec );

    // create nonlinear solver parameters
    AMP::shared_ptr<AMP::Solver::TrilinosNOXSolverParameters> nonlinearSolverParams(
        new AMP::Solver::TrilinosNOXSolverParameters( nonlinearSolver_db ) );

    // change the next line to get the correct communicator out
    nonlinearSolverParams->d_comm            = globalComm;
    nonlinearSolverParams->d_pOperator       = nonlinearOperator;
    nonlinearSolverParams->d_pLinearOperator = nonlinearOperator;
    nonlinearSolverParams->d_pInitialGuess   = solVec;

    // create nonlinear solver
    AMP::shared_ptr<AMP::Solver::TrilinosNOXSolver> nonlinearSolver(
        new AMP::Solver::TrilinosNOXSolver( nonlinearSolverParams ) );

    // create linear solver
    // nonlinearSolver->getKrylovSolver()->setPreconditioner(linearFlowPreconditioner);


    // don't use zero initial guess
    nonlinearSolver->setZeroInitialGuess( false );

    // solve
    nonlinearSolver->solve( rhsVec, solVec );
    nonlinearOperator->residual( rhsVec, solVec, resVec );

    // Compute the flow temperature
    int tempDOFsPerFace[3] = { 0, 0, 1 };
    AMP::Discretization::DOFManager::shared_ptr tempDOFManager =
        AMP::Discretization::structuredFaceDOFManager::create( subchannelMesh, tempDOFsPerFace, 0 );
    AMP::LinearAlgebra::Variable::shared_ptr tempVariable(
        new AMP::LinearAlgebra::Variable( "Temperature" ) );
    AMP::LinearAlgebra::Vector::shared_ptr tempVec =
        AMP::LinearAlgebra::createVector( tempDOFManager, tempVariable, true );
    face = xyFaceMesh->getIterator( AMP::Mesh::GeomType::Face, 0 );
    std::vector<size_t> tdofs;
    bool pass = true;
    for ( int i = 0; i < (int) face.size(); i++ ) {
        faceDOFManager->getDOFs( face->globalID(), dofs );
        tempDOFManager->getDOFs( face->globalID(), tdofs );
        double h = h_scale * solVec->getValueByGlobalID( dofs[0] );
        double P = P_scale * solVec->getValueByGlobalID( dofs[1] );
        std::map<std::string, AMP::shared_ptr<std::vector<double>>> temperatureArgMap;
        temperatureArgMap.insert( std::make_pair( "enthalpy", AMP::make_shared<doubleVec>( 1, h ) ) );
        temperatureArgMap.insert( std::make_pair( "pressure", AMP::make_shared<doubleVec>( 1, P ) ) );
        std::vector<double> temperatureResult( 1 );
        subchannelPhysicsModel->getProperty( "Temperature", temperatureResult, temperatureArgMap );
        tempVec->setValueByGlobalID( tdofs[0], temperatureResult[0] );
        // Check that we recover the enthalapy from the temperature
        std::map<std::string, AMP::shared_ptr<std::vector<double>>> enthalpyArgMap;
        enthalpyArgMap.insert(
            std::make_pair( "temperature", AMP::make_shared<doubleVec>( 1, temperatureResult[0] ) ) );
        enthalpyArgMap.insert( std::make_pair( "pressure", AMP::make_shared<doubleVec>( 1, P ) ) );
        std::vector<double> enthalpyResult( 1 );
        subchannelPhysicsModel->getProperty( "Enthalpy", enthalpyResult, enthalpyArgMap );
        double h2 = enthalpyResult[0];
        if ( !AMP::Utilities::approx_equal( h, h2, 1e-7 ) )
            pass = false;
        ++face;
    }
    if ( !pass )
        ut->failure( "failed to recover h" );

    // Print the Inlet/Outlet properties
    std::cout << std::endl << std::endl;
    face = xyFaceMesh->getIterator( AMP::Mesh::GeomType::Face, 0 );
    faceDOFManager->getDOFs( face->globalID(), dofs );
    tempDOFManager->getDOFs( face->globalID(), tdofs );
    double TinSol = tempVec->getValueByGlobalID( tdofs[0] );
    std::cout << "Inlet Computed Enthalpy = " << h_scale * solVec->getValueByGlobalID( dofs[0] )
              << std::endl;
    std::cout << "Inlet Computed Pressure = " << P_scale * solVec->getValueByGlobalID( dofs[1] )
              << std::endl;
    std::cout << "Inlet Computed Temperature = " << TinSol << std::endl;
    std::cout << std::endl;
    face = --( ( xyFaceMesh->getIterator( AMP::Mesh::GeomType::Face, 0 ) ).end() );
    faceDOFManager->getDOFs( face->globalID(), dofs );
    tempDOFManager->getDOFs( face->globalID(), tdofs );
    double ToutSol = tempVec->getValueByGlobalID( tdofs[0] );
    std::cout << "Outlet Computed Enthalpy = " << h_scale * solVec->getValueByGlobalID( dofs[0] )
              << std::endl;
    std::cout << "Outlet Computed Pressure = " << P_scale * solVec->getValueByGlobalID( dofs[1] )
              << std::endl;
    std::cout << "Outlet Computed Temperature = " << ToutSol << std::endl;

    // Compute the error
    AMP::LinearAlgebra::Vector::shared_ptr absErrorVec = solVec->cloneVector();
    absErrorVec->axpy( -1.0, solVec, manufacturedVec );
    AMP::LinearAlgebra::Vector::shared_ptr relErrorVec = solVec->cloneVector();
    relErrorVec->divide( absErrorVec, manufacturedVec );
    /*face  = xyFaceMesh->getIterator(AMP::Mesh::GeomType::Face, 0);
    for (int i=0; i<(int)face.size(); i++){
        faceDOFManager->getDOFs( face->globalID(), dofs );
        absErrorVec->setValueByGlobalID(dofs[1],0.0);   // We don't have the correct solution for
    the pressure yet
        relErrorVec->setValueByGlobalID(dofs[1],0.0);
        ++face;
    }*/
    double absErrorNorm = absErrorVec->L2Norm();
    double relErrorNorm = relErrorVec->L2Norm();

    // check that norm of relative error is less than tolerance
    double tol = input_db->getDoubleWithDefault( "TOLERANCE", 1e-6 );
    if ( relErrorNorm <= tol && fabs( Tin - TinSol ) < tol ) {
        ut->passes( exeName + ": manufactured solution test" );
    } else {
        ut->failure( exeName + ": manufactured solution test" );
    }

    // Print final solution
    face = xyFaceMesh->getIterator( AMP::Mesh::GeomType::Face, 0 );
    std::cout << std::endl;
    int N_print = std::max( 1, (int) face.size() / 10 );
    for ( int i = 0; i < (int) face.size(); i++ ) {
        if ( i % N_print == 0 ) {
            faceDOFManager->getDOFs( face->globalID(), dofs );
            std::cout << "Computed Enthalpy[" << i
                      << "] = " << h_scale * solVec->getValueByGlobalID( dofs[0] ) << std::endl;
            std::cout << "Solution Enthalpy[" << i
                      << "] = " << h_scale * manufacturedVec->getValueByGlobalID( dofs[0] )
                      << std::endl;
            std::cout << "Computed Pressure[" << i
                      << "] = " << P_scale * solVec->getValueByGlobalID( dofs[1] ) << std::endl;
            std::cout << "Solution Pressure[" << i
                      << "] = " << P_scale * manufacturedVec->getValueByGlobalID( dofs[1] )
                      << std::endl;
            std::cout << std::endl;
        }
        ++face;
    }
    std::cout << "L2 Norm of Absolute Error: " << absErrorNorm << std::endl;
    std::cout << "L2 Norm of Relative Error: " << relErrorNorm << std::endl;

    input_db.reset();

#ifdef USE_EXT_SILO
    // Rescale the solution to get the correct units
    AMP::LinearAlgebra::Vector::shared_ptr enthalpy, pressure;
    enthalpy = solVec->select( AMP::LinearAlgebra::VS_Stride( 0, 2 ), "H" );
    pressure = solVec->select( AMP::LinearAlgebra::VS_Stride( 1, 2 ), "P" );
    enthalpy->scale( h_scale );
    pressure->scale( P_scale );
    enthalpy = manufacturedVec->select( AMP::LinearAlgebra::VS_Stride( 0, 2 ), "H" );
    pressure = manufacturedVec->select( AMP::LinearAlgebra::VS_Stride( 1, 2 ), "P" );
    enthalpy->scale( h_scale );
    pressure->scale( P_scale );
    // Register the quantities to plot
    AMP::Utilities::Writer::shared_ptr siloWriter = AMP::Utilities::Writer::buildWriter( "Silo" );
    AMP::LinearAlgebra::Vector::shared_ptr subchannelEnthalpy =
        solVec->select( AMP::LinearAlgebra::VS_Stride( 0, 2 ), "H" );
    AMP::LinearAlgebra::Vector::shared_ptr subchannelPressure =
        solVec->select( AMP::LinearAlgebra::VS_Stride( 1, 2 ), "P" );
    subchannelEnthalpy->scale( h_scale );
    subchannelPressure->scale( P_scale );
    siloWriter->registerVector(
        manufacturedVec, xyFaceMesh, AMP::Mesh::GeomType::Face, "ManufacturedSolution" );
    siloWriter->registerVector( solVec, xyFaceMesh, AMP::Mesh::GeomType::Face, "ComputedSolution" );
    siloWriter->registerVector( subchannelEnthalpy, xyFaceMesh, AMP::Mesh::GeomType::Face, "Enthalpy" );
    siloWriter->registerVector( subchannelPressure, xyFaceMesh, AMP::Mesh::GeomType::Face, "Pressure" );
    siloWriter->registerVector( tempVec, xyFaceMesh, AMP::Mesh::GeomType::Face, "Temperature" );
    siloWriter->writeFile( exeName, 0 );
#endif
}

int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    std::vector<std::string> files( 2 );
    files[0] = "testSubchannelSolution-1";
    files[1] = "testSubchannelSolution-2";

    for ( size_t i = 0; i < files.size(); i++ )
        flowTest( &ut, files[i] );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
