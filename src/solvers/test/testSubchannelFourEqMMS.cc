#include <string>
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "ampmesh/Mesh.h"
#include "boost/shared_ptr.hpp"
#include "utils/InputDatabase.h"
#include "utils/Utilities.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"
#include "utils/Database.h"
#include "vectors/Variable.h"
#include "vectors/Vector.h"
#include "vectors/SimpleVector.h"
#include "operators/subchannel/SubchannelFourEqNonlinearOperator.h"
#include "operators/subchannel/SubchannelConstants.h"
#include "operators/OperatorBuilder.h"
#include "solvers/ColumnSolver.h"
#include "solvers/PetscKrylovSolverParameters.h"
#include "solvers/PetscKrylovSolver.h"
#include "solvers/PetscSNESSolverParameters.h"
#include "solvers/PetscSNESSolver.h"
#include "solvers/trilinos/TrilinosMLSolver.h"

#include "ampmesh/SiloIO.h"
#include "vectors/VectorBuilder.h"
#include "discretization/simpleDOF_Manager.h"
#include "discretization/MultiDOF_Manager.h"
#include "ampmesh/StructuredMeshHelper.h"


// Function to get the linear heat generation rate
double getLinearHeatGeneration( double Q, double H, double z )
{
    const double pi = 3.141592653589793;
    return 0.5*pi*Q/H*sin(pi*z/H);
}


// Function to get the enthalpy solution
// Note: this is only an approximation that assumes incompressible water and no friction
double getSolutionEnthalpy( double Q, double H, double m, double hin, double z )
{
    const double pi = 3.141592653589793;
    return hin + 0.5*Q/m*(1.0-cos(pi*z/H));
}


// Function to get the pressure solution
// Note: this is only an approximation for an incompressible fluid with a fixed density
double getSolutionPressure( AMP::Database::shared_ptr db, double H, double Pout, double p, double z )
{
    if ( db->keyExists("Inlet_Pressure") )
        return Pout + (1.-z/H)*(db->getDouble("Inlet_Pressure")-Pout);
    else
        return Pout + (H-z)*9.80665*p;
}


void flowTest(AMP::UnitTest *ut, std::string exeName )
{
    std::string input_file = "input_" + exeName;
    std::string log_file = "output_" + exeName;
    std::string silo_name = exeName;
    AMP::PIO::logAllNodes(log_file);
    AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

    // Read the input file
    boost::shared_ptr<AMP::InputDatabase>  input_db ( new AMP::InputDatabase ( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile ( input_file , input_db );

//=============================================================================
// mesh and dof manager
//=============================================================================

    // Get the Mesh database and create the mesh parameters
    AMP_INSIST(input_db->keyExists("Mesh"), "Key ''Mesh'' is missing!");
    boost::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase( "Mesh" );
    boost::shared_ptr<AMP::Mesh::MeshParameters> meshParams(new AMP::Mesh::MeshParameters(mesh_db));
    meshParams->setComm(globalComm);

    // Create the meshes from the input database
    boost::shared_ptr<AMP::Mesh::Mesh> subchannelMesh = AMP::Mesh::Mesh::buildMesh(meshParams);
    AMP::Mesh::Mesh::shared_ptr xyFaceMesh;
    xyFaceMesh = subchannelMesh->Subset( AMP::Mesh::StructuredMeshHelper::getXYFaceIterator( subchannelMesh , 0 ) );

    // create DOF manager
    AMP::Discretization::DOFManager::shared_ptr subchannelDOFManager;
    if ( subchannelMesh.get() != NULL ) {
      AMP::Mesh::MeshIterator axialFaces0 = AMP::Mesh::StructuredMeshHelper::getXYFaceIterator(subchannelMesh,0);
      AMP::Mesh::MeshIterator axialFaces1 = AMP::Mesh::StructuredMeshHelper::getXYFaceIterator(subchannelMesh,1);
      AMP::Mesh::MeshIterator gapFaces0 = AMP::Mesh::Mesh::getIterator( AMP::Mesh::Union, 
        AMP::Mesh::StructuredMeshHelper::getXZFaceIterator(subchannelMesh,0), 
        AMP::Mesh::StructuredMeshHelper::getYZFaceIterator(subchannelMesh,0) );
      AMP::Mesh::MeshIterator gapFaces1 = AMP::Mesh::Mesh::getIterator( AMP::Mesh::Union, 
        AMP::Mesh::StructuredMeshHelper::getXZFaceIterator(subchannelMesh,1), 
        AMP::Mesh::StructuredMeshHelper::getYZFaceIterator(subchannelMesh,1) );
      std::vector<AMP::Discretization::DOFManager::shared_ptr> subchannelChildrenDOFManagers(2);
      subchannelChildrenDOFManagers[0] = AMP::Discretization::simpleDOFManager::create( subchannelMesh, axialFaces1, axialFaces0, 3 );
      subchannelChildrenDOFManagers[1] = AMP::Discretization::simpleDOFManager::create( subchannelMesh, gapFaces1, gapFaces0, 1 );
      subchannelDOFManager = AMP::Discretization::DOFManager::shared_ptr( 
        new AMP::Discretization::multiDOFManager( subchannelMesh->getComm(), subchannelChildrenDOFManagers ) );
    }

//=============================================================================
// physics model, parameters, and operator creation
//=============================================================================

    // get subchannel physics model
    boost::shared_ptr<AMP::Database> subchannelPhysics_db = input_db->getDatabase("SubchannelPhysicsModel");
    boost::shared_ptr<AMP::Operator::ElementPhysicsModelParameters> params( new AMP::Operator::ElementPhysicsModelParameters(subchannelPhysics_db));
    boost::shared_ptr<AMP::Operator::SubchannelPhysicsModel>  subchannelPhysicsModel (new AMP::Operator::SubchannelPhysicsModel(params));

    // Create the SubchannelOperatorParameters
    boost::shared_ptr<AMP::Database> nonlinearOperator_db = input_db->getDatabase("SubchannelFourEqNonlinearOperator");
    boost::shared_ptr<AMP::Operator::SubchannelOperatorParameters> subchannelOpParams(new AMP::Operator::SubchannelOperatorParameters( nonlinearOperator_db ));
    subchannelOpParams->d_Mesh = subchannelMesh ;
    subchannelOpParams->d_subchannelPhysicsModel = subchannelPhysicsModel;
    subchannelOpParams->clad_x = input_db->getDatabase("CladProperties")->getDoubleArray("x");
    subchannelOpParams->clad_y = input_db->getDatabase("CladProperties")->getDoubleArray("y");
    subchannelOpParams->clad_d = input_db->getDatabase("CladProperties")->getDoubleArray("d");

    // create nonlinear operator
    boost::shared_ptr<AMP::Operator::SubchannelFourEqNonlinearOperator> nonlinearOperator (new AMP::Operator::SubchannelFourEqNonlinearOperator(subchannelOpParams));
    // reset the nonlinear operator
    nonlinearOperator->reset(subchannelOpParams);

    // pass creation test
    ut->passes(exeName+": creation");
    std::cout.flush();

    // get input and output variables
    AMP::LinearAlgebra::Variable::shared_ptr inputVariable  = nonlinearOperator->getInputVariable();
    AMP::LinearAlgebra::Variable::shared_ptr outputVariable = nonlinearOperator->getOutputVariable();

    // create solution, rhs, and residual vectors
    AMP::LinearAlgebra::Vector::shared_ptr manufacturedVec = AMP::LinearAlgebra::createVector( subchannelDOFManager, inputVariable, true );
    AMP::LinearAlgebra::Vector::shared_ptr solVec = AMP::LinearAlgebra::createVector( subchannelDOFManager, inputVariable,  true );
    AMP::LinearAlgebra::Vector::shared_ptr rhsVec = AMP::LinearAlgebra::createVector( subchannelDOFManager, outputVariable, true );
    AMP::LinearAlgebra::Vector::shared_ptr resVec = AMP::LinearAlgebra::createVector( subchannelDOFManager, outputVariable, true );

//=============================================================================
// compute manufactured solution
//=============================================================================

    // Get the problem parameters
    std::vector<double> box = subchannelMesh->getBoundingBox();
    AMP_ASSERT(box[4]==0.0);
    double H = box[5]-box[4];
    double m_in = nonlinearOperator_db->getDouble("Inlet_Mass_Flow_Rate");
    double w_in = nonlinearOperator_db->getDouble("Inlet_Lateral_Flow_Rate");
    double Q = nonlinearOperator_db->getDouble("Max_Rod_Power");
    double Pout = nonlinearOperator_db->getDouble("Exit_Pressure");
    double Tin = nonlinearOperator_db->getDouble("Inlet_Temperature");

    // compute inlet enthalpy
    double Pin = Pout;
    double hin = 0.0;
    double rho_in = 1000;
    // iterate to find inlet pressure and inlet enthalpy
    for (int i=0; i<3; i++) {
       // compute inlet enthalpy using inlet temperature and outlet pressure
       std::map<std::string, boost::shared_ptr<std::vector<double> > > enthalpyArgMap;
       enthalpyArgMap.insert(std::make_pair("temperature",new std::vector<double>(1,Tin)));
       enthalpyArgMap.insert(std::make_pair("pressure",   new std::vector<double>(1,Pin)));
       std::vector<double> enthalpyResult(1);
       subchannelPhysicsModel->getProperty("Enthalpy",enthalpyResult,enthalpyArgMap); 
       hin = enthalpyResult[0];
       // compute inlet density using computed inlet enthalpy and outlet pressure
       std::map<std::string, boost::shared_ptr<std::vector<double> > > volumeArgMap_plus;
       volumeArgMap_plus.insert(std::make_pair("enthalpy",new std::vector<double>(1,hin)));
       volumeArgMap_plus.insert(std::make_pair("pressure",new std::vector<double>(1,Pin)));
       std::vector<double> volumeResult_plus(1);
       subchannelPhysicsModel->getProperty("SpecificVolume",volumeResult_plus,volumeArgMap_plus); 
       rho_in = 1.0/volumeResult_plus[0];
       // compute inlet pressure
       Pin = getSolutionPressure(input_db,H,Pout,rho_in,0);
    }
    std::cout<< "Inlet density:"<< rho_in <<std::endl;
    std::cout<< "Enthalpy Solution:"<< hin <<std::endl;

    // Compute the manufactured solution
    AMP::Mesh::MeshIterator face = xyFaceMesh->getIterator(AMP::Mesh::Face, 0);
    std::vector<size_t> axialDofs;
    // Scale to change the input vector back to correct units
    const double h_scale = 1.0/AMP::Operator::Subchannel::scaleEnthalpy;
    const double P_scale = 1.0/AMP::Operator::Subchannel::scalePressure;
    const double m_scale = 1.0/AMP::Operator::Subchannel::scaleAxialMassFlowRate;
    const double w_scale = 1.0/AMP::Operator::Subchannel::scaleLateralMassFlowRate;
    // loop over axial faces
    for (int i=0; i<(int)face.size(); i++){
        subchannelDOFManager->getDOFs( face->globalID(), axialDofs );
        std::vector<double> coord = face->centroid();
        double z = coord[2];
        double h = getSolutionEnthalpy( Q, H, m_in, hin, z );
        double P = getSolutionPressure( input_db, H, Pout, rho_in, z );
        manufacturedVec->setValueByGlobalID(axialDofs[0], m_in/m_scale);
        manufacturedVec->setValueByGlobalID(axialDofs[1], h/h_scale);
        manufacturedVec->setValueByGlobalID(axialDofs[2], P/P_scale);
        ++face;
    }
    // get lateral face map
    std::map<std::vector<double>,AMP::Mesh::MeshElement> lateralFaceMap = nonlinearOperator->getLateralFaces(subchannelOpParams->d_Mesh);
    // loop over lateral faces
    for (face = face.begin(); face != face.end(); face++) {
       std::vector<double> faceCentroid = face->centroid();
       std::map<std::vector<double>,AMP::Mesh::MeshElement>::iterator lateralFaceIterator = lateralFaceMap.find(faceCentroid);
       if (lateralFaceIterator != lateralFaceMap.end()) {
          // get lateral face
          AMP::Mesh::MeshElement lateralFace = lateralFaceIterator->second;
          // get crossflow from solution vector
          std::vector<size_t> gapDofs;
          subchannelDOFManager->getDOFs(lateralFace.globalID(),gapDofs);
          double w = 0.0;
          manufacturedVec->setValueByGlobalID(gapDofs[0], w/w_scale);
       }
    }

//=============================================================================
// compute initial guess
//=============================================================================

    // Compute the initial guess solution
    face = xyFaceMesh->getIterator(AMP::Mesh::Face, 0);
    // loop over axial faces
    for (int i=0; i<(int)face.size(); i++){
        subchannelDOFManager->getDOFs( face->globalID(), axialDofs );
        solVec->setValueByGlobalID(axialDofs[0], m_in/m_scale);
        solVec->setValueByGlobalID(axialDofs[1], hin/h_scale);
        solVec->setValueByGlobalID(axialDofs[2], Pout/P_scale);
        ++face;
    }
    // loop over lateral faces
    for (face = face.begin(); face != face.end(); face++) {
       std::vector<double> faceCentroid = face->centroid();
       std::map<std::vector<double>,AMP::Mesh::MeshElement>::iterator lateralFaceIterator = lateralFaceMap.find(faceCentroid);
       if (lateralFaceIterator != lateralFaceMap.end()) {
          // get lateral face
          AMP::Mesh::MeshElement lateralFace = lateralFaceIterator->second;
          // get crossflow from solution vector
          std::vector<size_t> gapDofs;
          subchannelDOFManager->getDOFs(lateralFace.globalID(),gapDofs);
          solVec->setValueByGlobalID(gapDofs[0], w_in/w_scale);
       }
    }
    solVec->copyVector(manufacturedVec);

//=============================================================================
// solve
//=============================================================================

    // get nonlinear solver database
    boost::shared_ptr<AMP::Database> nonlinearSolver_db = input_db->getDatabase("NonlinearSolver"); 
  
    // get linear solver database
//    boost::shared_ptr<AMP::Database> linearSolver_db = nonlinearSolver_db->getDatabase("LinearSolver"); 
 
    // put manufactured RHS into resVec
//    linearOperator->reset(nonlinearOperator->getJacobianParameters(solVec));
//    linearOperator->apply(rhsVec, solVec, resVec, 1.0, -1.0);
   
    // create nonlinear solver parameters
    boost::shared_ptr<AMP::Solver::PetscSNESSolverParameters> nonlinearSolverParams(new AMP::Solver::PetscSNESSolverParameters(nonlinearSolver_db));

    // change the next line to get the correct communicator out
    nonlinearSolverParams->d_comm = globalComm;
    nonlinearSolverParams->d_pOperator = nonlinearOperator;
    nonlinearSolverParams->d_pInitialGuess = solVec;

    // create nonlinear solver
    boost::shared_ptr<AMP::Solver::PetscSNESSolver> nonlinearSolver(new AMP::Solver::PetscSNESSolver(nonlinearSolverParams));

    // create linear solver
//    boost::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver = nonlinearSolver->getKrylovSolver();

/*
    // create preconditioner
    boost::shared_ptr<AMP::Database> Preconditioner_db =  linearSolver_db->getDatabase("Preconditioner");
    boost::shared_ptr<AMP::Solver::SolverStrategyParameters> PreconditionerParams(new AMP::Solver::SolverStrategyParameters(Preconditioner_db));
    PreconditionerParams->d_pOperator = linearOperator;
    boost::shared_ptr<AMP::Solver::TrilinosMLSolver> linearFlowPreconditioner(new AMP::Solver::TrilinosMLSolver(PreconditionerParams));
    // set preconditioner
    linearSolver->setPreconditioner(linearFlowPreconditioner);
*/

    // don't use zero initial guess
    nonlinearSolver->setZeroInitialGuess(false);

    // solve
    nonlinearSolver->solve(rhsVec, solVec);
    nonlinearOperator->apply(rhsVec, solVec, resVec, 1.0, -1.0);

//=============================================================================
// examine solution
//=============================================================================

/*
    // Compute the flow temperature
    AMP::Discretization::DOFManager::shared_ptr tempDOFManager = AMP::Discretization::simpleDOFManager::create( subchannelMesh, 
        AMP::Mesh::StructuredMeshHelper::getXYFaceIterator(subchannelMesh,1), AMP::Mesh::StructuredMeshHelper::getXYFaceIterator(subchannelMesh,0), 1 );
    AMP::LinearAlgebra::Variable::shared_ptr  tempVariable( new AMP::LinearAlgebra::Variable("Temperature") );
    AMP::LinearAlgebra::Vector::shared_ptr tempVec = AMP::LinearAlgebra::createVector( tempDOFManager , tempVariable  , true );
    face  = xyFaceMesh->getIterator(AMP::Mesh::Face, 0);
    std::vector<size_t> tdofs;
    bool pass = true;
    for (int i=0; i<(int)face.size(); i++){
        subchannelDOFManager->getDOFs( face->globalID(), axialDofs );
        tempDOFManager->getDOFs( face->globalID(), tdofs );
        double h = h_scale*solVec->getValueByGlobalID(axialDofs[0]);
        double P = P_scale*solVec->getValueByGlobalID(axialDofs[1]);
        std::map<std::string, boost::shared_ptr<std::vector<double> > > temperatureArgMap;
        temperatureArgMap.insert(std::make_pair("enthalpy",new std::vector<double>(1,h)));
        temperatureArgMap.insert(std::make_pair("pressure",new std::vector<double>(1,P)));
        std::vector<double> temperatureResult(1);
        subchannelPhysicsModel->getProperty("Temperature", temperatureResult, temperatureArgMap); 
        tempVec->setValueByGlobalID(tdofs[0],temperatureResult[0]);
        // Check that we recover the enthalpy from the temperature
        std::map<std::string, boost::shared_ptr<std::vector<double> > > enthalpyArgMap;
        enthalpyArgMap.insert(std::make_pair("temperature",new std::vector<double>(1,temperatureResult[0])));
        enthalpyArgMap.insert(std::make_pair("pressure",   new std::vector<double>(1,P)));
        std::vector<double> enthalpyResult(1);
        subchannelPhysicsModel->getProperty("Enthalpy",enthalpyResult,enthalpyArgMap); 
        double h2 = enthalpyResult[0];
        if ( !AMP::Utilities::approx_equal(h,h2,1e-7) )
            pass = false;
        ++face;
    } 
    if ( !pass )
        ut->failure("failed to recover h");

    // Print the Inlet/Outlet properties
    std::cout << std::endl << std::endl;
    face  = xyFaceMesh->getIterator(AMP::Mesh::Face, 0);
    subchannelDOFManager->getDOFs( face->globalID(), axialDofs );
    tempDOFManager->getDOFs( face->globalID(), tdofs );
    std::cout<< "Inlet Computed Mass Flow Rate = " << m_scale*solVec->getValueByGlobalID(axialDofs[0]) << std::endl;
    std::cout<< "Inlet Computed Enthalpy = " << h_scale*solVec->getValueByGlobalID(axialDofs[1]) << std::endl;
    std::cout<< "Inlet Computed Pressure = " << P_scale*solVec->getValueByGlobalID(axialDofs[2]) << std::endl;
    std::cout<< "Inlet Computed Temperature = " << tempVec->getValueByGlobalID(tdofs[0]) << std::endl;
    std::cout << std::endl;
    face = --((xyFaceMesh->getIterator(AMP::Mesh::Face,0)).end());
    subchannelDOFManager->getDOFs( face->globalID(), axialDofs );
    tempDOFManager->getDOFs( face->globalID(), tdofs );
    std::cout<< "Outlet Computed Mass Flow Rate = " << m_scale*solVec->getValueByGlobalID(axialDofs[0]) << std::endl;
    std::cout<< "Outlet Computed Enthalpy = " << h_scale*solVec->getValueByGlobalID(axialDofs[1]) << std::endl;
    std::cout<< "Outlet Computed Pressure = " << P_scale*solVec->getValueByGlobalID(axialDofs[2]) << std::endl;
    std::cout<< "Outlet Computed Temperature = " << tempVec->getValueByGlobalID(tdofs[0]) << std::endl;

    // Compute the error
    AMP::LinearAlgebra::Vector::shared_ptr absErrorVec = solVec->cloneVector();
    absErrorVec->axpy(-1.0,solVec,manufacturedVec);
    AMP::LinearAlgebra::Vector::shared_ptr relErrorVec = solVec->cloneVector();
    relErrorVec->divide(absErrorVec,manufacturedVec);
    
    double absErrorNorm = absErrorVec->L2Norm();
    double relErrorNorm = relErrorVec->L2Norm();

    // check that norm of relative error is less than tolerance
    double tol = input_db->getDoubleWithDefault("TOLERANCE",1e-6);
    if(relErrorNorm > tol){
        ut->failure(exeName+": manufactured solution test");
    } else {
        ut->passes(exeName+": manufactured solution test");
    }

    // Print final solution
    face  = xyFaceMesh->getIterator(AMP::Mesh::Face, 0);
    std::cout<<std::endl;
    int N_print = std::max(1,(int)face.size()/10);
    for (int i=0; i<(int)face.size(); i++){
        if ( i%N_print==0 ) {
            subchannelDOFManager->getDOFs( face->globalID(), axialDofs );
            std::cout<< "Computed Mass Flow Rate["<<i<<"] = "<< m_scale*solVec->getValueByGlobalID(axialDofs[0]) << std::endl;
            std::cout<< "Solution Mass Flow Rate["<<i<<"] = "<< m_scale*manufacturedVec->getValueByGlobalID(axialDofs[0]) << std::endl;
            std::cout<< "Computed Enthalpy["<<i<<"] = "<< h_scale*solVec->getValueByGlobalID(axialDofs[1]) << std::endl;
            std::cout<< "Solution Enthalpy["<<i<<"] = "<< h_scale*manufacturedVec->getValueByGlobalID(axialDofs[1]) << std::endl;
            std::cout<< "Computed Pressure["<<i<<"] = "<< P_scale*solVec->getValueByGlobalID(axialDofs[2]) << std::endl;
            std::cout<< "Solution Pressure["<<i<<"] = "<< P_scale*manufacturedVec->getValueByGlobalID(axialDofs[2]) << std::endl;
            std::cout<<std::endl;
        }
        ++face;
    }
    std::cout<<"L2 Norm of Absolute Error: "<<absErrorNorm<<std::endl;
    std::cout<<"L2 Norm of Relative Error: "<<relErrorNorm<<std::endl;

    input_db.reset();

#ifdef USE_EXT_SILO
    // Rescale the solution to get the correct units
    AMP::LinearAlgebra::Vector::shared_ptr mass, enthalpy, pressure;
    mass = solVec->select( AMP::LinearAlgebra::VS_Stride(0,3), "M" );
    enthalpy = solVec->select( AMP::LinearAlgebra::VS_Stride(1,3), "H" );
    pressure = solVec->select( AMP::LinearAlgebra::VS_Stride(2,3), "P" );
    mass->scale(m_scale);
    enthalpy->scale(h_scale);
    pressure->scale(P_scale);
    mass = manufacturedVec->select( AMP::LinearAlgebra::VS_Stride(0,3), "M" );
    enthalpy = manufacturedVec->select( AMP::LinearAlgebra::VS_Stride(1,3), "H" );
    pressure = manufacturedVec->select( AMP::LinearAlgebra::VS_Stride(2,3), "P" );
    mass->scale(m_scale);
    enthalpy->scale(h_scale);
    pressure->scale(P_scale);
    // Register the quantities to plot
    AMP::Mesh::SiloIO::shared_ptr  siloWriter( new AMP::Mesh::SiloIO );
    AMP::LinearAlgebra::Vector::shared_ptr subchannelMass = solVec->select( AMP::LinearAlgebra::VS_Stride(0,3), "M" );
    AMP::LinearAlgebra::Vector::shared_ptr subchannelEnthalpy = solVec->select( AMP::LinearAlgebra::VS_Stride(1,3), "H" );
    AMP::LinearAlgebra::Vector::shared_ptr subchannelPressure = solVec->select( AMP::LinearAlgebra::VS_Stride(2,3), "P" );
    subchannelMass->scale(m_scale);
    subchannelEnthalpy->scale(h_scale);
    subchannelPressure->scale(P_scale);
    siloWriter->registerVector( manufacturedVec, xyFaceMesh, AMP::Mesh::Face, "ManufacturedSolution" );
    siloWriter->registerVector( solVec, xyFaceMesh, AMP::Mesh::Face, "ComputedSolution" );
    siloWriter->registerVector( subchannelMass, xyFaceMesh, AMP::Mesh::Face, "Axial Mass Flow Rate" );
    siloWriter->registerVector( subchannelEnthalpy, xyFaceMesh, AMP::Mesh::Face, "Enthalpy" );
    siloWriter->registerVector( subchannelPressure, xyFaceMesh, AMP::Mesh::Face, "Pressure" );
    siloWriter->registerVector( tempVec, xyFaceMesh, AMP::Mesh::Face, "Temperature" );
    siloWriter->writeFile( silo_name , 0 );
#endif
*/


}

int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

    std::vector<std::string> files(1);
    files[0] = "testSubchannelFourEqMMS-1";

    for (size_t i=0; i<files.size(); i++)
        flowTest(&ut,files[i]);

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}   

