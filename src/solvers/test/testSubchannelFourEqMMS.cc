#include <string>
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "ampmesh/Mesh.h"
#include "utils/shared_ptr.h"
#include "utils/InputDatabase.h"
#include "utils/Utilities.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"
#include "utils/Database.h"
#include "vectors/Variable.h"
#include "vectors/Vector.h"
#include "vectors/SimpleVector.h"
#include "operators/subchannel/SubchannelFourEqNonlinearOperator.h"
#include "operators/subchannel/SubchannelFourEqLinearOperator.h"
#include "operators/subchannel/SubchannelConstants.h"
#include "operators/subchannel/SubchannelHelpers.h"
#include "operators/OperatorBuilder.h"
#include "solvers/ColumnSolver.h"
#include "solvers/petsc/PetscKrylovSolverParameters.h"
#include "solvers/petsc/PetscKrylovSolver.h"
#include "solvers/petsc/PetscSNESSolverParameters.h"
#include "solvers/petsc/PetscSNESSolver.h"
#include "solvers/trilinos/TrilinosMLSolver.h"

#include "ampmesh/SiloIO.h"
#include "vectors/VectorBuilder.h"
#include "ampmesh/StructuredMeshHelper.h"
#include "discretization/structuredFaceDOFManager.h"
#include "discretization/simpleDOF_Manager.h"


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
    AMP::PIO::logAllNodes(log_file);
    AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

    // Read the input file
    AMP::shared_ptr<AMP::InputDatabase>  input_db ( new AMP::InputDatabase ( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile ( input_file , input_db );

//=============================================================================
// mesh and dof manager
//=============================================================================

    // Get the Mesh database and create the mesh parameters
    AMP_INSIST(input_db->keyExists("Mesh"), "Key ''Mesh'' is missing!");
    AMP::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase( "Mesh" );
    AMP::shared_ptr<AMP::Mesh::MeshParameters> meshParams(new AMP::Mesh::MeshParameters(mesh_db));
    meshParams->setComm(globalComm);

    // Create the meshes from the input database
    AMP::shared_ptr<AMP::Mesh::Mesh> subchannelMesh = AMP::Mesh::Mesh::buildMesh(meshParams);

    // get dof manager
    int DOFsPerFace[3]={1,1,3};
    AMP::Discretization::DOFManager::shared_ptr subchannelDOFManager = 
        AMP::Discretization::structuredFaceDOFManager::create( subchannelMesh, DOFsPerFace, 1 );

//=============================================================================
// physics model, parameters, and operator creation
//=============================================================================
    // get input and output variables
    AMP::LinearAlgebra::Variable::shared_ptr inputVariable  (new AMP::LinearAlgebra::Variable("flow"));
    AMP::LinearAlgebra::Variable::shared_ptr outputVariable (new AMP::LinearAlgebra::Variable("flow"));

    // create solution, rhs, and residual vectors
    AMP::LinearAlgebra::Vector::shared_ptr manufacturedVec = AMP::LinearAlgebra::createVector( subchannelDOFManager, inputVariable,  true );
    AMP::LinearAlgebra::Vector::shared_ptr solVec          = AMP::LinearAlgebra::createVector( subchannelDOFManager, inputVariable,  true );
    AMP::LinearAlgebra::Vector::shared_ptr rhsVec          = AMP::LinearAlgebra::createVector( subchannelDOFManager, outputVariable, true );
    AMP::LinearAlgebra::Vector::shared_ptr resVec          = AMP::LinearAlgebra::createVector( subchannelDOFManager, outputVariable, true );

    // get subchannel physics model
    AMP::shared_ptr<AMP::Database> subchannelPhysics_db = input_db->getDatabase("SubchannelPhysicsModel");
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModelParameters> params( new AMP::Operator::ElementPhysicsModelParameters(subchannelPhysics_db));
    AMP::shared_ptr<AMP::Operator::SubchannelPhysicsModel>  subchannelPhysicsModel (new AMP::Operator::SubchannelPhysicsModel(params));

    // Create the SubchannelOperatorParameters
    AMP::shared_ptr<AMP::Database> nonlinearOperator_db = input_db->getDatabase("SubchannelFourEqNonlinearOperator");
    AMP::shared_ptr<AMP::Database> linearOperator_db    = input_db->getDatabase("SubchannelFourEqLinearOperator");
    AMP::shared_ptr<AMP::Operator::SubchannelOperatorParameters> nonlinearOpParams(new AMP::Operator::SubchannelOperatorParameters( nonlinearOperator_db ));
    AMP::shared_ptr<AMP::Operator::SubchannelOperatorParameters> linearOpParams(   new AMP::Operator::SubchannelOperatorParameters( linearOperator_db ));
    nonlinearOpParams->d_Mesh = subchannelMesh ;
    nonlinearOpParams->d_subchannelPhysicsModel = subchannelPhysicsModel;
    nonlinearOpParams->clad_x = input_db->getDatabase("CladProperties")->getDoubleArray("x");
    nonlinearOpParams->clad_y = input_db->getDatabase("CladProperties")->getDoubleArray("y");
    nonlinearOpParams->clad_d = input_db->getDatabase("CladProperties")->getDoubleArray("d");
    linearOpParams->d_Mesh = subchannelMesh ;
    linearOpParams->d_subchannelPhysicsModel = subchannelPhysicsModel;
    linearOpParams->clad_x = input_db->getDatabase("CladProperties")->getDoubleArray("x");
    linearOpParams->clad_y = input_db->getDatabase("CladProperties")->getDoubleArray("y");
    linearOpParams->clad_d = input_db->getDatabase("CladProperties")->getDoubleArray("d");

    // create nonlinear operator
    AMP::shared_ptr<AMP::Operator::SubchannelFourEqNonlinearOperator> nonlinearOperator (new AMP::Operator::SubchannelFourEqNonlinearOperator(nonlinearOpParams));
    // reset the nonlinear operator
    nonlinearOperator->reset(nonlinearOpParams);

    // create linear operator
    AMP::shared_ptr<AMP::Operator::SubchannelFourEqLinearOperator> linearOperator (new AMP::Operator::SubchannelFourEqLinearOperator(linearOpParams));

    // pass creation test
    ut->passes(exeName+": creation");
    std::cout.flush();

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
    size_t N_subchannels = AMP::Operator::Subchannel::getNumberOfSubchannels( subchannelMesh );
    m_in = m_in/N_subchannels;

    // compute inlet enthalpy
    double Pin = Pout;
    double hin = 0.0;
    double rho_in = 1000;
    // iterate to find inlet pressure and inlet enthalpy
    for (int i=0; i<3; i++) {
       // compute inlet enthalpy using inlet temperature and outlet pressure
       std::map<std::string, AMP::shared_ptr<std::vector<double> > > enthalpyArgMap;
       enthalpyArgMap.insert(std::make_pair("temperature",AMP::shared_ptr<std::vector<double> >(new std::vector<double>(1,Tin))));
       enthalpyArgMap.insert(std::make_pair("pressure",   AMP::shared_ptr<std::vector<double> >(new std::vector<double>(1,Pin))));
       std::vector<double> enthalpyResult(1);
       subchannelPhysicsModel->getProperty("Enthalpy",enthalpyResult,enthalpyArgMap); 
       hin = enthalpyResult[0];
       // compute inlet density using computed inlet enthalpy and outlet pressure
       std::map<std::string, AMP::shared_ptr<std::vector<double> > > volumeArgMap_plus;
       volumeArgMap_plus.insert(std::make_pair("enthalpy",AMP::shared_ptr<std::vector<double> >(new std::vector<double>(1,hin))));
       volumeArgMap_plus.insert(std::make_pair("pressure",AMP::shared_ptr<std::vector<double> >(new std::vector<double>(1,Pin))));
       std::vector<double> volumeResult_plus(1);
       subchannelPhysicsModel->getProperty("SpecificVolume",volumeResult_plus,volumeArgMap_plus); 
       rho_in = 1.0/volumeResult_plus[0];
       // compute inlet pressure
       Pin = getSolutionPressure(input_db,H,Pout,rho_in,0);
    }
    std::cout<< "Inlet density:"<< rho_in <<std::endl;
    std::cout<< "Enthalpy Solution:"<< hin <<std::endl;

    // Compute the manufactured solution
    AMP::Mesh::Mesh::shared_ptr xyFaceMesh = 
        subchannelMesh->Subset( AMP::Mesh::StructuredMeshHelper::getXYFaceIterator( subchannelMesh , 0 ) );
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
    std::map<std::vector<double>,AMP::Mesh::MeshElement> interiorLateralFaceMap;
    std::map<std::vector<double>,AMP::Mesh::MeshElement> exteriorLateralFaceMap;
    nonlinearOperator->getLateralFaces(nonlinearOpParams->d_Mesh,interiorLateralFaceMap,exteriorLateralFaceMap);
    // loop over lateral faces
    for (face = face.begin(); face != face.end(); face++) {
       std::vector<double> faceCentroid = face->centroid();
       std::map<std::vector<double>,AMP::Mesh::MeshElement>::iterator lateralFaceIterator = interiorLateralFaceMap.find(faceCentroid);
       if (lateralFaceIterator != interiorLateralFaceMap.end()) {
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
       std::map<std::vector<double>,AMP::Mesh::MeshElement>::iterator lateralFaceIterator = interiorLateralFaceMap.find(faceCentroid);
       if (lateralFaceIterator != interiorLateralFaceMap.end()) {
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
    AMP::shared_ptr<AMP::Database> nonlinearSolver_db = input_db->getDatabase("NonlinearSolver"); 
  
    // get linear solver database
    AMP::shared_ptr<AMP::Database> linearSolver_db = nonlinearSolver_db->getDatabase("LinearSolver"); 
 
    // put manufactured RHS into resVec
    nonlinearOperator->reset(nonlinearOpParams);
    linearOperator->reset(nonlinearOperator->getJacobianParameters(solVec));
    linearOperator->residual(rhsVec, solVec, resVec);
   
    // create nonlinear solver parameters
    AMP::shared_ptr<AMP::Solver::PetscSNESSolverParameters> nonlinearSolverParams(new AMP::Solver::PetscSNESSolverParameters(nonlinearSolver_db));

    // change the next line to get the correct communicator out
    nonlinearSolverParams->d_comm = globalComm;
    nonlinearSolverParams->d_pOperator = nonlinearOperator;
    nonlinearSolverParams->d_pInitialGuess = solVec;

    // create nonlinear solver
    AMP::shared_ptr<AMP::Solver::PetscSNESSolver> nonlinearSolver(new AMP::Solver::PetscSNESSolver(nonlinearSolverParams));

    // create linear solver
    AMP::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver = nonlinearSolver->getKrylovSolver();

    // create preconditioner
    AMP::shared_ptr<AMP::Database> Preconditioner_db =  linearSolver_db->getDatabase("Preconditioner");
    AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> PreconditionerParams(new AMP::Solver::SolverStrategyParameters(Preconditioner_db));
    PreconditionerParams->d_pOperator = linearOperator;
    AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> linearFlowPreconditioner(new AMP::Solver::TrilinosMLSolver(PreconditionerParams));
    // set preconditioner
    linearSolver->setPreconditioner(linearFlowPreconditioner);

    // don't use zero initial guess
    nonlinearSolver->setZeroInitialGuess(false);

    // solve
/*
    std::cout<<" the initial matrix entries are: "<<std::endl;
    std::cout<< *linearOperator->getMatrix() <<std::endl;
    std::cout<<" those were the initial matrix entries. "<<std::endl;
    std::cout<<" the source is... "<<rhsVec->getGlobalSize() <<std::endl;
    std::cout<< *rhsVec <<std::endl;
    std::cout<<" ... was the source "<<std::endl;
*/
    nonlinearOperator->residual(rhsVec, solVec, resVec);
/*
    std::cout<<" the initial residual is... "<<rhsVec->getGlobalSize() <<std::endl;
    std::cout<< *resVec <<std::endl;
    std::cout<<" ... was the initial residual "<<std::endl;
    std::cout<<" the initial solution is... "<<rhsVec->getGlobalSize() <<std::endl;
    std::cout<< *solVec <<std::endl;
    std::cout<<" ... was the initial solution "<<std::endl;
*/
    nonlinearSolver->solve(rhsVec, solVec);
/*
    std::cout<<" the solution is... "<<rhsVec->getGlobalSize() <<std::endl;
    std::cout<< *solVec <<std::endl;
    std::cout<<" ... was the solution "<<std::endl;
*/
    nonlinearOperator->residual(rhsVec, solVec, resVec);
/*
    std::cout<<" the residual is... "<<rhsVec->getGlobalSize() <<std::endl;
    std::cout<< *resVec <<std::endl;
    std::cout<<" ... was the residual "<<std::endl;
    std::cout<<" the final matrix entries are: "<<std::endl;
    std::cout<< *linearOperator->getMatrix() <<std::endl;
    std::cout<<" those were the final matrix entries. "<<std::endl;
*/

//=============================================================================
// examine solution
//=============================================================================

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
        double h = h_scale*solVec->getValueByGlobalID(axialDofs[1]);
        double P = P_scale*solVec->getValueByGlobalID(axialDofs[2]);
        std::map<std::string, AMP::shared_ptr<std::vector<double> > > temperatureArgMap;
        temperatureArgMap.insert(std::make_pair("enthalpy",AMP::shared_ptr<std::vector<double> >(new std::vector<double>(1,h))));
        temperatureArgMap.insert(std::make_pair("pressure",AMP::shared_ptr<std::vector<double> >(new std::vector<double>(1,P))));
        std::vector<double> temperatureResult(1);
        subchannelPhysicsModel->getProperty("Temperature", temperatureResult, temperatureArgMap); 
        tempVec->setValueByGlobalID(tdofs[0],temperatureResult[0]);
        // Check that we recover the enthalpy from the temperature
        std::map<std::string, AMP::shared_ptr<std::vector<double> > > enthalpyArgMap;
        enthalpyArgMap.insert(std::make_pair("temperature",AMP::shared_ptr<std::vector<double> >(new std::vector<double>(1,temperatureResult[0]))));
        enthalpyArgMap.insert(std::make_pair("pressure",   AMP::shared_ptr<std::vector<double> >(new std::vector<double>(1,P))));
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
    double TinSol  = tempVec->getValueByGlobalID(tdofs[0]);
    std::cout<< "Inlet Computed Mass Flow Rate = " << m_scale*solVec->getValueByGlobalID(axialDofs[0]) << std::endl;
    std::cout<< "Inlet Computed Enthalpy = " << h_scale*solVec->getValueByGlobalID(axialDofs[1]) << std::endl;
    std::cout<< "Inlet Computed Pressure = " << P_scale*solVec->getValueByGlobalID(axialDofs[2]) << std::endl;
    std::cout<< "Inlet Computed Temperature = " << TinSol << std::endl;
    std::cout << std::endl;
    face = --((xyFaceMesh->getIterator(AMP::Mesh::Face,0)).end());
    subchannelDOFManager->getDOFs( face->globalID(), axialDofs );
    tempDOFManager->getDOFs( face->globalID(), tdofs );
    double ToutSol = tempVec->getValueByGlobalID(tdofs[0]);
    std::cout<< "Outlet Computed Mass Flow Rate = " << m_scale*solVec->getValueByGlobalID(axialDofs[0]) << std::endl;
    std::cout<< "Outlet Computed Enthalpy = " << h_scale*solVec->getValueByGlobalID(axialDofs[1]) << std::endl;
    std::cout<< "Outlet Computed Pressure = " << P_scale*solVec->getValueByGlobalID(axialDofs[2]) << std::endl;
    std::cout<< "Outlet Computed Temperature = " << ToutSol << std::endl;

    // Compute the error
    AMP::LinearAlgebra::Vector::shared_ptr absErrorVec = solVec->cloneVector();
    absErrorVec->axpy(-1.0,solVec,manufacturedVec);
    AMP::LinearAlgebra::Vector::shared_ptr relErrorVec = solVec->cloneVector();
    relErrorVec->divide(absErrorVec,manufacturedVec);
    for (size_t i=0; i<solVec->getLocalSize(); i++) {
        if ( manufacturedVec->getValueByLocalID(i) == 0 ) {
            double val = solVec->getValueByLocalID(i);
            relErrorVec->setValueByLocalID(i,fabs(val));
        }
    }
    double absErrorNorm = absErrorVec->L2Norm();
    double relErrorNorm = relErrorVec->L2Norm();

    // check that norm of relative error is less than tolerance
    double tol = input_db->getDoubleWithDefault("TOLERANCE",1e-6);
    if( relErrorNorm<=tol && fabs(Tin-TinSol)<tol ){
        ut->passes(exeName+": manufactured solution test");
    } else {
        ut->failure(exeName+": manufactured solution test");
    }

    // Print final solution
    AMP::Mesh::Mesh::shared_ptr channel0 = AMP::Operator::Subchannel::subsetForSubchannel( subchannelMesh, 0, 0 );
    face = AMP::Mesh::StructuredMeshHelper::getXYFaceIterator(channel0,0);
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
    std::cout<<"Delta T: " << ToutSol-TinSol << std::endl << std::endl;
    std::cout<<"L2 Norm of Absolute Error: "<<absErrorNorm<<std::endl;
    std::cout<<"L2 Norm of Relative Error: "<<relErrorNorm<<std::endl;

    input_db.reset();

/*
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

    std::vector<std::string> files;
    if ( argc >= 2 ) {
        files.resize(argc-1);
        for (int i=0; i<argc-1; i++)
            files[i] = std::string(argv[i+1]);
    } else {
        files.resize(2);
        files[0] = "testSubchannelFourEqMMS-1";
        files[1] = "testSubchannelFourEqMMS-2";
    }

    for (size_t i=0; i<files.size(); i++)
        flowTest(&ut,files[i]);

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}   

