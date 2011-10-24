#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include <string>
#include "ampmesh/MeshManager.h"
#include "materials/Material.h"
#include "boost/shared_ptr.hpp"
#include "utils/InputDatabase.h"
#include "utils/Utilities.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"
#include "utils/Database.h"
#include "operators/NeutronicsRhs.h"
#include "ampmesh/MeshAdapter.h"
#include "vectors/Variable.h"
#include "operators/VolumeIntegralOperator.h"

#include "vectors/Vector.h"
#include "vectors/Variable.h"
#include "ampmesh/SiloIO.h"
#include "matrices/Matrix.h"
#include "matrices/trilinos/EpetraMatrix.h"
//#include <math>

#include "ampmesh/SiloIO.h"
#include "vectors/Vector.h"

#include "flow/ConsMomentumGalWFLinearElement.h"
#include "flow/ConsMomentumGalWFLinearFEOperator.h"
#include "flow/ConsMassGalWFLinearElement.h"
#include "flow/ConsMassGalWFLinearFEOperator.h"

#include "operators/EpetraMatrixOperator.h"
#include "operators/EpetraMatrixOperatorParameters.h"
#include "operators/LinearBVPOperator.h"
#include "operators/OperatorBuilder.h"

#include "operators/BlockOperator.h"
#include "operators/TrilinosMatrixShellOperator.h"
#include "operators/PetscMatrixShellOperator.h"

#include "solvers/TrilinosMLSolver.h"
#include "solvers/ColumnSolver.h"
#include "solvers/PetscKrylovSolverParameters.h"
#include "solvers/PetscKrylovSolver.h"
#include "solvers/PetscSNESSolverParameters.h"
#include "solvers/PetscSNESSolver.h"


#define ITFAILS ut.failure(__LINE__);
#define UNIT_TEST(a) if (!(a)) ut.failure(__LINE__);

double h_Size = 1.46286; 

void myTest(AMP::UnitTest *ut, std::string exeName)

{

	std::string input_file = "input_" + exeName;
        std::string silo_name = exeName;

        boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
	AMP::AMP_MPI globalComm = AMP::AMP_MPI(AMP_COMM_WORLD);
        AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
        input_db->printClassData(AMP::plog);

	AMP_INSIST(input_db->keyExists("NumberOfMeshes"), "Key does not exist");
	// int numMeshes = input_db->getInteger("NumberOfMeshes");

	AMP::Mesh::MeshManagerParameters::shared_ptr meshmgrParams (new AMP::Mesh::MeshManagerParameters ( input_db ) );
	AMP::Mesh::MeshManager::shared_ptr manager (new AMP::Mesh::MeshManager ( meshmgrParams ) );
	AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter = manager->getMesh ("cube");
	
  /////////////////////////////////////////////////
  //   CREATE THE Conservation of Mass Operator  //
  /////////////////////////////////////////////////

	boost::shared_ptr<AMP::Operator::ElementPhysicsModel> FlowTransportModel; 
	AMP_INSIST ( input_db->keyExists("ConsMassLinearFEOperator"),"key missing!");
	boost::shared_ptr<AMP::Operator::LinearBVPOperator> ConsMassOperator = boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter, "ConsMassLinearBVPOperator", input_db, FlowTransportModel));

        AMP::pout << "Finished creating Mass Operator" << std::endl;
  /////////////////////////////////////////////////
  //   CREATE THE Conservation of Momentum Operator  //
  /////////////////////////////////////////////////

	AMP_INSIST ( input_db->keyExists("ConsMomentumLinearFEOperator"),"key missing!");
	boost::shared_ptr<AMP::Operator::LinearBVPOperator> ConsMomentumOperator = boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter, "ConsMomentumLinearBVPOperator", input_db, FlowTransportModel));
	
        //ease of use shared pointer
        AMP::LinearAlgebra::Vector::shared_ptr nullVec;

        AMP::LinearAlgebra::Matrix::shared_ptr FMat = ConsMomentumOperator->getMatrix();
        AMP::LinearAlgebra::Matrix::shared_ptr BMat = ConsMassOperator->getMatrix();

        AMP::LinearAlgebra::Matrix::shared_ptr BtMat =  BMat->transpose();
        AMP::LinearAlgebra::Matrix::shared_ptr zeroMat = meshAdapter->createMatrix ( ConsMassOperator->getOutputVariable(), ConsMassOperator->getOutputVariable() ) ; 
        zeroMat->zero();

        AMP::LinearAlgebra::Variable::shared_ptr velocityVar   = ConsMomentumOperator->getOutputVariable();
        AMP::LinearAlgebra::Variable::shared_ptr pressureVar   = ConsMassOperator->getOutputVariable();

        boost::shared_ptr<AMP::LinearAlgebra::MultiVariable> globalMultiVar(new AMP::LinearAlgebra::MultiVariable("inputVariable"));
        globalMultiVar->add(velocityVar);
        globalMultiVar->add(pressureVar);

        AMP::LinearAlgebra::Vector::shared_ptr globalSolVec = manager->createVector ( globalMultiVar );
        AMP::LinearAlgebra::Vector::shared_ptr globalRhsVec = manager->createVector ( globalMultiVar );
//        AMP::LinearAlgebra::Vector::shared_ptr globalResVec = manager->createVector ( globalMultiVar );

	boost::shared_ptr<AMP::Database> dummy_db;
        boost::shared_ptr<AMP::Operator::EpetraMatrixOperatorParameters> dummyParams1(new AMP::Operator::EpetraMatrixOperatorParameters( dummy_db ));
        dummyParams1->d_Matrix = &(boost::dynamic_pointer_cast<AMP::LinearAlgebra::EpetraMatrix>(BtMat)->getEpetra_CrsMatrix());
        boost::shared_ptr<AMP::Operator::EpetraMatrixOperator> bTOperator ( new AMP::Operator::EpetraMatrixOperator (dummyParams1) );
        bTOperator->setVariables(ConsMassOperator->getOutputVariable(), ConsMassOperator->getInputVariable());

        boost::shared_ptr<AMP::Operator::EpetraMatrixOperatorParameters> dummyParams2(new AMP::Operator::EpetraMatrixOperatorParameters( dummy_db ));
        dummyParams2->d_Matrix = &(boost::dynamic_pointer_cast<AMP::LinearAlgebra::EpetraMatrix>(zeroMat)->getEpetra_CrsMatrix());
         boost::shared_ptr<AMP::Operator::EpetraMatrixOperator> zeroOperator ( new AMP::Operator::EpetraMatrixOperator (dummyParams2) );
        zeroOperator->setVariables(ConsMassOperator->getOutputVariable(), ConsMassOperator->getOutputVariable());

        AMP_INSIST(input_db->keyExists("LinearSolver"),   "Key ''LinearSolver'' is missing!");
        boost::shared_ptr<AMP::Database>  mlSolver_db  = input_db->getDatabase("LinearSolver"); 
       
        boost::shared_ptr<AMP::Solver::SolverStrategyParameters> convdiffSolverParams(new AMP::Solver::SolverStrategyParameters(mlSolver_db));
        convdiffSolverParams->d_pOperator = ConsMomentumOperator;
        boost::shared_ptr<AMP::Solver::TrilinosMLSolver>  convdiffSolver(new AMP::Solver::TrilinosMLSolver(convdiffSolverParams));
        convdiffSolver->setZeroInitialGuess(false);

        AMP::LinearAlgebra::Vector::shared_ptr diagonalVec = FMat->extractDiagonal();
        AMP::LinearAlgebra::Vector::shared_ptr diagonalInvVec = diagonalVec->cloneVector();
        diagonalInvVec->reciprocal(diagonalVec);

        AMP::LinearAlgebra::Matrix::shared_ptr DMat = FMat->cloneMatrix() ;
        DMat->zero();
        DMat->setDiagonal(diagonalVec);

        AMP::LinearAlgebra::Matrix::shared_ptr DInvMat = FMat->cloneMatrix() ;
        DInvMat->zero();
        DInvMat->setDiagonal(diagonalInvVec);

        AMP::LinearAlgebra::Matrix::shared_ptr schurMat = zeroMat->cloneMatrix() ;

        AMP::LinearAlgebra::Matrix::shared_ptr DInvBtMat = AMP::LinearAlgebra::Matrix::matMultiply(DInvMat, BtMat);

        schurMat = AMP::LinearAlgebra::Matrix::matMultiply(BtMat, DInvBtMat);

        boost::shared_ptr<AMP::Operator::EpetraMatrixOperatorParameters> dummyParams3(new AMP::Operator::EpetraMatrixOperatorParameters( dummy_db ));
        dummyParams3->d_Matrix = &(boost::dynamic_pointer_cast<AMP::LinearAlgebra::EpetraMatrix>(schurMat)->getEpetra_CrsMatrix());
        boost::shared_ptr<AMP::Operator::EpetraMatrixOperator> schurMatOperator ( new AMP::Operator::EpetraMatrixOperator (dummyParams3) );
        schurMatOperator->setVariables(ConsMassOperator->getOutputVariable(), ConsMassOperator->getOutputVariable());

        boost::shared_ptr<AMP::Solver::SolverStrategyParameters> schurMatSolverParams(new AMP::Solver::SolverStrategyParameters(mlSolver_db));
        schurMatSolverParams->d_pOperator = schurMatOperator;
        boost::shared_ptr<AMP::Solver::TrilinosMLSolver>  schurMatSolver(new AMP::Solver::TrilinosMLSolver(schurMatSolverParams));
        schurMatSolver->setZeroInitialGuess(false);

        AMP::LinearAlgebra::Vector::shared_ptr velocityRhsVec = globalRhsVec->subsetVectorForVariable(velocityVar); 
        AMP::LinearAlgebra::Vector::shared_ptr pressureRhsVec = globalRhsVec->subsetVectorForVariable(pressureVar);

        AMP::LinearAlgebra::Vector::shared_ptr velocitySolVec = globalSolVec->subsetVectorForVariable(velocityVar); 
        AMP::LinearAlgebra::Vector::shared_ptr pressureSolVec = globalSolVec->subsetVectorForVariable(pressureVar);
       
        AMP::LinearAlgebra::Vector::shared_ptr pressureUpdateVec = pressureSolVec->cloneVector();
        AMP::LinearAlgebra::Vector::shared_ptr velocityUpdateVec = velocitySolVec->cloneVector();
        
//        AMP::LinearAlgebra::Vector::shared_ptr pressurePrimeVec = pressureSolVec->cloneVector();
        AMP::LinearAlgebra::Vector::shared_ptr velocityPrimeVec = velocitySolVec->cloneVector();

        //SIMPLE(Semi Implicit Method for Pressure Linked Equations) ALGORITHM

        // STEP 1 :
        BtMat->mult( pressureSolVec , velocityRhsVec );
        convdiffSolver->solve( velocityRhsVec , velocityPrimeVec );

        // STEP 2 :
        BMat->mult( velocityPrimeVec , pressureRhsVec );
        schurMatSolver->solve( pressureRhsVec, pressureUpdateVec );

        // STEP 3 :
        DInvBtMat->mult(pressureUpdateVec, velocityUpdateVec); 
        velocitySolVec->subtract( velocityPrimeVec, velocityUpdateVec);

        // STEP 4 :
        pressureSolVec->add(velocityPrimeVec, pressureUpdateVec );

}

int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  AMP::UnitTest ut;
  std::vector<std::string> exeNames;
  exeNames.push_back("testLinearFlow-1");

  for(unsigned int i=0; i < exeNames.size(); i++)
  {
    try
    {
      myTest(&ut, exeNames[i]);
    } catch (std::exception &err) {
      std::cout << "ERROR: While testing " << argv[0] << err.what() << std::endl;
      ut.failure("ERROR: While testing");
    } catch( ... ){
      std::cout << "ERROR: While testing " << argv[0] << "An unknown exception was thrown." << std::endl;
      ut.failure("ERROR: While testing");
    }

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;

  }

}




























