#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include <iostream>
#include <string>

#include "boost/shared_ptr.hpp"

#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/PIO.h"
#include "materials/Material.h"

#include "ampmesh/MeshVariable.h"

#include "libmesh.h"

#include "operators/mechanics/VonMisesElastoPlasticModel.h"
#include "operators/mechanics/MechanicsNonlinearElement.h"
#include "operators/mechanics/MechanicsLinearElement.h"
#include "operators/mechanics/MechanicsLinearFEOperator.h"
#include "operators/mechanics/MechanicsNonlinearFEOperator.h"

#include "operators/boundary/DirichletMatrixCorrection.h"
#include "operators/boundary/DirichletVectorCorrection.h"

#include "operators/BVPOperatorParameters.h"
#include "operators/LinearBVPOperator.h"
#include "operators/NonlinearBVPOperator.h"


void myTest(AMP::UnitTest *ut)
{
  std::string exeName("testNonlinearMechanics-reset-1");
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

  AMP::PIO::logOnlyNodeZero(log_file);

  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  AMP_INSIST(input_db->keyExists("Mesh"), "Key ''Mesh'' is missing!");
  std::string mesh_file = input_db->getString("Mesh");

  AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter = AMP::Mesh::MeshManager::Adapter::shared_ptr ( new AMP::Mesh::MeshManager::Adapter () );
  meshAdapter->readExodusIIFile ( mesh_file.c_str() );

//  AMP_INSIST(input_db->keyExists("NumberOfLoadingSteps"), "Key ''NumberOfLoadingSteps'' is missing!");
//  int NumberOfLoadingSteps = input_db->getInteger("NumberOfLoadingSteps");

  //Material model shared by both the linear and nonlinear operators
  AMP_INSIST(input_db->keyExists("VonMises_Model"), "Key ''VonMises_Model'' is missing!");
  boost::shared_ptr<AMP::Database> matModel_db = input_db->getDatabase("VonMises_Model");
  boost::shared_ptr<AMP::Operator::MechanicsMaterialModelParameters> matModelParams(new
      AMP::Operator::MechanicsMaterialModelParameters( matModel_db ) );
  boost::shared_ptr<AMP::Operator::VonMisesElastoPlasticModel> matModel (new AMP::Operator::VonMisesElastoPlasticModel(matModelParams));

  for(int useReduced = 0; useReduced < 2; useReduced++) {

    std::string mechNonlinElemDbStr;
    std::string mechLinElemDbStr;
    if(useReduced) {
      AMP_INSIST(input_db->keyExists("Mechanics_Nonlinear_Element_Reduced"), "Key ''Mechanics_Nonlinear_Element_Reduced'' is missing!");
      AMP_INSIST(input_db->keyExists("Mechanics_Linear_Element_Reduced"), "Key ''Mechanics_Linear_Element_Reduced'' is missing!");
      mechNonlinElemDbStr = "Mechanics_Nonlinear_Element_Reduced";
      mechLinElemDbStr = "Mechanics_Linear_Element_Reduced";
    } else {
      AMP_INSIST(input_db->keyExists("Mechanics_Nonlinear_Element_Normal"), "Key ''Mechanics_Nonlinear_Element_Normal'' is missing!");
      AMP_INSIST(input_db->keyExists("Mechanics_Linear_Element_Normal"), "Key ''Mechanics_Linear_Element_Normal'' is missing!");
      mechNonlinElemDbStr = "Mechanics_Nonlinear_Element_Normal";
      mechLinElemDbStr = "Mechanics_Linear_Element_Normal";
    }

    boost::shared_ptr<AMP::Database> nonLinElemOp_db = input_db->getDatabase(mechNonlinElemDbStr);
    boost::shared_ptr<AMP::Operator::ElementOperationParameters> nonLinElemOpParams (new AMP::Operator::ElementOperationParameters( nonLinElemOp_db ));
    boost::shared_ptr<AMP::Operator::MechanicsNonlinearElement> mechNonlinElem (new AMP::Operator::MechanicsNonlinearElement( nonLinElemOpParams ));

    AMP_INSIST(input_db->keyExists("Mechanics_Nonlinear_Assembly"), "Key ''Mechanics_Nonlinear_Assembly'' is missing!");
    boost::shared_ptr<AMP::Database> mechNonlinAssembly_db = input_db->getDatabase("Mechanics_Nonlinear_Assembly");
    boost::shared_ptr<AMP::Operator::MechanicsNonlinearFEOperatorParameters> mechNonlinOpParams(new
        AMP::Operator::MechanicsNonlinearFEOperatorParameters( mechNonlinAssembly_db ));
    mechNonlinOpParams->d_materialModel = matModel;
    mechNonlinOpParams->d_elemOp = mechNonlinElem;
    mechNonlinOpParams->d_MeshAdapter = meshAdapter;
    boost::shared_ptr<AMP::Operator::MechanicsNonlinearFEOperator> mechNonlinOp (new AMP::Operator::MechanicsNonlinearFEOperator( mechNonlinOpParams ));
    mechNonlinOp->init();

    AMP::LinearAlgebra::Variable::shared_ptr displacementVariable = mechNonlinOp->getInputVariable(AMP::Operator::Mechanics::DISPLACEMENT);
    AMP::LinearAlgebra::Variable::shared_ptr residualVariable = mechNonlinOp->getOutputVariable();

    boost::shared_ptr<AMP::Database> linElemOp_db = input_db->getDatabase(mechLinElemDbStr);
    boost::shared_ptr<AMP::Operator::ElementOperationParameters> linElemOpParams (new AMP::Operator::ElementOperationParameters( linElemOp_db ));
    boost::shared_ptr<AMP::Operator::MechanicsLinearElement> mechLinElem (new AMP::Operator::MechanicsLinearElement( linElemOpParams ));

    AMP_INSIST(input_db->keyExists("Mechanics_Linear_Assembly"), "Key ''Mechanics_Linear_Assembly'' is missing!");
    boost::shared_ptr<AMP::Database> mechLinAssembly_db = input_db->getDatabase("Mechanics_Linear_Assembly");
    boost::shared_ptr<AMP::Operator::MechanicsLinearFEOperatorParameters> mechLinOpParams(new
        AMP::Operator::MechanicsLinearFEOperatorParameters( mechLinAssembly_db ));
    mechLinOpParams->d_materialModel = matModel;
    mechLinOpParams->d_elemOp = mechLinElem;
    mechLinOpParams->d_MeshAdapter = meshAdapter;
    boost::shared_ptr<AMP::Operator::MechanicsLinearFEOperator> mechLinOp (new AMP::Operator::MechanicsLinearFEOperator( mechLinOpParams ));

    AMP_INSIST(input_db->keyExists("Displacement_Boundary"), "Key ''Displacement_Boundary'' is missing!");
    boost::shared_ptr<AMP::Database> disp_db = input_db->getDatabase("Displacement_Boundary");
    boost::shared_ptr<AMP::Operator::DirichletMatrixCorrectionParameters> dirichletOpParams (new 
        AMP::Operator::DirichletMatrixCorrectionParameters( disp_db ) );
    dirichletOpParams->d_inputMatrix = mechLinOp->getMatrix();
    //This is just the variable used to extract the dof_map.
    //This boundary operator itself has an empty input and output variable
    dirichletOpParams->d_variable = residualVariable;
    dirichletOpParams->d_MeshAdapter = meshAdapter;
    boost::shared_ptr<AMP::Operator::DirichletMatrixCorrection> dirichletMatOp (new 
        AMP::Operator::DirichletMatrixCorrection( dirichletOpParams ) );

    boost::shared_ptr<AMP::Operator::DirichletVectorCorrectionParameters> dirichletDispInVecParams(new 
        AMP::Operator::DirichletVectorCorrectionParameters(disp_db));
    //This has an in-place apply. So, it has an empty input variable and
    //the output variable is the same as what it is operating on. 
    dirichletDispInVecParams->d_variable = displacementVariable;
    dirichletDispInVecParams->d_MeshAdapter = meshAdapter;
    boost::shared_ptr<AMP::Operator::DirichletVectorCorrection> dirichletDispInVecOp(new 
        AMP::Operator::DirichletVectorCorrection(dirichletDispInVecParams));

    boost::shared_ptr<AMP::Operator::DirichletVectorCorrectionParameters> dirichletDispOutVecParams(new 
        AMP::Operator::DirichletVectorCorrectionParameters(disp_db));
    //This has an in-place apply. So, it has an empty input variable and
    //the output variable is the same as what it is operating on. 
    dirichletDispOutVecParams->d_variable = residualVariable;
    dirichletDispOutVecParams->d_MeshAdapter = meshAdapter;
    boost::shared_ptr<AMP::Operator::DirichletVectorCorrection> dirichletDispOutVecOp (new 
        AMP::Operator::DirichletVectorCorrection(dirichletDispOutVecParams));

    AMP_INSIST(input_db->keyExists("LinearBVPOperator"), "Key ''LinearBVPOperator'' is missing!");
    boost::shared_ptr<AMP::Database> linBvpOp_db = input_db->getDatabase("LinearBVPOperator");
    boost::shared_ptr<AMP::Operator::BVPOperatorParameters> linBvpOperatorParams (new AMP::Operator::BVPOperatorParameters(linBvpOp_db));
    linBvpOperatorParams->d_volumeOperator = mechLinOp;
    linBvpOperatorParams->d_boundaryOperator = dirichletMatOp;
    boost::shared_ptr<AMP::Operator::LinearBVPOperator> linBvpOperator (new AMP::Operator::LinearBVPOperator(linBvpOperatorParams));

    AMP_INSIST(input_db->keyExists("NonlinearBVPOperator"), "Key ''NonlinearBVPOperator'' is missing!");
    boost::shared_ptr<AMP::Database> nonlinBvpOp_db = input_db->getDatabase("NonlinearBVPOperator");
    boost::shared_ptr<AMP::Operator::BVPOperatorParameters> nonlinBvpOperatorParams (new AMP::Operator::BVPOperatorParameters(nonlinBvpOp_db));
    nonlinBvpOperatorParams->d_volumeOperator = mechNonlinOp;
    nonlinBvpOperatorParams->d_boundaryOperator = dirichletDispOutVecOp;
    boost::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinBvpOperator (new AMP::Operator::NonlinearBVPOperator(nonlinBvpOperatorParams));

    AMP_INSIST(input_db->keyExists("Load_Boundary"), "Key ''Load_Boundary'' is missing!");
    boost::shared_ptr<AMP::Database> load_db = input_db->getDatabase("Load_Boundary");
    boost::shared_ptr<AMP::Operator::DirichletVectorCorrectionParameters> dirichletLoadVecParams (new 
        AMP::Operator::DirichletVectorCorrectionParameters(load_db));
    //This has an in-place apply. So, it has an empty input variable and
    //the output variable is the same as what it is operating on. 
    dirichletLoadVecParams->d_variable = residualVariable;
    dirichletLoadVecParams->d_MeshAdapter = meshAdapter;
    boost::shared_ptr<AMP::Operator::DirichletVectorCorrection> dirichletLoadVecOp (new 
        AMP::Operator::DirichletVectorCorrection(dirichletLoadVecParams));

    AMP::LinearAlgebra::Vector::shared_ptr nullVec;

    AMP::LinearAlgebra::Vector::shared_ptr mechNlSolVec = meshAdapter->createVector( displacementVariable );
    AMP::LinearAlgebra::Vector::shared_ptr mechNlRhsVec = meshAdapter->createVector( residualVariable );
    AMP::LinearAlgebra::Vector::shared_ptr mechNlResVec = meshAdapter->createVector( residualVariable );

    mechNlRhsVec->setToScalar(0.0);
    dirichletLoadVecOp->apply(nullVec, nullVec, mechNlRhsVec, 1.0, 0.0);

    for(int i = 0; i < 3; i++) {
      //Initial guess for NL solver must satisfy the displacement boundary
      //conditions
      mechNlSolVec->setRandomValues();
      dirichletDispInVecOp->apply(nullVec, nullVec, mechNlSolVec, 1.0, 0.0);

      nonlinBvpOperator->apply(mechNlRhsVec, mechNlSolVec, mechNlResVec, 1.0, -1.0);
    }//end for i

    mechNonlinOp->reset(mechNonlinOpParams);

    ut->passes(exeName + " : " + mechNonlinElemDbStr);

  }//end for useReduced

}

int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

    try {
        myTest(&ut);
    } catch (std::exception &err) {
        std::cout << "ERROR: While testing "<<argv[0] << err.what() << std::endl;
        ut.failure("ERROR: While testing");
    } catch( ... ) {
        std::cout << "ERROR: While testing "<<argv[0] << "An unknown exception was thrown." << std::endl;
        ut.failure("ERROR: While testing");
    }
   
    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}   


