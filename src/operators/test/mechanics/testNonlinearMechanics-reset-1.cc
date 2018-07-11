
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/materials/Material.h"
#include "AMP/operators/BVPOperatorParameters.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/NonlinearBVPOperator.h"
#include "AMP/operators/boundary/DirichletMatrixCorrection.h"
#include "AMP/operators/boundary/DirichletVectorCorrection.h"
#include "AMP/operators/mechanics/MechanicsLinearElement.h"
#include "AMP/operators/mechanics/MechanicsLinearFEOperator.h"
#include "AMP/operators/mechanics/MechanicsNonlinearElement.h"
#include "AMP/operators/mechanics/MechanicsNonlinearFEOperator.h"
#include "AMP/operators/mechanics/VonMisesElastoPlasticModel.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/InputDatabase.h"
#include "AMP/utils/InputManager.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/shared_ptr.h"
#include "AMP/vectors/VectorBuilder.h"

#include <iostream>
#include <string>


void myTest( AMP::UnitTest *ut )
{
    std::string exeName( "testNonlinearMechanics-reset-1" );
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::PIO::logOnlyNodeZero( log_file );

    AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );
    input_db->printClassData( AMP::plog );

    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    AMP::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase( "Mesh" );
    AMP::shared_ptr<AMP::Mesh::MeshParameters> meshParams(
        new AMP::Mesh::MeshParameters( mesh_db ) );
    meshParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    AMP::Mesh::Mesh::shared_ptr meshAdapter = AMP::Mesh::Mesh::buildMesh( meshParams );

    AMP::Discretization::DOFManager::shared_ptr dofMap =
        AMP::Discretization::simpleDOFManager::create(
            meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 3, true );

    //  AMP_INSIST(input_db->keyExists("NumberOfLoadingSteps"), "Key ''NumberOfLoadingSteps'' is
    //  missing!");
    //  int NumberOfLoadingSteps = input_db->getInteger("NumberOfLoadingSteps");

    // Material model shared by both the linear and nonlinear operators
    AMP_INSIST( input_db->keyExists( "VonMises_Model" ), "Key ''VonMises_Model'' is missing!" );
    AMP::shared_ptr<AMP::Database> matModel_db = input_db->getDatabase( "VonMises_Model" );
    AMP::shared_ptr<AMP::Operator::MechanicsMaterialModelParameters> matModelParams(
        new AMP::Operator::MechanicsMaterialModelParameters( matModel_db ) );
    AMP::shared_ptr<AMP::Operator::VonMisesElastoPlasticModel> matModel(
        new AMP::Operator::VonMisesElastoPlasticModel( matModelParams ) );

    for ( int useReduced = 0; useReduced < 2; useReduced++ ) {

        std::string mechNonlinElemDbStr;
        std::string mechLinElemDbStr;
        if ( useReduced ) {
            AMP_INSIST( input_db->keyExists( "Mechanics_Nonlinear_Element_Reduced" ),
                        "Key ''Mechanics_Nonlinear_Element_Reduced'' is missing!" );
            AMP_INSIST( input_db->keyExists( "Mechanics_Linear_Element_Reduced" ),
                        "Key ''Mechanics_Linear_Element_Reduced'' is missing!" );
            mechNonlinElemDbStr = "Mechanics_Nonlinear_Element_Reduced";
            mechLinElemDbStr    = "Mechanics_Linear_Element_Reduced";
        } else {
            AMP_INSIST( input_db->keyExists( "Mechanics_Nonlinear_Element_Normal" ),
                        "Key ''Mechanics_Nonlinear_Element_Normal'' is missing!" );
            AMP_INSIST( input_db->keyExists( "Mechanics_Linear_Element_Normal" ),
                        "Key ''Mechanics_Linear_Element_Normal'' is missing!" );
            mechNonlinElemDbStr = "Mechanics_Nonlinear_Element_Normal";
            mechLinElemDbStr    = "Mechanics_Linear_Element_Normal";
        }

        AMP::shared_ptr<AMP::Database> nonLinElemOp_db =
            input_db->getDatabase( mechNonlinElemDbStr );
        AMP::shared_ptr<AMP::Operator::ElementOperationParameters> nonLinElemOpParams(
            new AMP::Operator::ElementOperationParameters( nonLinElemOp_db ) );
        AMP::shared_ptr<AMP::Operator::MechanicsNonlinearElement> mechNonlinElem(
            new AMP::Operator::MechanicsNonlinearElement( nonLinElemOpParams ) );

        AMP_INSIST( input_db->keyExists( "Mechanics_Nonlinear_Assembly" ),
                    "Key ''Mechanics_Nonlinear_Assembly'' is missing!" );
        AMP::shared_ptr<AMP::Database> mechNonlinAssembly_db =
            input_db->getDatabase( "Mechanics_Nonlinear_Assembly" );
        AMP::shared_ptr<AMP::Operator::MechanicsNonlinearFEOperatorParameters> mechNonlinOpParams(
            new AMP::Operator::MechanicsNonlinearFEOperatorParameters( mechNonlinAssembly_db ) );
        mechNonlinOpParams->d_materialModel                                  = matModel;
        mechNonlinOpParams->d_elemOp                                         = mechNonlinElem;
        mechNonlinOpParams->d_Mesh                                           = meshAdapter;
        mechNonlinOpParams->d_dofMap[AMP::Operator::Mechanics::DISPLACEMENT] = dofMap;
        AMP::shared_ptr<AMP::Operator::MechanicsNonlinearFEOperator> mechNonlinOp(
            new AMP::Operator::MechanicsNonlinearFEOperator( mechNonlinOpParams ) );

        AMP::LinearAlgebra::Variable::shared_ptr var = mechNonlinOp->getOutputVariable();

        AMP::shared_ptr<AMP::Database> linElemOp_db = input_db->getDatabase( mechLinElemDbStr );
        AMP::shared_ptr<AMP::Operator::ElementOperationParameters> linElemOpParams(
            new AMP::Operator::ElementOperationParameters( linElemOp_db ) );
        AMP::shared_ptr<AMP::Operator::MechanicsLinearElement> mechLinElem(
            new AMP::Operator::MechanicsLinearElement( linElemOpParams ) );

        AMP_INSIST( input_db->keyExists( "Mechanics_Linear_Assembly" ),
                    "Key ''Mechanics_Linear_Assembly'' is missing!" );
        AMP::shared_ptr<AMP::Database> mechLinAssembly_db =
            input_db->getDatabase( "Mechanics_Linear_Assembly" );
        AMP::shared_ptr<AMP::Operator::MechanicsLinearFEOperatorParameters> mechLinOpParams(
            new AMP::Operator::MechanicsLinearFEOperatorParameters( mechLinAssembly_db ) );
        mechLinOpParams->d_materialModel = matModel;
        mechLinOpParams->d_elemOp        = mechLinElem;
        mechLinOpParams->d_Mesh          = meshAdapter;
        mechLinOpParams->d_inDofMap      = dofMap;
        mechLinOpParams->d_outDofMap     = dofMap;
        AMP::shared_ptr<AMP::Operator::MechanicsLinearFEOperator> mechLinOp(
            new AMP::Operator::MechanicsLinearFEOperator( mechLinOpParams ) );

        AMP_INSIST( input_db->keyExists( "Displacement_Boundary" ),
                    "Key ''Displacement_Boundary'' is missing!" );
        AMP::shared_ptr<AMP::Database> disp_db = input_db->getDatabase( "Displacement_Boundary" );
        AMP::shared_ptr<AMP::Operator::DirichletMatrixCorrectionParameters> dirichletOpParams(
            new AMP::Operator::DirichletMatrixCorrectionParameters( disp_db ) );
        dirichletOpParams->d_inputMatrix = mechLinOp->getMatrix();
        // This is just the variable used to extract the dof_map.
        // This boundary operator itself has an empty input and output variable
        dirichletOpParams->d_variable = var;
        dirichletOpParams->d_Mesh     = meshAdapter;
        AMP::shared_ptr<AMP::Operator::DirichletMatrixCorrection> dirichletMatOp(
            new AMP::Operator::DirichletMatrixCorrection( dirichletOpParams ) );

        AMP::shared_ptr<AMP::Operator::DirichletVectorCorrectionParameters>
            dirichletDispInVecParams(
                new AMP::Operator::DirichletVectorCorrectionParameters( disp_db ) );
        // This has an in-place apply. So, it has an empty input variable and
        // the output variable is the same as what it is operating on.
        dirichletDispInVecParams->d_variable = var;
        dirichletDispInVecParams->d_Mesh     = meshAdapter;
        AMP::shared_ptr<AMP::Operator::DirichletVectorCorrection> dirichletDispInVecOp(
            new AMP::Operator::DirichletVectorCorrection( dirichletDispInVecParams ) );

        AMP::shared_ptr<AMP::Operator::DirichletVectorCorrectionParameters>
            dirichletDispOutVecParams(
                new AMP::Operator::DirichletVectorCorrectionParameters( disp_db ) );
        // This has an in-place apply. So, it has an empty input variable and
        // the output variable is the same as what it is operating on.
        dirichletDispOutVecParams->d_variable = var;
        dirichletDispOutVecParams->d_Mesh     = meshAdapter;
        AMP::shared_ptr<AMP::Operator::DirichletVectorCorrection> dirichletDispOutVecOp(
            new AMP::Operator::DirichletVectorCorrection( dirichletDispOutVecParams ) );

        AMP_INSIST( input_db->keyExists( "LinearBVPOperator" ),
                    "Key ''LinearBVPOperator'' is missing!" );
        AMP::shared_ptr<AMP::Database> linBvpOp_db = input_db->getDatabase( "LinearBVPOperator" );
        AMP::shared_ptr<AMP::Operator::BVPOperatorParameters> linBvpOperatorParams(
            new AMP::Operator::BVPOperatorParameters( linBvpOp_db ) );
        linBvpOperatorParams->d_volumeOperator   = mechLinOp;
        linBvpOperatorParams->d_boundaryOperator = dirichletMatOp;
        AMP::shared_ptr<AMP::Operator::LinearBVPOperator> linBvpOperator(
            new AMP::Operator::LinearBVPOperator( linBvpOperatorParams ) );

        AMP_INSIST( input_db->keyExists( "NonlinearBVPOperator" ),
                    "Key ''NonlinearBVPOperator'' is missing!" );
        AMP::shared_ptr<AMP::Database> nonlinBvpOp_db =
            input_db->getDatabase( "NonlinearBVPOperator" );
        AMP::shared_ptr<AMP::Operator::BVPOperatorParameters> nonlinBvpOperatorParams(
            new AMP::Operator::BVPOperatorParameters( nonlinBvpOp_db ) );
        nonlinBvpOperatorParams->d_volumeOperator   = mechNonlinOp;
        nonlinBvpOperatorParams->d_boundaryOperator = dirichletDispOutVecOp;
        AMP::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinBvpOperator(
            new AMP::Operator::NonlinearBVPOperator( nonlinBvpOperatorParams ) );

        AMP_INSIST( input_db->keyExists( "Load_Boundary" ), "Key ''Load_Boundary'' is missing!" );
        AMP::shared_ptr<AMP::Database> load_db = input_db->getDatabase( "Load_Boundary" );
        AMP::shared_ptr<AMP::Operator::DirichletVectorCorrectionParameters> dirichletLoadVecParams(
            new AMP::Operator::DirichletVectorCorrectionParameters( load_db ) );
        // This has an in-place apply. So, it has an empty input variable and
        // the output variable is the same as what it is operating on.
        dirichletLoadVecParams->d_variable = var;
        dirichletLoadVecParams->d_Mesh     = meshAdapter;
        AMP::shared_ptr<AMP::Operator::DirichletVectorCorrection> dirichletLoadVecOp(
            new AMP::Operator::DirichletVectorCorrection( dirichletLoadVecParams ) );

        AMP::LinearAlgebra::Vector::shared_ptr nullVec;

        AMP::LinearAlgebra::Vector::shared_ptr mechNlSolVec =
            AMP::LinearAlgebra::createVector( dofMap, var, true );
        AMP::LinearAlgebra::Vector::shared_ptr mechNlRhsVec = mechNlSolVec->cloneVector();
        AMP::LinearAlgebra::Vector::shared_ptr mechNlResVec = mechNlSolVec->cloneVector();

        mechNlRhsVec->setToScalar( 0.0 );
        dirichletLoadVecOp->apply( nullVec, mechNlRhsVec );

        for ( int i = 0; i < 3; i++ ) {
            // Initial guess for NL solver must satisfy the displacement boundary
            // conditions
            mechNlSolVec->setRandomValues();
            dirichletDispInVecOp->apply( nullVec, mechNlSolVec );
            mechNlSolVec->makeConsistent( AMP::LinearAlgebra::Vector::ScatterType::CONSISTENT_SET );
            nonlinBvpOperator->residual( mechNlRhsVec, mechNlSolVec, mechNlResVec );
        } // end for i

        mechNonlinOp->reset( mechNonlinOpParams );

        ut->passes( exeName + " : " + mechNonlinElemDbStr );

    } // end for useReduced
}

int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    myTest( &ut );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
