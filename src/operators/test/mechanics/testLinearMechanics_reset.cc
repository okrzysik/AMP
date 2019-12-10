
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/operators/boundary/DirichletMatrixCorrection.h"
#include "AMP/operators/boundary/DirichletVectorCorrection.h"
#include "AMP/operators/mechanics/IsotropicElasticModel.h"
#include "AMP/operators/mechanics/MechanicsLinearElement.h"
#include "AMP/operators/mechanics/MechanicsLinearFEOperator.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/shared_ptr.h"
#include "AMP/vectors/VectorBuilder.h"

#include <iostream>
#include <string>


static void myTest( AMP::UnitTest *ut )
{
    std::string exeName( "testLinearMechanics-reset-1" );
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::PIO::logOnlyNodeZero( log_file );


    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    AMP::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase( "Mesh" );
    AMP::shared_ptr<AMP::Mesh::MeshParameters> meshParams(
        new AMP::Mesh::MeshParameters( mesh_db ) );
    meshParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    AMP::Mesh::Mesh::shared_ptr meshAdapter = AMP::Mesh::Mesh::buildMesh( meshParams );

    AMP_INSIST( input_db->keyExists( "Isotropic_Model" ), "Key ''Isotropic_Model'' is missing!" );
    AMP::shared_ptr<AMP::Database> matModel_db = input_db->getDatabase( "Isotropic_Model" );
    AMP::shared_ptr<AMP::Operator::MechanicsMaterialModelParameters> matModelParams(
        new AMP::Operator::MechanicsMaterialModelParameters( matModel_db ) );
    AMP::shared_ptr<AMP::Operator::IsotropicElasticModel> isotropicModel(
        new AMP::Operator::IsotropicElasticModel( matModelParams ) );

    for ( int useReduced = 0; useReduced < 2; useReduced++ ) {

        std::string mechElemDbStr;
        if ( useReduced ) {
            AMP_INSIST( input_db->keyExists( "Mechanics_Linear_Element_Reduced" ),
                        "Key ''Mechanics_Linear_Element_Reduced'' is missing!" );
            mechElemDbStr = "Mechanics_Linear_Element_Reduced";
        } else {
            AMP_INSIST( input_db->keyExists( "Mechanics_Linear_Element_Normal" ),
                        "Key ''Mechanics_Linear_Element_Normal'' is missing!" );
            mechElemDbStr = "Mechanics_Linear_Element_Normal";
        }
        AMP::shared_ptr<AMP::Database> elemOp_db = input_db->getDatabase( mechElemDbStr );
        AMP::shared_ptr<AMP::Operator::ElementOperationParameters> elemOpParams(
            new AMP::Operator::ElementOperationParameters( elemOp_db ) );
        AMP::shared_ptr<AMP::Operator::MechanicsLinearElement> mechLinElem(
            new AMP::Operator::MechanicsLinearElement( elemOpParams ) );

        AMP::Discretization::DOFManager::shared_ptr dofMap =
            AMP::Discretization::simpleDOFManager::create(
                meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 3, true );

        AMP_INSIST( input_db->keyExists( "Mechanics_Assembly" ),
                    "Key ''Mechanics_Assembly'' is missing!" );
        AMP::shared_ptr<AMP::Database> mechAssembly_db =
            input_db->getDatabase( "Mechanics_Assembly" );
        AMP::shared_ptr<AMP::Operator::MechanicsLinearFEOperatorParameters> mechOpParams(
            new AMP::Operator::MechanicsLinearFEOperatorParameters( mechAssembly_db ) );
        mechOpParams->d_materialModel = isotropicModel;
        mechOpParams->d_elemOp        = mechLinElem;
        mechOpParams->d_Mesh          = meshAdapter;
        mechOpParams->d_inDofMap      = dofMap;
        mechOpParams->d_outDofMap     = dofMap;
        AMP::shared_ptr<AMP::Operator::MechanicsLinearFEOperator> mechOp(
            new AMP::Operator::MechanicsLinearFEOperator( mechOpParams ) );

        AMP::LinearAlgebra::Variable::shared_ptr mechVariable = mechOp->getOutputVariable();

        AMP_INSIST( input_db->keyExists( "Displacement_Boundary" ),
                    "Key ''Displacement_Boundary'' is missing!" );
        AMP::shared_ptr<AMP::Database> disp_db = input_db->getDatabase( "Displacement_Boundary" );
        AMP::shared_ptr<AMP::Operator::DirichletMatrixCorrectionParameters> dirichletOpParams(
            new AMP::Operator::DirichletMatrixCorrectionParameters( disp_db ) );
        dirichletOpParams->d_inputMatrix = mechOp->getMatrix();
        // This is just the variable used to extract the dof_map.
        // This boundary operator itself has an empty input and output variable
        dirichletOpParams->d_variable = mechVariable;
        dirichletOpParams->d_Mesh     = meshAdapter;
        AMP::shared_ptr<AMP::Operator::DirichletMatrixCorrection> dirichletMatOp(
            new AMP::Operator::DirichletMatrixCorrection( dirichletOpParams ) );

        AMP::LinearAlgebra::Vector::shared_ptr mechSolVec =
            AMP::LinearAlgebra::createVector( dofMap, mechVariable, true );
        AMP::LinearAlgebra::Vector::shared_ptr mechRhsVec = mechSolVec->cloneVector();
        AMP::LinearAlgebra::Vector::shared_ptr mechResVec = mechSolVec->cloneVector();

        for ( int i = 0; i < 3; i++ ) {
            mechSolVec->setRandomValues();
            mechRhsVec->setRandomValues();
            mechResVec->setRandomValues();
            mechOp->residual( mechRhsVec, mechSolVec, mechResVec );
        } // end for i

        mechOp->reset( mechOpParams );

        ut->passes( exeName + " : " + mechElemDbStr );

    } // end for useReduced
}

int testLinearMechanics_reset( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    myTest( &ut );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
