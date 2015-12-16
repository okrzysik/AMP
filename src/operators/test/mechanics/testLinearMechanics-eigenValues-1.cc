
#include "utils/AMPManager.h"
#include "utils/AMP_MPI.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include <iostream>
#include <string>

#include "utils/shared_ptr.h"

#include "utils/AMPManager.h"
#include "utils/AMP_MPI.h"
#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"

#include "ampmesh/Mesh.h"
#include "ampmesh/libmesh/libMesh.h"

#include "discretization/simpleDOF_Manager.h"

#include "libmesh/mesh_generation.h"

#include "operators/mechanics/IsotropicElasticModel.h"
#include "operators/mechanics/MechanicsLinearElement.h"
#include "operators/mechanics/MechanicsLinearFEOperator.h"

#include "operators/boundary/DirichletMatrixCorrection.h"
#include "operators/boundary/DirichletVectorCorrection.h"


void myTest( AMP::UnitTest *ut )
{
    std::string exeName( "testLinearMechanics-eigenValues-1" );
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::PIO::logOnlyNodeZero( log_file );

    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
    int npes = globalComm.getSize();

    if ( npes == 1 ) {
        AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );
        AMP::InputManager::getManager()->parseInputFile( input_file, input_db );
        input_db->printClassData( AMP::plog );

        AMP_INSIST( input_db->keyExists( "OutputFileName" ), "Key ''OutputFileName'' is missing!" );
        std::string outFileName = input_db->getString( "OutputFileName" );

        FILE *fp;
        fp = fopen( outFileName.c_str(), "w" );
        fprintf( fp, "clc; \n clear; \n A = zeros(24, 24); \n \n" );

        AMP_INSIST( input_db->keyExists( "DISTORT_ELEMENT" ),
                    "Key ''DISTORT_ELEMENT'' is missing!" );
        bool distortElement = input_db->getBool( "DISTORT_ELEMENT" );

        AMP::shared_ptr<::Mesh> mesh( new ::Mesh( 3 ) );
        MeshTools::Generation::build_cube(
            ( *( mesh.get() ) ), 1, 1, 1, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, HEX8, false );

        if ( distortElement ) {
            ::Elem *elemPtr = mesh->elem( 0 );

            ( elemPtr->point( 0 ) )( 0 ) -= 0.1;
            ( elemPtr->point( 0 ) )( 1 ) -= 0.2;
            ( elemPtr->point( 0 ) )( 2 ) -= 0.3;

            ( elemPtr->point( 3 ) )( 0 ) -= 0.2;
            ( elemPtr->point( 3 ) )( 1 ) += 0.2;
            ( elemPtr->point( 3 ) )( 2 ) -= 0.1;

            ( elemPtr->point( 6 ) )( 0 ) += 0.3;
            ( elemPtr->point( 6 ) )( 1 ) += 0.2;
            ( elemPtr->point( 6 ) )( 2 ) += 0.1;
        }

        AMP::Mesh::Mesh::shared_ptr meshAdapter =
            AMP::Mesh::Mesh::shared_ptr( new AMP::Mesh::libMesh( mesh, "TestMesh" ) );

        AMP_INSIST( input_db->keyExists( "Isotropic_Model" ),
                    "Key ''Isotropic_Model'' is missing!" );
        AMP::shared_ptr<AMP::Database> matModel_db = input_db->getDatabase( "Isotropic_Model" );
        AMP::shared_ptr<AMP::Operator::MechanicsMaterialModelParameters> matModelParams(
            new AMP::Operator::MechanicsMaterialModelParameters( matModel_db ) );
        AMP::shared_ptr<AMP::Operator::IsotropicElasticModel> isotropicModel(
            new AMP::Operator::IsotropicElasticModel( matModelParams ) );

        AMP_INSIST( input_db->keyExists( "Mechanics_Linear_Element" ),
                    "Key ''Mechanics_Linear_Element'' is missing!" );
        AMP::shared_ptr<AMP::Database> elemOp_db =
            input_db->getDatabase( "Mechanics_Linear_Element" );
        AMP::shared_ptr<AMP::Operator::ElementOperationParameters> elemOpParams(
            new AMP::Operator::ElementOperationParameters( elemOp_db ) );
        AMP::shared_ptr<AMP::Operator::MechanicsLinearElement> mechLinElem(
            new AMP::Operator::MechanicsLinearElement( elemOpParams ) );

        AMP::Discretization::DOFManager::shared_ptr dofMap =
            AMP::Discretization::simpleDOFManager::create(
                meshAdapter, AMP::Mesh::Vertex, 1, 3, true );

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

        AMP::shared_ptr<AMP::LinearAlgebra::Matrix> mechMat = mechOp->getMatrix();

        for ( int i = 0; i < 24; ++i ) {
            std::vector<unsigned int> matCols;
            std::vector<double> matVals;
            mechMat->getRowByGlobalID( i, matCols, matVals );
            for ( unsigned int j = 0; j < matCols.size(); j++ ) {
                fprintf(
                    fp, "A(%d, %d) = %.15f ; \n", ( i + 1 ), (int) ( matCols[j] + 1 ), matVals[j] );
            } // end for j
            fprintf( fp, "\n" );
        } // end for i

        AMP::Mesh::MeshIterator nd     = meshAdapter->getIterator( AMP::Mesh::Vertex, 0 );
        AMP::Mesh::MeshIterator end_nd = nd.end();

        for ( int i = 0; nd != end_nd; ++nd, ++i ) {
            std::vector<double> pt = nd->coord();
            fprintf( fp, "nd = %d, x = %.15f, y = %.15f, z = %.15f \n", i, pt[0], pt[1], pt[2] );
            std::vector<size_t> globalIds;
            dofMap->getDOFs( nd->globalID(), globalIds );
            fprintf( fp,
                     "nd = %d, d0 = %d, d1 = %d, d2 = %d \n",
                     i,
                     (int) globalIds[0],
                     (int) globalIds[1],
                     (int) globalIds[2] );
        }

        fprintf( fp, "format long e; \n\n" );
        fprintf( fp, "sort(eig(A)) \n\n" );
        fclose( fp );
    }

    ut->passes( exeName );
}

int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::shared_ptr<AMP::Mesh::initializeLibMesh> libmeshInit(
        new AMP::Mesh::initializeLibMesh( AMP_COMM_WORLD ) );

    AMP::UnitTest ut;
    try {
        myTest( &ut );
    }
    catch ( std::exception &err ) {
        std::cout << "ERROR: While testing " << argv[0] << err.what() << std::endl;
        ut.failure( "ERROR: While testing" );
    }
    catch ( ... ) {
        std::cout << "ERROR: While testing " << argv[0] << "An unknown exception was thrown."
                  << std::endl;
        ut.failure( "ERROR: While testing" );
    }
    ut.report();
    int num_failed = ut.NumFailGlobal();

    libmeshInit.reset();
    AMP::AMPManager::shutdown();
    return num_failed;
}
