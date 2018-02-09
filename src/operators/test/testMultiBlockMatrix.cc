#include "AMP/materials/Material.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/InputDatabase.h"
#include "AMP/utils/InputManager.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/shared_ptr.h"
#include "AMP/vectors/Variable.h"
#include <string>

#include "AMP/ampmesh/Mesh.h"
#include "AMP/utils/AsciiWriter.h"
#include "AMP/utils/Writer.h"

#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"

#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/libmesh/MassLinearFEOperator.h"

#include "AMP/operators/LinearBVPOperator.h"


#define ITFAILS ut.failure( __LINE__ );
#define UNIT_TEST( a ) \
    if ( !( a ) )      \
        ut.failure( __LINE__ );

void LinearTimeOperatorTest( AMP::UnitTest *ut )
{
    std::string input_file = "input_testMultiBlockMatrix";
    std::string log_file   = "log_testMultiBlockMatrix";

    AMP::PIO::logOnlyNodeZero( log_file );

    AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );
    input_db->printClassData( AMP::plog );

    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    AMP::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase( "Mesh" );
    AMP::shared_ptr<AMP::Mesh::MeshParameters> mgrParams(
        new AMP::Mesh::MeshParameters( mesh_db ) );
    mgrParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    AMP::shared_ptr<AMP::Mesh::Mesh> meshAdapter = AMP::Mesh::Mesh::buildMesh( mgrParams );

    //--------------------------------------------------
    // Create a DOF manager for a nodal vector
    //--------------------------------------------------
    int DOFsPerNode     = 1;
    int nodalGhostWidth = 1;
    bool split          = true;

    AMP::Discretization::DOFManager::shared_ptr nodalDofMap =
        AMP::Discretization::simpleDOFManager::create(
            meshAdapter, AMP::Mesh::GeomType::Vertex, nodalGhostWidth, DOFsPerNode, split );

    //----------------------------------------------------------------------------------------
    // create a linear BVP operator
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> elementModel;

    AMP::shared_ptr<AMP::Operator::LinearBVPOperator> linearOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "LinearOperator", input_db, elementModel ) );

    // ---------------------------------------------------------------------------------------
    // create a mass linear BVP operator
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> massElementModel;
    AMP::shared_ptr<AMP::Operator::LinearBVPOperator> massOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "MassLinearOperator", input_db, massElementModel ) );

    // ---------------------------------------------------------------------------------------
    AMP::LinearAlgebra::Matrix::shared_ptr fullMat = linearOperator->getMatrix();
    AMP::Utilities::AsciiWriter fullMatWriter;
    fullMatWriter.registerMatrix( fullMat );
    fullMatWriter.writeFile( "FullMat", 0 );

    AMP::LinearAlgebra::Matrix::shared_ptr massMat = massOperator->getMatrix();
    AMP::LinearAlgebra::Matrix::shared_ptr diffMat = linearOperator->getMatrix();

    AMP::LinearAlgebra::Matrix::shared_ptr sinMat = diffMat->cloneMatrix();
    sinMat->makeConsistent();
    sinMat->zero();

    sinMat->axpy( 1.0, diffMat );
    sinMat->makeConsistent();
    //    diffMat->axpy(1.0/0.01, massMat);
    //    diffMat->makeConsistent();

    AMP::Utilities::AsciiWriter sinMatWriter;
    sinMatWriter.registerMatrix( sinMat );
    sinMatWriter.writeFile( "SinMat", 0 );

    ut->passes( "Ran to completion" );
}


//---------------------------------------------------------------------------//

int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    LinearTimeOperatorTest( &ut );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}

//---------------------------------------------------------------------------//
//                        end of SundialsVectorTest.cc
//---------------------------------------------------------------------------//
