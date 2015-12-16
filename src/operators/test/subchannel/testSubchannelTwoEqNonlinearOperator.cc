#include "utils/AMPManager.h"
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

#include "vectors/VectorBuilder.h"

#include "operators/OperatorBuilder.h"
#include "operators/subchannel/SubchannelConstants.h"
#include "operators/subchannel/SubchannelOperatorParameters.h"
#include "operators/subchannel/SubchannelPhysicsModel.h"
#include "operators/subchannel/SubchannelTwoEqNonlinearOperator.h"

#include "ampmesh/StructuredMeshHelper.h"
#include "discretization/simpleDOF_Manager.h"
#include "discretization/structuredFaceDOFManager.h"

void Test( AMP::UnitTest *ut, std::string exeName )
{
    // create input and output file names
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;
    AMP::PIO::logOnlyNodeZero( log_file );

    // get input database from input file
    AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );
    input_db->printClassData( AMP::plog );

    // create mesh
    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    AMP::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase( "Mesh" );
    AMP::shared_ptr<AMP::Mesh::MeshParameters> meshParams(
        new AMP::Mesh::MeshParameters( mesh_db ) );
    meshParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    AMP::shared_ptr<AMP::Mesh::Mesh> subchannelMesh = AMP::Mesh::Mesh::buildMesh( meshParams );
    AMP::Mesh::Mesh::shared_ptr xyFaceMesh;
    xyFaceMesh = subchannelMesh->Subset(
        AMP::Mesh::StructuredMeshHelper::getXYFaceIterator( subchannelMesh, 0 ) );

    // get dof manager
    int DofsPerFace[3] = { 0, 0, 2 };
    AMP::Discretization::DOFManager::shared_ptr faceDOFManager =
        AMP::Discretization::structuredFaceDOFManager::create( subchannelMesh, DofsPerFace, 1 );

    // get input and output variables
    AMP::LinearAlgebra::Variable::shared_ptr inputVariable(
        new AMP::LinearAlgebra::Variable( "flow" ) );
    AMP::LinearAlgebra::Variable::shared_ptr outputVariable(
        new AMP::LinearAlgebra::Variable( "flow" ) );

    // create solution, rhs, and residual vectors
    AMP::LinearAlgebra::Vector::shared_ptr SolVec =
        AMP::LinearAlgebra::createVector( faceDOFManager, inputVariable, true );
    AMP::LinearAlgebra::Vector::shared_ptr RhsVec =
        AMP::LinearAlgebra::createVector( faceDOFManager, outputVariable, true );
    AMP::LinearAlgebra::Vector::shared_ptr ResVec =
        AMP::LinearAlgebra::createVector( faceDOFManager, outputVariable, true );

    // create subchannel physics model
    AMP::shared_ptr<AMP::Database> subchannelPhysics_db =
        input_db->getDatabase( "SubchannelPhysicsModel" );
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModelParameters> params(
        new AMP::Operator::ElementPhysicsModelParameters( subchannelPhysics_db ) );
    AMP::shared_ptr<AMP::Operator::SubchannelPhysicsModel> subchannelPhysicsModel(
        new AMP::Operator::SubchannelPhysicsModel( params ) );

    // get nonlinear operator database
    AMP::shared_ptr<AMP::Database> subchannelOperator_db =
        input_db->getDatabase( "SubchannelTwoEqNonlinearOperator" );
    // set operator parameters
    AMP::shared_ptr<AMP::Operator::SubchannelOperatorParameters> subchannelOpParams(
        new AMP::Operator::SubchannelOperatorParameters( subchannelOperator_db ) );
    subchannelOpParams->d_Mesh                   = subchannelMesh;
    subchannelOpParams->d_subchannelPhysicsModel = subchannelPhysicsModel;
    subchannelOpParams->clad_x = input_db->getDatabase( "CladProperties" )->getDoubleArray( "x" );
    subchannelOpParams->clad_y = input_db->getDatabase( "CladProperties" )->getDoubleArray( "y" );
    subchannelOpParams->clad_d = input_db->getDatabase( "CladProperties" )->getDoubleArray( "d" );
    AMP::shared_ptr<AMP::Operator::SubchannelTwoEqNonlinearOperator> subchannelOperator(
        new AMP::Operator::SubchannelTwoEqNonlinearOperator( subchannelOpParams ) );

    // report successful creation
    ut->passes( exeName + ": creation" );
    std::cout.flush();

    // reset the nonlinear operator
    subchannelOperator->reset( subchannelOpParams );
    std::vector<size_t> dofs;
    AMP::Mesh::MeshIterator face     = xyFaceMesh->getIterator( AMP::Mesh::Face, 0 );
    AMP::Mesh::MeshIterator end_face = face.end();

    const double h_scale = AMP::Operator::Subchannel::scaleEnthalpy;
    const double P_scale = AMP::Operator::Subchannel::scalePressure;

    {
        // Test apply with known residual evaluation
        face = xyFaceMesh->getIterator( AMP::Mesh::Face, 0 );
        faceDOFManager->getDOFs( face->globalID(), dofs );
        ++face;
        SolVec->setValueByGlobalID( dofs[0], h_scale * 700.0e3 );
        SolVec->setValueByGlobalID( dofs[1], P_scale * 12.4e6 );
        faceDOFManager->getDOFs( face->globalID(), dofs );
        ++face;
        SolVec->setValueByGlobalID( dofs[0], h_scale * 900.0e3 );
        SolVec->setValueByGlobalID( dofs[1], P_scale * 12.3e6 );
        faceDOFManager->getDOFs( face->globalID(), dofs );
        ++face;
        SolVec->setValueByGlobalID( dofs[0], h_scale * 800.0e3 );
        SolVec->setValueByGlobalID( dofs[1], P_scale * 16.2e6 );
        faceDOFManager->getDOFs( face->globalID(), dofs );
        ++face;
        SolVec->setValueByGlobalID( dofs[0], h_scale * 650.0e3 );
        SolVec->setValueByGlobalID( dofs[1], P_scale * 14.1e5 );
        faceDOFManager->getDOFs( face->globalID(), dofs );
        ++face;
        SolVec->setValueByGlobalID( dofs[0], h_scale * 367.4e3 );
        SolVec->setValueByGlobalID( dofs[1], P_scale * 31.5e5 );
        faceDOFManager->getDOFs( face->globalID(), dofs );
        ++face;
        SolVec->setValueByGlobalID( dofs[0], h_scale * 657.2e3 );
        SolVec->setValueByGlobalID( dofs[1], P_scale * 12.5e6 );
        faceDOFManager->getDOFs( face->globalID(), dofs );
        ++face;
        SolVec->setValueByGlobalID( dofs[0], h_scale * 788.5e3 );
        SolVec->setValueByGlobalID( dofs[1], P_scale * 12.7e6 );
        faceDOFManager->getDOFs( face->globalID(), dofs );
        ++face;
        SolVec->setValueByGlobalID( dofs[0], h_scale * 235.7e2 );
        SolVec->setValueByGlobalID( dofs[1], P_scale * 17.8e6 );
        faceDOFManager->getDOFs( face->globalID(), dofs );
        ++face;
        SolVec->setValueByGlobalID( dofs[0], h_scale * 673.1e3 );
        SolVec->setValueByGlobalID( dofs[1], P_scale * 13.6e6 );
        faceDOFManager->getDOFs( face->globalID(), dofs );
        ++face;
        SolVec->setValueByGlobalID( dofs[0], h_scale * 385.2e3 );
        SolVec->setValueByGlobalID( dofs[1], P_scale * 16.3e6 );

        subchannelOperator->apply( SolVec, ResVec );
        bool passedKnownTest = true;
        double known[20]     = {
            -565469.235701075, -8.34485290710296, 148482.09384907,   343.485735068754,
            -91004.5857613742, -1298.94510822818, -136965.213005729, 153.631662269915,
            -243658.201123872, 822.509798734236,  193771.812083809,  18.4080923950774,
            73891.3070728496,  448.83209941079,   -608747.573661467, -368.184663947199,
            484044.660140265,  237.580675983725,  -225841.184839454, 800000
        };
        face  = xyFaceMesh->getIterator( AMP::Mesh::Face, 0 );
        int i = 0;
        for ( ; face != end_face; ++face, ++i ) {
            faceDOFManager->getDOFs( face->globalID(), dofs );
            double h_val = ResVec->getValueByGlobalID( dofs[0] ) / h_scale;
            double p_val = ResVec->getValueByGlobalID( dofs[1] ) / P_scale;
            if ( !AMP::Utilities::approx_equal( h_val, known[2 * i], 0.01 ) ) {
                passedKnownTest = false;
                AMP::pout << "Calculated: " << h_val << ", Known: " << known[2 * i] << "\n";
            }
            if ( !AMP::Utilities::approx_equal( p_val, known[2 * i + 1], 0.01 ) ) {
                passedKnownTest = false;
                AMP::pout << "Calculated: " << p_val << ", Known: " << known[2 * i + 1] << "\n";
            }
        }
        if ( passedKnownTest )
            ut->passes( exeName + ": known value test" );
        else
            ut->failure( exeName + ": known residual test" );
    }
}

int main( int argc, char *argv[] )
{
    AMP::AMPManagerProperties startup_properties;
    startup_properties.use_MPI_Abort = false;
    AMP::AMPManager::startup( argc, argv, startup_properties );

    AMP::UnitTest ut;

    const int NUMFILES          = 1;
    std::string files[NUMFILES] = { "testSubchannelTwoEqNonlinearOperator" };

    for ( auto &file : files ) {
        try {
            Test( &ut, file );
        } catch ( std::exception &err ) {
            std::cout << "ERROR: While testing " << argv[0] << err.what() << std::endl;
            ut.failure( "ERROR: While testing: " + file );
        } catch ( ... ) {
            std::cout << "ERROR: While testing " << argv[0] << "An unknown exception was thrown."
                      << std::endl;
            ut.failure( "ERROR: While testing: " + file );
        }
    }

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
