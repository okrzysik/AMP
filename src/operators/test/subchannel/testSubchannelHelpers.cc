#include "utils/AMPManager.h"
#include "utils/AMP_MPI.h"
#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"

#include "ampmesh/Mesh.h"
#include "ampmesh/StructuredMeshHelper.h"
#include "discretization/simpleDOF_Manager.h"
#include "discretization/structuredFaceDOFManager.h"
#include "operators/subchannel/SubchannelConstants.h"
#include "operators/subchannel/SubchannelHelpers.h"
#include "vectors/VectorBuilder.h"
#include "vectors/VectorSelector.h"


#include <iostream>
#include <string>

void testSubchannelHelpers( AMP::UnitTest *ut, std::string input_file )
{

    const double pi = 3.1415926535897932;

    // Read the input file
    AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );
    input_db->printClassData( AMP::plog );

    // Get the Mesh database and create the mesh parameters
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
    AMP::shared_ptr<AMP::Database> database = input_db->getDatabase( "Mesh" );
    AMP::shared_ptr<AMP::Mesh::MeshParameters> meshParams(
        new AMP::Mesh::MeshParameters( database ) );
    meshParams->setComm( globalComm );

    // Get the meshes and clad properties
    AMP::Mesh::Mesh::shared_ptr manager        = AMP::Mesh::Mesh::buildMesh( meshParams );
    AMP::Mesh::Mesh::shared_ptr cladMesh       = manager->Subset( "clad" );
    AMP::Mesh::Mesh::shared_ptr subchannelMesh = manager->Subset( "subchannel" );
    std::vector<double> clad_x, clad_y, clad_d;
    AMP::Operator::Subchannel::getCladProperties( globalComm, cladMesh, clad_x, clad_y, clad_d );
    // std::vector<double> box = cladMesh->getBoundingBox();
    bool pass = true;
    for ( size_t i = 0; i < clad_x.size(); i++ ) {
        if ( !AMP::Utilities::approx_equal( clad_d[i], clad_d[0] ) ) pass    = false;
        if ( !AMP::Utilities::approx_equal( clad_d[0], 2 * 0.004705 ) ) pass = false;
    }
    if ( pass )
        ut->passes( "Get clad properties passes basic sanity check" );
    else
        ut->failure( "Get clad properties passes basic sanity check" );

    // Get the subchannel properties
    std::vector<double> x, y, z, area, subchannel_diam, rod_diameter, channel_fraction;
    AMP::Mesh::StructuredMeshHelper::getXYZCoordinates( subchannelMesh, x, y, z );
    AMP::Operator::Subchannel::getSubchannelProperties( subchannelMesh,
                                                        clad_x,
                                                        clad_y,
                                                        clad_d,
                                                        x,
                                                        y,
                                                        area,
                                                        subchannel_diam,
                                                        rod_diameter,
                                                        channel_fraction );
    size_t N_subchannels = ( x.size() - 1 ) * ( y.size() - 1 );
    AMP_ASSERT( area.size() == N_subchannels );
    pass = true;
    for ( size_t i = 0; i < N_subchannels; i++ ) {
        if ( subchannel_diam[i] == 0 ) pass   = false;
        if ( channel_fraction[i] > 1.0 ) pass = false;
        if ( !AMP::Utilities::approx_equal( rod_diameter[i], clad_d[0] ) ) pass = false;
    }
    if ( pass )
        ut->passes( "Get subchannel properties passes basic sanity check" );
    else
        ut->failure( "Get subchannel properties passes basic sanity check" );

    // First test getHeatFluxGeneration
    size_t N = z.size() - 1;
    std::vector<double> dz( N );
    for ( size_t i = 0; i < N; i++ ) dz[i] = z[i + 1] - z[i];
    double diam                            = 0.00822;
    double Q_tot                           = 66.5e3;
    std::vector<double> flux1 =
        AMP::Operator::Subchannel::getHeatFluxGeneration( "Sinusoidal", z, diam, Q_tot );
    std::vector<double> flux2 =
        AMP::Operator::Subchannel::getHeatFluxGeneration( "Flat", z, diam, Q_tot );
    double Q_tot1  = 0.0;
    double Q_tot2  = 0.0;
    bool pass_sine = true;
    bool pass_flat = true;
    for ( size_t i = 0; i < N; i++ ) {
        Q_tot1 += flux1[i] * pi * diam * dz[i];
        Q_tot2 += flux2[i] * pi * diam * dz[i];
        if ( flux2[i] != flux2[0] ) pass_flat = false;
    }
    if ( !AMP::Utilities::approx_equal( Q_tot, Q_tot1 ) ) pass_sine = false;
    if ( !AMP::Utilities::approx_equal( Q_tot, Q_tot2 ) ) pass_flat = false;
    if ( pass_sine )
        ut->passes( "Sinusoidal shape gives correct flux" );
    else
        ut->failure( "Sinusoidal shape gives correct flux" );
    if ( pass_flat )
        ut->passes( "Flat shape gives correct flux" );
    else
        ut->failure( "Flat shape gives correct flux" );

    // Test getHeatFluxClad
    AMP::shared_ptr<AMP::Database> subchannel_db =
        input_db->getDatabase( "SubchannelPhysicsModel" );
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModelParameters> params(
        new AMP::Operator::ElementPhysicsModelParameters( subchannel_db ) );
    AMP::shared_ptr<AMP::Operator::SubchannelPhysicsModel> subchannelPhysicsModel(
        new AMP::Operator::SubchannelPhysicsModel( params ) );
    double reynolds = subchannel_db->getDatabase( "Defaults" )->getDouble( "reynolds" );
    double prandtl  = subchannel_db->getDatabase( "Defaults" )->getDouble( "prandtl" );
    AMP::LinearAlgebra::Vector::shared_ptr flowVec, cladTemp;
    if ( subchannelMesh.get() != NULL ) {
        int DOFsPerFace[3] = { 0, 0, 2 };
        AMP::Discretization::DOFManager::shared_ptr flowDOF =
            AMP::Discretization::structuredFaceDOFManager::create( subchannelMesh, DOFsPerFace, 1 );
        DOFsPerFace[2] = 1;
        AMP::Discretization::DOFManager::shared_ptr cladDOF =
            AMP::Discretization::structuredFaceDOFManager::create( subchannelMesh, DOFsPerFace, 1 );
        AMP::LinearAlgebra::Variable::shared_ptr flowVariable(
            new AMP::LinearAlgebra::Variable( "Flow" ) );
        AMP::LinearAlgebra::Variable::shared_ptr thermalVariable(
            new AMP::LinearAlgebra::Variable( "Temperature" ) );
        flowVec          = AMP::LinearAlgebra::createVector( flowDOF, flowVariable );
        cladTemp         = AMP::LinearAlgebra::createVector( cladDOF, thermalVariable );
        double clad_temp = 632;
        double flow_temp = 570;
        double pressure  = 1.6e7;
        std::map<std::string, AMP::shared_ptr<std::vector<double>>> enthalpyArgMap;
        enthalpyArgMap.insert( std::make_pair(
            "temperature",
            AMP::shared_ptr<std::vector<double>>( new std::vector<double>( 1, flow_temp ) ) ) );
        enthalpyArgMap.insert( std::make_pair(
            "pressure",
            AMP::shared_ptr<std::vector<double>>( new std::vector<double>( 1, pressure ) ) ) );
        std::vector<double> enthalpyResult( 1 );
        subchannelPhysicsModel->getProperty( "Enthalpy", enthalpyResult, enthalpyArgMap );
        AMP::LinearAlgebra::Vector::shared_ptr subchannelEnthalpy =
            flowVec->select( AMP::LinearAlgebra::VS_Stride( 0, 2 ), "H" );
        AMP::LinearAlgebra::Vector::shared_ptr subchannelPressure =
            flowVec->select( AMP::LinearAlgebra::VS_Stride( 1, 2 ), "P" );
        subchannelEnthalpy->setToScalar( AMP::Operator::Subchannel::scaleEnthalpy *
                                         enthalpyResult[0] );
        subchannelPressure->setToScalar( AMP::Operator::Subchannel::scalePressure * pressure );
        cladTemp->setToScalar( clad_temp );
    }
    pass = true;
    for ( size_t i = 0; i < N_subchannels; i++ ) {
        AMP_INSIST( N_subchannels == 1, "This is currently only programmed for one subchannel" );
        AMP::Mesh::MeshIterator it =
            AMP::Mesh::StructuredMeshHelper::getXYFaceIterator( subchannelMesh, 0 );
        std::vector<AMP::Mesh::MeshElementID> face_ids( it.size() );
        std::vector<double> z_face( it.size() );
        for ( size_t j = 0; j < it.size(); j++ ) {
            face_ids[j] = it->globalID();
            z_face[j]   = it->centroid()[2];
            ++it;
        }
        double Q_tot3 = 0.0;
        std::vector<double> flux3 =
            AMP::Operator::Subchannel::getHeatFluxClad( z_face,
                                                        face_ids,
                                                        subchannel_diam[i],
                                                        reynolds,
                                                        prandtl,
                                                        channel_fraction[i],
                                                        subchannelPhysicsModel,
                                                        flowVec,
                                                        cladTemp );
        AMP_ASSERT( flux3.size() == N );
        for ( size_t i = 0; i < N; i++ ) Q_tot3 += flux3[i] * pi * diam * dz[i];
        if ( fabs( Q_tot3 - Q_tot ) / Q_tot > 0.2 ) pass = false;
    }
    if ( pass )
        ut->passes( "Clad shape gives ~ correct flux" );
    else
        ut->failure( "Clad shape gives ~ correct flux" );
}


int main( int argc, char *argv[] )
{
    AMP::AMPManagerProperties startup_properties;
    startup_properties.use_MPI_Abort = false;
    AMP::AMPManager::startup( argc, argv, startup_properties );

    AMP::UnitTest ut;

    std::string input = "input_testSubchannelUtilities";
    testSubchannelHelpers( &ut, input );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
