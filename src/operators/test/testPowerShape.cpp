#include "AMP/IO/PIO.h"
#include "AMP/IO/Writer.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/mesh/Mesh.h"
#include "AMP/mesh/MeshParameters.h"
#include "AMP/operators/libmesh/PowerShape.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"

#include <memory>
#include <string>


static void test_with_shape( AMP::UnitTest *ut, const std::string &exeName )
{
    //  Read Input File
    std::cout << "Testing " << exeName << std::endl;
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::logAllNodes( log_file );

    auto input_db = AMP::Database::parseInputFile( input_file );

    //   Create the Mesh
    auto mesh_db   = input_db->getDatabase( "Mesh" );
    auto mgrParams = std::make_shared<AMP::Mesh::MeshParameters>( mesh_db );
    mgrParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    auto meshAdapter = AMP::Mesh::Mesh::buildMesh( mgrParams );

    //  Construct PowerShape.
    AMP_INSIST( input_db->keyExists( "MyPowerShape" ), "Key ''MyPowerShape'' is missing!" );
    auto shape_db        = input_db->getDatabase( "MyPowerShape" );
    auto shape_params    = std::make_shared<AMP::Operator::PowerShapeParameters>( shape_db );
    shape_params->d_Mesh = meshAdapter;
    auto shape           = std::make_shared<AMP::Operator::PowerShape>( shape_params );

    // Create a DOF manager for a gauss point vector
    int DOFsPerNode = 8;
    int ghostWidth  = 1;
    bool split      = true;
    auto dof_map    = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Volume, ghostWidth, DOFsPerNode, split );

    // Create a shared pointer to a Variable - Power - Output because it will be used in the
    // "residual" location of apply
    auto SpecificPowerShapeVar =
        std::make_shared<AMP::LinearAlgebra::Variable>( "SpecificPowerInWattsPerKg" );

    // Create a vector associated with the Variable
    auto SpecificPowerShapeVec =
        AMP::LinearAlgebra::createVector( dof_map, SpecificPowerShapeVar, split );
    auto SpecificPowerMagnitudeVec = SpecificPowerShapeVec->cloneVector();
    SpecificPowerMagnitudeVec->setToScalar( 4157. );

    // Set the initial value for all nodes of SpecificPowerVec to zero
    SpecificPowerShapeVec->setToScalar( 0.0 );
    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    try {
        shape->apply( SpecificPowerMagnitudeVec, SpecificPowerShapeVec );
    } catch ( std::exception const &a ) {
        std::cout << a.what() << std::endl;
        ut->failure( exeName + ": exception" );
    }

    ut->passes( exeName + ": PowerShape gets past apply with a non-flat power shape." );

#ifdef USE_EXT_SILO
    auto siloWriter = AMP::IO::Writer::buildWriter( "Silo" );
    siloWriter->registerMesh( meshAdapter );
    siloWriter->registerVector( SpecificPowerShapeVec,
                                meshAdapter,
                                AMP::Mesh::GeomType::Volume,
                                "SpecificPowerInWattsPerKg" );
    siloWriter->writeFile( input_file, 0 );
#endif

    AMP::pout << "SpecificPowerShapeVec->max()"
              << " : " << SpecificPowerShapeVec->min() << " : " << SpecificPowerShapeVec->max()
              << std::endl;
    // Check that the data is non-negative
    bool itpasses  = true;
    auto elem      = meshAdapter->getIterator( AMP::Mesh::GeomType::Volume, ghostWidth );
    auto end_elems = elem.end();

    for ( ; elem != end_elems; ++elem ) {
        for ( int i = 0; i < DOFsPerNode; i++ ) {
            std::vector<size_t> ndx;
            dof_map->getDOFs( elem->globalID(), ndx );
            int offset = ndx[i];
            if ( SpecificPowerShapeVec->getValueByGlobalID( offset ) < 0.0 ) {
                if ( !itpasses )
                    ut->failure( exeName + ": PowerShape error" );
                itpasses = false;
            }
        } // end for gauss-points
    }     // end for elements

    if ( itpasses )
        ut->passes( exeName + ": PowerShape produces a non-negative power shape." );

    //  Testing the new legendre function. valLegendre(int n, double x)
    double pn = shape->evalLegendre( 3, 2.0 );
    if ( pn != 17. )
        ut->failure( exeName + ": PowerShape error" );
}


int testPowerShape( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    auto exeNames = { "testPowerShape-1",
                      "testPowerShape-2",
                      "testPowerShape-3",
                      "testPowerShape-5",
                      "testPowerShape-diffusion" };
    for ( const auto &exeName : exeNames )
        test_with_shape( &ut, exeName );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
