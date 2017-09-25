#include "utils/AMPManager.h"
#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "utils/Utilities.h"
#include "utils/shared_ptr.h"
#include "vectors/Variable.h"
#include <string>

#include "utils/Writer.h"
#include "vectors/Vector.h"

#include "ampmesh/Mesh.h"
#include "discretization/simpleDOF_Manager.h"
#include "vectors/VectorBuilder.h"

#include "operators/libmesh/PowerShape.h"


#define ITFAILS ut.failure( __LINE__ );
#define UNIT_TEST( a ) \
    if ( !( a ) )      \
        ut.failure( __LINE__ );

void test_with_shape( AMP::UnitTest *ut )
{

    //--------------------------------------------------
    //  Read Input File.
    //--------------------------------------------------
    std::string exeName( "testPowerShape-3" );
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::PIO::logAllNodes( log_file );

    AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );

    //--------------------------------------------------
    //   Create the Mesh.
    //--------------------------------------------------
    AMP::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase( "Mesh" );
    AMP::shared_ptr<AMP::Mesh::MeshParameters> mgrParams(
        new AMP::Mesh::MeshParameters( mesh_db ) );
    mgrParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    AMP::shared_ptr<AMP::Mesh::Mesh> meshAdapter = AMP::Mesh::Mesh::buildMesh( mgrParams );

    //--------------------------------------------------
    //  Construct PowerShape.
    //--------------------------------------------------
    AMP_INSIST( input_db->keyExists( "MyPowerShape" ), "Key ''MyPowerShape'' is missing!" );
    AMP::shared_ptr<AMP::Database> shape_db = input_db->getDatabase( "MyPowerShape" );
    AMP::shared_ptr<AMP::Operator::PowerShapeParameters> shape_params(
        new AMP::Operator::PowerShapeParameters( shape_db ) );
    shape_params->d_Mesh = meshAdapter;
    AMP::shared_ptr<AMP::Operator::PowerShape> shape(
        new AMP::Operator::PowerShape( shape_params ) );

    // Create a DOF manager for a gauss point vector
    int DOFsPerNode = 8;
    int ghostWidth  = 0;
    bool split      = true;
    AMP::Discretization::DOFManager::shared_ptr dof_map =
        AMP::Discretization::simpleDOFManager::create(
            meshAdapter, AMP::Mesh::GeomType::Volume, ghostWidth, DOFsPerNode, split );

    // Create a shared pointer to a Variable - Power - Output because it will be used in the
    // "residual" location of
    // apply.
    AMP::LinearAlgebra::Variable::shared_ptr SpecificPowerShapeVar(
        new AMP::LinearAlgebra::Variable( "SpecificPowerInWattsPerKg" ) );

    // Create a vector associated with the Variable.
    AMP::LinearAlgebra::Vector::shared_ptr SpecificPowerShapeVec =
        AMP::LinearAlgebra::createVector( dof_map, SpecificPowerShapeVar, split );
    AMP::LinearAlgebra::Vector::shared_ptr SpecificPowerMagnitudeVec =
        SpecificPowerShapeVec->cloneVector();
    SpecificPowerMagnitudeVec->setToScalar( 0.1 );

    // Set the initial value for all nodes of SpecificPowerVec to zero.
    SpecificPowerShapeVec->setToScalar( 0.0 );
    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    try {
        shape->apply( SpecificPowerMagnitudeVec, SpecificPowerShapeVec );
    } catch ( std::exception const &a ) {
        std::cout << a.what() << std::endl;
        ut->failure( "error" );
    }

    ut->passes( "PowerShape gets past apply with a non-flat power shape." );

    // Check that the data is non-negative
    bool itpasses = 1;
    AMP::Mesh::MeshIterator elem =
        meshAdapter->getIterator( AMP::Mesh::GeomType::Volume, ghostWidth );
    AMP::Mesh::MeshIterator end_elems = elem.end();

    for ( ; elem != end_elems; ++elem ) {
        for ( int i = 0; i < DOFsPerNode; i++ ) {
            std::vector<size_t> ndx;
            dof_map->getDOFs( elem->globalID(), ndx );
            int offset = ndx[i];
            if ( SpecificPowerShapeVec->getValueByGlobalID( offset ) < 0.0 ) {
                if ( !itpasses )
                    ut->failure( "PowerShape error" );
                itpasses = 0;
            }
        } // end for gauss-points
    }     // end for elements

    if ( itpasses )
        ut->passes( "PowerShape produces a non-negative power shape." );

    //-----------------------------------------------------------------
    //  Testing the new legendre function. valLegendre(int n, double x)
    //-----------------------------------------------------------------
    double pn = shape->evalLegendre( 3, 2.0 );
    if ( pn != 17. )
        ut->failure( "PowerShape error" );
}


int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    test_with_shape( &ut );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
