#include "AMP/IO/PIO.h"
#include "AMP/IO/Writer.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/mesh/Mesh.h"
#include "AMP/mesh/MeshFactory.h"
#include "AMP/mesh/MeshParameters.h"
#include "AMP/operators/libmesh/GradientOperator.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/vectors/VectorBuilder.h"

#include <memory>


static void testGrad( std::shared_ptr<AMP::Operator::GradientOperator> op,
                      std::shared_ptr<AMP::LinearAlgebra::Vector> rhsVec,
                      std::shared_ptr<AMP::LinearAlgebra::Vector> solVec,
                      std::shared_ptr<AMP::LinearAlgebra::Vector> ansVec,
                      std::function<double( AMP::Mesh::Point )> f,
                      std::function<AMP::Mesh::Point( AMP::Mesh::Point )> g )
{
    auto mesh = op->getMesh();

    // Initialize rhs/ans
    rhsVec->zero();
    solVec->zero();
    auto rhsDOFManager = rhsVec->getDOFManager();
    auto ansDOFManager = ansVec->getDOFManager();
    std::vector<size_t> dofs1, dofs3;
    for ( auto &node : mesh->getIterator( AMP::Mesh::GeomType::Vertex ) ) {
        auto id = node.globalID();
        auto p  = node.coord();
        auto f0 = f( p );
        auto g0 = g( p );
        rhsDOFManager->getDOFs( id, dofs1 );
        ansDOFManager->getDOFs( id, dofs3 );
        AMP_ASSERT( dofs1.size() == 1 );
        AMP_ASSERT( dofs3.size() == 3 );
        rhsVec->setValuesByGlobalID( dofs1.size(), dofs1.data(), &f0 );
        ansVec->setValuesByGlobalID( dofs3.size(), dofs3.data(), g0.data() );
    }
    rhsVec->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );
    ansVec->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );

    // Apply
    op->apply( rhsVec, solVec );
}


static void run( const std::string &input_file, AMP::UnitTest &ut )
{
    // Load the input file
    auto input_db = AMP::Database::parseInputFile( input_file );

    // Get the Mesh database and create the mesh parameters
    auto database   = input_db->getDatabase( "Mesh" );
    auto meshParams = std::make_shared<AMP::Mesh::MeshParameters>( database );
    meshParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );

    // Create the meshes from the input database
    auto mesh = AMP::Mesh::MeshFactory::create( meshParams );

    // Create the gradient operator
    auto grad_db   = input_db->getDatabase( "GradientOperator" );
    auto params    = std::make_shared<AMP::Operator::OperatorParameters>( grad_db );
    params->d_Mesh = mesh;
    auto op        = std::make_shared<AMP::Operator::GradientOperator>( params );

    // Create the vectors
    auto inputVar    = op->getInputVariable();
    auto outputVar   = op->getOutputVariable();
    auto nodalDofMap = AMP::Discretization::simpleDOFManager::create(
        mesh, AMP::Mesh::GeomType::Vertex, 1, 1, true );
    auto nodalVectorDofMap = AMP::Discretization::simpleDOFManager::create(
        mesh, AMP::Mesh::GeomType::Vertex, 1, 3, true );
    auto rhsVec = AMP::LinearAlgebra::createVector( nodalDofMap, inputVar, true );
    auto solVec = AMP::LinearAlgebra::createVector( nodalVectorDofMap, outputVar, true );
    auto ansVec = AMP::LinearAlgebra::createVector( nodalVectorDofMap, outputVar, true );
    auto errVec = AMP::LinearAlgebra::createVector( nodalVectorDofMap, outputVar, true );

    // Create the writer
    auto writer = AMP::IO::Writer::buildWriter( "silo" );
    writer->registerMesh( mesh );
    writer->registerVector( rhsVec, mesh, AMP::Mesh::GeomType::Vertex, "f" );
    writer->registerVector( solVec, mesh, AMP::Mesh::GeomType::Vertex, "g" );
    writer->registerVector( ansVec, mesh, AMP::Mesh::GeomType::Vertex, "g0" );
    writer->registerVector( errVec, mesh, AMP::Mesh::GeomType::Vertex, "err" );

    // Create the functions to test
    using AMP::Mesh::Point;
    std::vector<std::string> name;
    std::vector<std::function<double( Point )>> f;
    std::vector<std::function<Point( Point )>> g;
    std::vector<double> tol;

    // Add a linear function
    name.emplace_back( "linear" );
    f.emplace_back( []( Point p ) { return 1.1 + 2.2 * p.x() + 3.3 * p.y() + 4.4 * p.z(); } );
    g.emplace_back( []( Point ) { return Point( 2.2, 3.3, 4.4 ); } );
    tol.push_back( 1e-12 );

    // Add a simple quadratic function
    name.emplace_back( "quadratic" );
    f.emplace_back( []( Point p ) { return p.x() * p.x() + p.y() * p.y() + p.z() * p.z(); } );
    g.emplace_back( []( Point p ) { return Point( 2 * p.x(), 2 * p.y(), 2 * p.z() ); } );
    tol.push_back( 0.05 );

    // Run the tests
    for ( size_t i = 0; i < f.size(); i++ ) {
        testGrad( op, rhsVec, solVec, ansVec, f[i], g[i] );
        errVec->subtract( *solVec, *ansVec );
        double err = static_cast<double>( errVec->L2Norm() / ansVec->L2Norm() );
        if ( err < tol[i] ) {
            ut.passes( name[i] );
        } else {
            printf( "Error for %s exceeded tol: %e (%e)\n", name[i].data(), err, tol[i] );
            std::cout << "  norm(errVec) = " << errVec->L2Norm() << std::endl;
            std::cout << "  norm(ansVec) = " << ansVec->L2Norm() << std::endl;
            ut.failure( name[i] );
        }
        writer->writeFile( "testGradientOperator", i, 0 );
    }
}


int testGradientOperator( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;
    AMP_ASSERT( argc >= 2 );

    for ( int i = 1; i < argc; i++ )
        run( argv[i], ut );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
