//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   operators/test/testMoabBasedOperator.cc
 * \brief  This tests the Moab iMesh interface with MoabBasedOperator
 */
//---------------------------------------------------------------------------//

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// AMP Moab Includes
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/operators/moab/MoabBasedOperator.h"
#include "AMP/operators/moab/MoabMapOperator.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/Writer.h"
#include "AMP/vectors/VectorBuilder.h"

// MOAB Includes
#include "Coupler.hpp"
#include "MBParallelComm.hpp"
#include "MBRange.hpp"
#include "MBiMesh.hpp"
#include "iMesh.h"
#include "moab/Interface.hpp"


// Helper Class
class MoabDummyOperator : public AMP::Operator::MoabBasedOperator
{
public:
    explicit MoabDummyOperator(
        std::shared_ptr<AMP::Operator::MoabBasedOperatorParameters> &moabParams )
        : AMP::Operator::MoabBasedOperator( moabParams )
    {
        // Create iMesh instance
        iMesh_Instance mbMesh;
        std::string options;
        int ierr;
        iMesh_newMesh( options.c_str(), &mbMesh, &ierr, options.length() );

        // Convert iMesh to MBiMesh
        MBiMesh *mbimesh = reinterpret_cast<MBiMesh *>( mbMesh );

        // Get Moab interface from MBiMesh
        d_moabInterface = mbimesh->mbImpl;

        // Get root set
        iBase_EntitySetHandle root;
        iMesh_createEntSet( mbMesh, 0, &root, &ierr );

        // Build ParallelComm and get index
        moab::ParallelComm *mbComm = new moab::ParallelComm( d_moabInterface );
        int index                  = mbComm->get_id();

        // Set read options and load file
        std::string newReadOpts;
        std::ostringstream extraOpt;
        extraOpt << ";PARALLEL_COMM=" << index;
        std::string readOpts =
            "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARTITION_DISTRIBUTE;PARALLEL_RESOLVE_"
            "SHARED_ENTS;PARALLEL_GHOSTS=3.0.1;CPUTIME";
        newReadOpts          = readOpts + extraOpt.str();
        std::string filename = moabParams->d_db->getString( "moabMeshName" );
        AMP::plog << "Moab mesh name is " << filename << std::endl;
        moab::ErrorCode result = d_moabInterface->load_file(
            filename.c_str(), (EntityHandle *) &root, newReadOpts.c_str() );

        AMP_INSIST( result == MB_SUCCESS, "File not loaded correctly" );

        // Extract nodes
        iBase_EntityHandle *nodes;
        int nodes_alloc = 0;
        int nodes_size;
        iMesh_getEntities( mbMesh,
                           root,
                           iBase_VERTEX,
                           iMesh_ALL_TOPOLOGIES,
                           &nodes,
                           &nodes_alloc,
                           &nodes_size,
                           &ierr );

        // Get list of node coordinates
        std::vector<double> myCoords;
        d_moabInterface->get_vertex_coordinates( myCoords );
        int num_nodes = myCoords.size() / 3;

        AMP_INSIST( num_nodes == nodes_size,
                    "Number of nodes must match number of vertex coordinates" );

        AMP::pout << "Found " << num_nodes << " nodes" << std::endl;

        // Put data on mesh
        std::vector<double> myTemps( num_nodes, -1.0 );
        for ( int i = 0; i < num_nodes; ++i ) {
            // Set temperature to T = x + z
            myTemps[i] = myCoords[i] + myCoords[i + 2 * num_nodes];
        }

        // Add temperature tag
        iBase_TagHandle tempTagHandle;
        std::string tempTagName = "TEMPERATURE";
        iMesh_createTag( mbMesh,
                         tempTagName.c_str(),
                         1,
                         MB_TYPE_DOUBLE,
                         &tempTagHandle,
                         &ierr,
                         tempTagName.length() );

        // Assign data to tag
        iMesh_setDblArrData(
            mbMesh, nodes, nodes_size, tempTagHandle, &myTemps[0], myTemps.size(), &ierr );

        // Free the comm
        delete mbComm;
    }

    void apply( AMP::LinearAlgebra::Vector::const_shared_ptr f,
                AMP::LinearAlgebra::Vector::const_shared_ptr u,
                AMP::LinearAlgebra::Vector::shared_ptr r,
                double a,
                double b )
    {
        /* Don't need an apply for this operator */
    }
    MoabOpParams void finalize(){ /* ... */ };
};


//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

static void moabInterface( AMP::UnitTest *ut )
{
    // Print out AMP banner
    AMP::Utilities::printBanner();

    // Log all nodes
    std::string exeName     = "testMoabBasedOperator";
    std::string input_file  = "input_" + exeName;
    std::string output_file = "output_" + exeName;
    AMP::logAllNodes( output_file );

    // Read Input File.
    auto input_db = AMP::Database::parseInputFile( input_file );

    // Create the Mesh.
    auto mesh_db   = input_db->getDatabase( "Mesh" );
    auto mgrParams = std::make_shared<AMP::Mesh::MeshParameters>( mesh_db );
    mgrParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    std::shared_ptr<AMP::Mesh::Mesh> mesh = AMP::Mesh::Mesh::buildMesh( mgrParams );

    // Put moab mesh filename onto DB
    std::string moabMeshFile = "input.h5m";
    input_db->putScalar( "moabMeshName", moabMeshFile );

    // Build operator params
    AMP::pout << "Building Moab Operator Parameters" << std::endl;
    auto moabParams = std::make_shared<AMP::Operator::MoabBasedOperatorParameters>( input_db );

    // Build operator
    AMP::pout << "Building Moab Operator" << std::endl;
    auto moabOp = std::make_shared<MoabDummyOperator>( moabParams );

    // Call apply
    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    moabOp->apply( nullVec, nullVec, nullVec, 0.0, 0.0 );

    // Create Parameters for Map Operator
    AMP::pout << "Creating map operator" << std::endl;
    input_db->putScalar( "MoabMapVariable", "TEMPERATURE" );
    auto mapParams = std::make_shared<AMP::Operator::MoabMapOperatorParameters>( input_db );
    mapParams->setMoabOperator( moabOp );
    mapParams->setMesh( mesh );

    // Create DOF manager
    size_t DOFsPerNode  = 1;
    int nodalGhostWidth = 1;
    bool split          = true;
    auto nodalDofMap    = AMP::Discretization::simpleDOFManager::create(
        mesh, AMP::Mesh::GeomType::Vertex, nodalGhostWidth, DOFsPerNode, split );
    auto nodalVar = = std::make_shared<AMP::LinearAlgebra::Variable>( "nodalPressure" );
    auto nodalVec   = AMP::LinearAlgebra::createVector( nodalDofMap, nodalVar, true );

    AMP::pout << "Nodal Vector size: " << nodalVec->getGlobalSize() << std::endl;

    // Now create Moab map operator
    AMP::pout << "Creating Node-Based Moab Map Operator" << std::endl;
    input_db->putScalar( "InterpolateToType", "GeomType::Vertex" );
    auto moabNodeMap = std::make_shared<AMP::Operator::MoabMapOperator>( mapParams );

    // Do interpolation
    moabNodeMap->apply( nullVec, nullVec, nodalVec, 0.0, 0.0 );

    // Check to make sure we didn't just get a vector of zeros
    bool nonZero = false;
    for ( auto v : nodalVec ) {
        if ( v != 0.0 )
            nonZero = true;
    }

    if ( nonZero )
        ut->passes( "Nodal vector is not identically zero" );
    else
        ut->failure( "Nodal vector is identically zero" );

    // Now let's see if the interpolated values are what we expected (should be equal to sum of x-
    // and z-coordinates)
    int offset        = 0;
    int numMismatched = 0;

    // loop over all meshes to create the preprocessor database for that mesh
    auto meshIDs = mesh->getBaseMeshIDs();
    for ( size_t meshIndex = 0; meshIndex < meshIDs.size(); meshIndex++ ) {
        // this is an accessor to all the mesh info.
        auto currentMesh = mesh->Subset( meshIDs[meshIndex] );
        if ( currentMesh.get() == NULL )
            continue;

        std::string meshCoords = "Mesh_Coords";
        auto thisMeshCoords    = currentMesh->getPositionVector( meshCoords );

        for ( unsigned int i = 0; i < currentMesh->numLocalElements( AMP::Mesh::GeomType::Vertex );
              ++i ) {
            double val1 = 100.0 * ( thisMeshCoords->getValueByLocalID( 3 * i ) +
                                    thisMeshCoords->getValueByLocalID(
                                        3 * i + 2 ) ); // AMP  coordinates are in meters
            double val2 = nodalVec->getValueByLocalID( offset + i ); // Moab coordinates are in cm

            // Linear interpolation should be 'exact' because we prescribed a linear function
            // Can't use approx_equal here because it fails for small values (it compares relative
            // rather than absolute
            // difference)
            if ( std::abs( val1 - val2 ) > 1.0e-10 ) {
                numMismatched++;

                // If value didn't match, print out it's index and the values
                AMP::pout << "Mismatch at index " << i << ": " << val1 << " " << val2 << std::endl;
            }
        }

        offset += currentMesh->numLocalElements( AMP::Mesh::GeomType::Vertex );
    }

    if ( numMismatched == 0 )
        ut->passes( "Values interpolated correctly" );
    else
        ut->failure( "Values not interpolated correctly" );

        // How about some output?
        // Useful for making sure everything looks right

#ifdef USE_EXT_SILO
    auto siloWriter = AMP::Utilities::Writer::buildWriter( "Silo" );
    siloWriter->registerMesh( mesh );
    siloWriter->registerVector( nodalVec, mesh, AMP::Mesh::GeomType::Vertex, "Temperatures" );
    siloWriter->writeFile( "Moab_Temp", 0 );
#endif

    if ( ut->NumPassGlobal() == 0 )
        ut->failure( "if it doesn't pass, it must have failed." );
}


int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    moabInterface( &ut );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
