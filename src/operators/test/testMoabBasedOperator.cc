//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   operators/test/testMoabOperator.cc
 * \brief  This tests the Moab iMesh interface with MoabBasedOperator
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>

#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "utils/PIO.h"
#include "utils/Database.h"

#include "ampmesh/SiloIO.h"
#include "discretization/simpleDOF_Manager.h"
#include "vectors/VectorBuilder.h"

// AMP Moab Includes
#include "operators/moab/MoabMapOperator.h"
#include "operators/moab/MoabBasedOperator.h"

// MOAB Includes
#include "moab/Interface.hpp"
#include "MBParallelComm.hpp"
#include "MBRange.hpp"
#include "Coupler.hpp"
#include "iMesh.h"
#include "MBiMesh.hpp"

//---------------------------------------------------------------------------//
// Helper Class
//---------------------------------------------------------------------------//

typedef AMP::Operator::MoabBasedOperator            MoabBasedOp;
typedef boost::shared_ptr< MoabBasedOp >            SP_MoabBasedOp;

typedef AMP::Operator::MoabBasedOperatorParameters  MoabOpParams;
typedef boost::shared_ptr<MoabOpParams>             SP_MoabOpParams;

typedef AMP::LinearAlgebra::Vector                  AMPVec;
typedef AMP::LinearAlgebra::Vector::shared_ptr      SP_AMPVec;


class MoabDummyOperator : public MoabBasedOp
{
    public : 

        MoabDummyOperator( SP_MoabOpParams &moabParams )
            : MoabBasedOp( moabParams )
        {
            // Create iMesh instance
            iMesh_Instance mbMesh;
            std::string options;
            int ierr;
            iMesh_newMesh( options.c_str(), &mbMesh, &ierr, options.length() );

            // Convert iMesh to MBiMesh
            MBiMesh *mbimesh = reinterpret_cast<MBiMesh *> (mbMesh);

            // Get Moab interface from MBiMesh
            d_moabInterface = mbimesh->mbImpl;

            // Get root set
            iBase_EntitySetHandle root;
            iMesh_createEntSet( mbMesh, 0, &root, &ierr );

            // Build ParallelComm and get index
            moab::ParallelComm *mbComm = new moab::ParallelComm( d_moabInterface );
            int index = mbComm->get_id();

            // Set read options and load file
            std::string newReadOpts;
            std::ostringstream extraOpt;
            extraOpt << ";PARALLEL_COMM=" << index;
            std::string readOpts = "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARTITION_DISTRIBUTE;PARALLEL_RESOLVE_SHARED_ENTS;PARALLEL_GHOSTS=3.0.1;CPUTIME";
            newReadOpts = readOpts + extraOpt.str();
            std::string filename = moabParams->d_db->getString("moabMeshName");
            AMP::plog << "Moab mesh name is " << filename << std::endl;
            moab::ErrorCode result = d_moabInterface->load_file( filename.c_str(), (EntityHandle *) &root, newReadOpts.c_str() );

            AMP_INSIST( result == MB_SUCCESS, "File not loaded correctly" );

            // Extract nodes
            iBase_EntityHandle *nodes;
            int nodes_alloc = 0;
            int nodes_size;
            iMesh_getEntities( mbMesh, root, iBase_VERTEX, iMesh_ALL_TOPOLOGIES, &nodes, &nodes_alloc, &nodes_size, &ierr );

            // Get list of node coordinates
            std::vector<double> myCoords;
            d_moabInterface->get_vertex_coordinates( myCoords );
            int num_nodes = myCoords.size()/3;

            AMP_INSIST( num_nodes == nodes_size, "Number of nodes must match number of vertex coordinates" );

            AMP::pout << "Found " << num_nodes << " nodes" << std::endl;

            // Put data on mesh
            std::vector<double> myTemps(num_nodes,-1.0);
            for( int i=0; i<num_nodes; ++i )
            {
                // Set temperature to T = x + z
                myTemps[i] = myCoords[i] + myCoords[i + 2*num_nodes];
            }

            // Add temperature tag
            iBase_TagHandle tempTagHandle;
            std::string tempTagName = "TEMPERATURE";
            iMesh_createTag(mbMesh,tempTagName.c_str(),1,MB_TYPE_DOUBLE,&tempTagHandle,&ierr,tempTagName.length());

            // Assign data to tag
            iMesh_setDblArrData(mbMesh,nodes,nodes_size,tempTagHandle,&myTemps[0],myTemps.size(),&ierr);

        }

        void apply( const SP_AMPVec &f, const SP_AMPVec &u, SP_AMPVec &r, double a, double b )
        {
            /* Don't need an apply for this operator */
        }

        void finalize() { /* ... */ };

};


//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void moabInterface(AMP::UnitTest *ut)
{

    // Print out AMP banner
    AMP::Utilities::printBanner();

    // Log all nodes
    AMP::PIO::logAllNodes( "output_testMoabBasedOperator" );

    // Build new database
    AMP::pout << "Building Input Database" << std::endl;
    boost::shared_ptr< AMP::InputDatabase > moabDB( new AMP::InputDatabase("Moab_DB") );

    std::string ampMeshFile = "pellet_1x.e";
    std::string moabMeshFile = "input.h5m";
    moabDB->putString("moabMeshName",moabMeshFile);

    // Build operator params
    typedef AMP::Operator::MoabBasedOperatorParameters MoabOpParams;
    typedef boost::shared_ptr< MoabOpParams >          SP_MoabOpParams;

    AMP::pout << "Building Moab Operator Parameters" << std::endl;
    SP_MoabOpParams moabParams( new MoabOpParams( moabDB ) );

    // Build operator
    typedef AMP::Operator::MoabBasedOperator MoabBasedOp;
    typedef boost::shared_ptr< MoabBasedOp > SP_MoabBasedOp;

    AMP::pout << "Building Moab Operator" << std::endl;
    SP_MoabBasedOp moabOp( new MoabDummyOperator( moabParams ) );

    // Call apply
    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    moabOp->apply( nullVec, nullVec, nullVec, 0.0, 0.0 );

    // Read AMP pellet mesh from file
    moabDB->putInteger("NumberOfMeshes",2);


    AMP::AMP_MPI globalComm(AMP_COMM_WORLD);
    if( globalComm.getSize() == 1 )
        moabDB->putInteger("DomainDecomposition",0);
    else
        moabDB->putInteger("DomainDecomposition",1);

    // Create the multimesh database
    boost::shared_ptr<AMP::Database> meshDB = moabDB->putDatabase("Mesh");
    meshDB->putString("MeshName","PelletMeshes");
    meshDB->putString("MeshType","Multimesh");
    meshDB->putString("MeshDatabasePrefix","Mesh_");
    meshDB->putString("MeshArrayDatabasePrefix","MeshArray_");
    // Create the mesh array database
    boost::shared_ptr<AMP::Database> meshArrayDatabase = meshDB->putDatabase("MeshArray_1");
    int N_meshes=2;
    meshArrayDatabase->putInteger("N",N_meshes);
    meshArrayDatabase->putString("iterator","%i");
    std::vector<int> indexArray(N_meshes);
    for (int i=0; i<N_meshes; i++)
            indexArray[i] = i+1;
    meshArrayDatabase->putIntegerArray("indicies",indexArray);
    meshArrayDatabase->putString("MeshName","pellet_%i");
    meshArrayDatabase->putString("FileName","pellet_1x.e");
    meshArrayDatabase->putString("MeshType","libMesh");
    meshArrayDatabase->putInteger("dim",3);
    meshArrayDatabase->putDouble("x_offset",0.0);
    meshArrayDatabase->putDouble("y_offset",0.0);
    std::vector<double> offsetArray(N_meshes);
    for (int i=0; i<N_meshes; i++)
            offsetArray[i] = ((double) i)*0.0105;
    meshArrayDatabase->putDoubleArray("z_offset",offsetArray);
    meshArrayDatabase->putInteger("NumberOfElements",300);

    // Create Mesh Manager
    AMP::pout << "Creating mesh manager" << std::endl;
    typedef AMP::Mesh::MeshParameters           MeshMgrParams;
    typedef boost::shared_ptr< MeshMgrParams >  SP_MeshMgrParams;

    typedef AMP::Mesh::Mesh                     MeshMgr;
    typedef AMP::Mesh::Mesh::shared_ptr         SP_MeshMgr;

    SP_MeshMgrParams mgrParams( new MeshMgrParams( meshDB ) );
    mgrParams->setComm( AMP::AMP_MPI(AMP_COMM_WORLD) );
    SP_MeshMgr manager = AMP::Mesh::Mesh::buildMesh( mgrParams );
    
    // Create Parameters for Map Operator
    AMP::pout << "Creating map operator" << std::endl;
    typedef AMP::Operator::MoabMapOperatorParameters    MoabMapParams;
    typedef boost::shared_ptr< MoabMapParams >          SP_MoabMapParams;

    typedef AMP::Operator::MoabMapOperator              MoabMap;
    typedef boost::shared_ptr< MoabMap>                 SP_MoabMap;

    moabDB->putString("MoabMapVariable","TEMPERATURE");
    SP_MoabMapParams mapParams( new MoabMapParams( moabDB ) );
    mapParams->setMoabOperator( moabOp );
    mapParams->setMeshManager( manager );

    // Create variable to hold pressure data
    typedef AMP::LinearAlgebra::MultiVariable AMPMultiVar;
    typedef boost::shared_ptr< AMPMultiVar >  SP_AMPMultiVar;

    SP_AMPMultiVar allGPPressures( new AMPMultiVar( "AllPressures" ) );
    SP_AMPMultiVar allNodePressures( new AMPMultiVar( "AllPressures" ) );

    // Create DOF manager
    size_t DOFsPerNode = 1;
    int nodalGhostWidth = 0;
    bool split = true;
    AMP::Discretization::DOFManager::shared_ptr nodalDofMap = AMP::Discretization::simpleDOFManager::create(manager, AMP::Mesh::Vertex, nodalGhostWidth, DOFsPerNode, split);

    // Have mesh manager create vector over all meshes
    AMP::LinearAlgebra::Vector::shared_ptr r_node = AMP::LinearAlgebra::createVector( nodalDofMap, allNodePressures );
    AMP::pout << "Nodal MultiVector size: " << r_node->getGlobalSize() << std::endl; 


    // Now create Moab map operator
    AMP::pout << "Creating Node-Based Moab Map Operator" << std::endl;
    moabDB->putString("InterpolateToType","Vertex");
    SP_MoabMap moabNodeMap( new MoabMap( mapParams ) );

    // Do interpolation
    moabNodeMap->apply( nullVec, nullVec, r_node, 0.0, 0.0 );

    // Check to make sure we didn't just get a vector of zeros
    AMPVec::iterator myIter;
    int ctr=0;
    bool nonZero = false;
    for( myIter  = r_node->begin();
         myIter != r_node->end();
         myIter++ )
    {
        if( *myIter != 0.0 )
            nonZero = true;

        ctr++;
    }

    if( nonZero )
        ut->passes("Nodal vector is not identically zero");
    else
        ut->failure("Nodal vector is identically zero");

    // Now let's see if the interpolated values are what we expected (should be equal to x-coordinate)
    int offset = 0;
    int numMismatched=0;

    // loop over all meshes to create the preprocessor database for that mesh
    std::vector<AMP::Mesh::MeshID> meshIDs = manager->getBaseMeshIDs();
    
    for( size_t meshIndex=0; meshIndex<meshIDs.size(); meshIndex++ )
    {
        // this is an accessor to all the mesh info.
        AMP::Mesh::Mesh::shared_ptr currentMesh = manager->Subset( meshIDs[meshIndex] );
        if( currentMesh.get() == NULL ) continue;

        std::string meshCoords = "Mesh_Coords";
        SP_AMPVec thisMeshCoords = (*currentMesh)->getPositionVector( meshCoords );

        for( unsigned int i=0; i<(*currentMesh)->numLocalNodes(); ++i )
        {
            double val1 = 100.0 * ( thisMeshCoords->getValueByLocalID(3*i)
                                  + thisMeshCoords->getValueByLocalID(3*i+2) ); // AMP  coordinates are in meters
            double val2 = r_node->getValueByLocalID(offset+i);                  // Moab coordinates are in cm

            // Linear interpolation should be 'exact' because we prescribed a linear function
            // Can't use approx_equal here because it fails for small values
            if( std::abs(val1-val2) > 1.0e-10 )
            {
                numMismatched++;

                // If value didn't match, print out it's index and the values
                AMP::pout << "Mismatch at index " << i << ": " << val1 << " " << val2 << std::endl;
            }
        }

        offset += (*currentMesh)->numLocalNodes();
    }

    if( numMismatched == 0 )
        ut->passes("Values interpolated correctly");
    else
        ut->failure("Values not interpolated correctly");

    // How about some output?
    // Useful for making sure everything looks right
    
#ifdef USE_SILO
    AMP::Mesh::SiloIO::shared_ptr  siloWriter( new AMP::Mesh::SiloIO);
    siloWriter->registerVector( r_node, manager, AMP::Mesh::Vertex, "Temperatures" );
    siloWriter->writeFile( "Moab_Temp", 0 );
#endif

    if (ut->NumPassGlobal() == 0) ut->failure("if it doesn't pass, it must have failed.");
} 


int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;
    
    try {
        moabInterface(&ut);
	      ut.passes("Moab interface used correctly.");
    } catch (std::exception &err) {
        std::cout << "ERROR: While testing "<<argv[0] << err.what() << std::endl;
        ut.failure("ERROR: While testing");
    } catch( ... ) {
        std::cout << "ERROR: While testing "<<argv[0] << "An unknown exception was thrown." << std::endl;
        ut.failure("ERROR: While testing");
    }
   
    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}

//---------------------------------------------------------------------------//
//                 end of testMoabOperator.cc
//---------------------------------------------------------------------------//
