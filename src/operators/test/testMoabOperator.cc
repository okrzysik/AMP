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

#include "ampmesh/SiloIO.h"

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
                myTemps[i] = myCoords[i];
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

    // Log all nodes
    AMP::PIO::logAllNodes( "output_testMoabOperator" );

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
    moabDB->putInteger("NumberOfMeshes",1);
    boost::shared_ptr<AMP::Database> meshDB = moabDB->putDatabase("Mesh_1");
    meshDB->putString("Filename",ampMeshFile);
    meshDB->putString("MeshName","fuel");
    meshDB->putDouble("x_offset",0.0);
    meshDB->putDouble("y_offset",0.0);
    meshDB->putDouble("z_offset",0.0);

    // Create Mesh Manager
    AMP::pout << "Creating mesh manager" << std::endl;
    typedef AMP::Mesh::MeshManagerParameters    MeshMgrParams;
    typedef boost::shared_ptr< MeshMgrParams >  SP_MeshMgrParams;

    typedef AMP::Mesh::MeshManager              MeshMgr;
    typedef boost::shared_ptr< MeshMgr >        SP_MeshMgr;

    SP_MeshMgrParams mgrParams( new MeshMgrParams( moabDB ) );
    SP_MeshMgr       manager(   new MeshMgr( mgrParams ) );
    

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
    typedef AMP::LinearAlgebra::Variable      AMPVar;
    typedef boost::shared_ptr< AMPVar >       SP_AMPVar;

    typedef AMP::LinearAlgebra::MultiVariable AMPMultiVar;
    typedef boost::shared_ptr< AMPMultiVar >  SP_AMPMultiVar;

    SP_AMPMultiVar allGPPressures( new AMPMultiVar( "AllPressures" ) );
    SP_AMPMultiVar allNodePressures( new AMPMultiVar( "AllPressures" ) );

    MeshMgr::MeshIterator currentMesh;
    for( currentMesh  = manager->beginMeshes();
         currentMesh != manager->endMeshes();
         currentMesh++ )
    {
        // Make variable on this mesh
        SP_AMPVar thisNodeVar( new AMP::Mesh::NodalScalarVariable("MoabNodeTemperature", *currentMesh));

        // Add variable on this mesh to multivariable
        allNodePressures->add( thisNodeVar );
    }

    // Have mesh manager create vector over all meshes
    AMP::LinearAlgebra::Vector::shared_ptr r_node = manager->createVector( allNodePressures );
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
    for( currentMesh  = manager->beginMeshes();
         currentMesh != manager->endMeshes();
         currentMesh++ )
    {
        std::string meshCoords = "Mesh_Coords";
        SP_AMPVec thisMeshCoords = (*currentMesh)->getPositionVector( meshCoords );

        for( int i=0; i<(*currentMesh)->numLocalNodes(); ++i )
        {
            double val1 = 100.0 * thisMeshCoords->getValueByLocalID(3*i); // AMP  coordinates are in meters
            double val2 = r_node->getValueByLocalID(offset+i);            // Moab coordinates are in cm

            // Can't use approx_equal here because it fails for small values
            if( std::abs(val1-val2) > 1.0e-8 )
            {
                numMismatched++;

                // If value didn't match, print out it's index and the values
                AMP::pout << "Index " << i << ": " << 100.0*thisMeshCoords->getValueByLocalID(3*i) << " " << r_node->getValueByLocalID(offset+i) << std::endl;
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
    manager->registerVectorAsData( r_node );
    manager->writeFile<AMP::Mesh::SiloIO>( "Moab_Temp", 0 );
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
