//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   operators/test/testMoabInterface.cc
 * \brief  This tests the Moab iMesh interface
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

// MOAB Includes
#include "moab/Interface.hpp"
#include "MBParallelComm.hpp"
#include "MBRange.hpp"
#include "Coupler.hpp"
#include "iMesh.h"
#include "MBiMesh.hpp"


extern "C" {
    void getmoabmeshdata_( void **, void ** );
}
//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void moabInterface(AMP::UnitTest *ut)
{

#ifdef USE_MOAB
    
    // Create new iMesh instance
    AMP::pout << "Creating iMesh instance" << std::endl;
    iMesh_Instance mesh;
    int ierr;
    std::string options;
    iMesh_newMesh(options.c_str(), &mesh, &ierr, 0);

    // Get Root Set
    AMP::pout << "Getting root set" << std::endl;
    iBase_EntitySetHandle rootSet;
    iMesh_getRootSet( mesh, &rootSet, &ierr );

    // Load mesh from file
    AMP::pout << "Loading mesh file" << std::endl;
    std::string meshFile = "pellet_1x.e";
    iMesh_load(mesh, rootSet, meshFile.c_str(),options.c_str(),&ierr,meshFile.length(),options.length());

    iBase_EntityHandle *ents;
    int ents_alloc = 0, ents_size;
    iMesh_getEntities(mesh, 0, iBase_REGION,
                      iMesh_ALL_TOPOLOGIES,
                      &ents, &ents_alloc, 
                      &ents_size, &ierr);

    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    int totalSize;
    totalSize = globalComm.sumReduce( ents_size );
    AMP::pout << "Mesh size is " << totalSize << std::endl;

    if( totalSize == 3705 )
        ut->passes("Mesh is the right size");
    else
        ut->failure("Mesh is not the right size");

    // accessing MBInterface
    AMP::pout << "Casting to MBInterface" << std::endl;
    //moab::Interface * moabInterface = reinterpret_cast<moab::Interface *> (mesh);
    //

    iMesh_Instance instance = mesh;



    MBiMesh *mbimesh = reinterpret_cast<MBiMesh*>(mesh);
    moab::Interface *moabInterface = mbimesh->mbImpl;

    AMP::pout << "Getting dimension" << std::endl;
    int moabDim;
    MBErrorCode result = moabInterface->get_dimension( moabDim );
    AMP::pout << "Dimension is " << moabDim << std::endl;

    AMP::pout << "Getting vertex coordinates" << std::endl;
    std::vector<double> moabMeshCoords;
    result = moabInterface->get_vertex_coordinates( moabMeshCoords );


    AMP::pout << "Retrieved " << moabMeshCoords.size() << " coordinates" << std::endl;
    


    // create MBParallelComm
    AMP::pout << "Creating MBParallelComm" << std::endl;
    int moabCommOut = 0;
/*    MBParallelComm *moabCommunicator = new MBParallelComm( moabInterface, 
                                                           myMpiComm, 
                                                           &moabCommOut );

    // Access the Range on the source mesh.
    AMP::pout << "Getting MGRange object" << std::endl;
    moab::Range moabRange;
    int problemDimension = 3;
    moab::ErrorCode moabError = moabCommunicator->get_part_entities( moabRange, problemDimension);

    // create MBCoupler 
    AMP::pout << "Creating MBCoupler" << std::endl;
    int moabCouplerID = 0;
    moab::Coupler moabCoupler( moabInterface, 
                               moabCommunicator, 
                               moabRange,
                               moabCouplerID  );

                               */

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
#ifdef USE_NEK     
    //NEK_END();
#endif
    return num_failed;
}

//---------------------------------------------------------------------------//
//                 end of test_ORIGEN_utils.cc
//---------------------------------------------------------------------------//
