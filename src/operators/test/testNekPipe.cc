//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   operators/test/testNekPipe.cc
 * \brief  This tests the header file that accesses Nek-5000; runs Nek exmaple moab/pipe.
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

// Nek includes
#include "nek/Nek5000_API2.h"

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

void nekPipe(AMP::UnitTest *ut)
{
    // this test is based on testSNES-B-TM-4


    AMP::pout << "Starting to run Nek-500 for the pipe problem"<< std::endl;

#ifdef USE_NEK     
	  std::cout<<"Preparing NEK Init"<< std::endl;
    // must agree with first line of SESSION.NAME
    std::string nekName = "pipe";
    NEK_PREPARE( nekName );
	  std::cout<<"Calling NEK Init"<< std::endl;
    int myMpiComm;
    NEK_INIT( &myMpiComm );
	  std::cout<<"NEK Init succeeded"<<std::endl;
    NEK_SOLVE();
    ut->passes("Nek has created additional problems.");
    
    void *mesh_ptr;
    void *tag;
    getmoabmeshdata_( &mesh_ptr, &tag );

    iMesh_Instance mesh = (iMesh_Instance) mesh_ptr;

    iBase_EntityHandle *ents;
    int ents_alloc = 0, ents_size;
    int ierr;
    iMesh_getEntities(mesh, 0, iBase_REGION,
                      iMesh_ALL_TOPOLOGIES,
                      &ents, &ents_alloc, 
                      &ents_size, &ierr);

    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    int totalSize;
    totalSize = globalComm.sumReduce( ents_size );
    AMP::pout << "Mesh size is " << totalSize << std::endl;

    if( totalSize == 5496 )
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
    std::vector<double> nekMeshCoords;
    result = moabInterface->get_vertex_coordinates( nekMeshCoords );

    AMP::pout << "Retrieved " << nekMeshCoords.size() << " coordinates" << std::endl;
    
    // create MBParallelComm
    AMP::pout << "Creating MBParallelComm" << std::endl;
    int moabCommOut = 0;
    MBParallelComm *moabCommunicator = new MBParallelComm( moabInterface, 
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

    // Create list of points
    int numCoords = 3;
    std::vector<double> myCoords(3*numCoords);

    // First point
    myCoords[0] = 0.0;
    myCoords[1] = 0.0;
    myCoords[2] = 0.0;

    // Second point
    myCoords[3] = 0.1;
    myCoords[4] = 0.2;
    myCoords[5] = 0.3;
    
    // Third point
    myCoords[6] = -0.1;
    myCoords[7] =  0.2;
    myCoords[8] = -0.3;

    // Input coords to coupler
    moabError = moabCoupler.locate_points( &myCoords[0], numCoords );
    AMP::pout<< "we've located the points." << std::endl;

    // Perform interpolation
    MBTag tempTag = reinterpret_cast<MBTag>( tag );
    ///AMP::pout<< "the tag name is "<< tempTag << std::endl;


    std::vector<double> interpTemps(numCoords,0.0);
    moabError = moabCoupler.interpolate( moab::Coupler::LINEAR_FE, tempTag, &interpTemps[0] );

    AMP::pout<< "the temperature is " << interpTemps[0] << std::endl;

    // We are done.
    NEK_END();
    ut->passes("Nek has cleaned itself up.");
#else
    ut->passes("Nek was not used.");
#endif

    if (ut->NumPassGlobal() == 0) ut->failure("if it doesn't pass, it must have failed.");
} 


int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;
    
    try {
        nekPipe(&ut);
	      ut.passes("Nek ran pipe to completion.");
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
//                 end of test_ORIGEN_utils.cc
//---------------------------------------------------------------------------//
