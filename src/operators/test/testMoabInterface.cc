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

    // Load mesh from file
    AMP::pout << "Loading mesh file" << std::endl;
    std::string meshFile = "pellet_1x.e";
    iMesh_load(mesh, 0, meshFile.c_str(),options.c_str(),&ierr,meshFile.length(),options.length());

    // Add a temperature tag
    //iBase_TagHandle tempTagHandle;
    //std::string tempTagName = "Temperature";
    //iMesh_createTag(mesh,tempTagName.c_str(),tempTagHandle,1,MB_TYPE_DOUBLE,&tempTagHandle,&ierr,tempTagName.length());


    //iMesh_setEntSetDblData( 

    iBase_EntityHandle *ents;
    int ents_alloc = 0, ents_size;
    iMesh_getEntities(mesh, 0, iBase_REGION,
                      iMesh_ALL_TOPOLOGIES,
                      &ents, &ents_alloc, 
                      &ents_size, &ierr);

    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    AMP::pout << "Mesh size is " << ents_size << std::endl;

    // Right now each processor is reading the mesh independently 
    //  so the mesh size is the same as the global size on all procs
    if( ents_size == 3705 )
        ut->passes("Mesh is the right size");
    else
        ut->failure("Mesh is not the right size");

    // accessing MBInterface
    AMP::pout << "Casting to MBInterface" << std::endl;


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
    MBParallelComm *moabCommunicator = new MBParallelComm( moabInterface, 
                                                           globalComm.getCommunicator(),
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

    // Perform interpolation
    MBTag tempTag;
    std::vector<double> interpTemps(numCoords,0.0);
    moabError = moabCoupler.interpolate( moab::Coupler::LINEAR_FE, tempTag, &interpTemps[0] );

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
