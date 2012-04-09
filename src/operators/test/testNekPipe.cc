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
#include "nek/Nek5000_API2.h"

#include "iMesh.h"

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
#ifdef USE_NEK     
    //NEK_END();
#endif
    return num_failed;
}

//---------------------------------------------------------------------------//
//                 end of test_ORIGEN_utils.cc
//---------------------------------------------------------------------------//
