#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "utils/PIO.h"

#include "ampmesh/MeshManager.h"

#include "test_MeshGenerators.h"
#include "test_MeshAdapterLoop.h"
#include "test_MeshManagerTests.h"




char **fix_libmesh_cmdline_cout ( int argc , char **argv )
{
  char ** retVal = new char * [ argc + 1 ];
  for ( int i = 0 ; i != argc ; i++ )
  {
    retVal[i] = new char [ strlen ( argv[i] ) + 1 ];
    strcpy ( retVal[i] , argv[i] );
  }
  retVal[argc] = new char [12];
  strcpy ( retVal[argc] , "--keep-cout" );
  return retVal;
}


template <typename ADAPTER>
void  mesh_manager_test ( const std::string &inputName , ParallelUnitTest &ut )
{
  boost::shared_ptr<AMP::InputDatabase>  input_db ( new AMP::InputDatabase ( "input_db" ) );
  AMP::InputManager::getManager()->parseInputFile ( inputName , input_db );
  input_db->printClassData (AMP::plog);

  AMP::Mesh::MeshManagerParameters::shared_ptr  meshmgrParams ( new AMP::Mesh::MeshManagerParameters ( input_db ) );
  typename AMP::Mesh::MeshManagerTmpl<ADAPTER>::shared_ptr  manager ( new AMP::Mesh::MeshManagerTmpl<ADAPTER> ( meshmgrParams ) );

  AdapterFromManagerTmpl<ADAPTER>::d_Manager = manager;
  for ( size_t i = 1 ; i <= manager->numMeshes() ; i++ )
  {
    AdapterFromManagerTmpl<ADAPTER>::d_which = i;
    MeshAdapterLoop <AdapterFromManagerTmpl<ADAPTER> >::test_mesh ( ut );
  }

  AdapterFromManagerTmpl<ADAPTER>::d_Manager.reset();
}


int main ( int argc , char **argv )
{

    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;
    AMP::PIO::logOnlyNodeZero ( "output_MeshManager" );
    mesh_manager_test<AMP::Mesh::LibMeshAdapter> ( "input_MeshManager-1" , ut );
    ut.report();
    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}

