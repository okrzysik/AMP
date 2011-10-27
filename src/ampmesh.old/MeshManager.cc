
#include "MeshManager.h"


namespace AMP { 
namespace Mesh {

  MeshManagerParameters::MeshManagerParameters ( const boost::shared_ptr<AMP::Database> &db )
  {
    d_db = db;
    d_argc = 2;
    d_argv = new char * [2];
    d_argv[0] = new char [4];
    strcpy ( d_argv[0] , "amp" );
    d_argv[1] = new char [12];
    strcpy ( d_argv[1] , "--keep-cout" );
  }

  MeshManagerParameters::~MeshManagerParameters ()
  {
    for ( int i = 0 ; i != d_argc ; i++ )
      delete [] d_argv[i];
    delete [] d_argv;
  }

}
}

