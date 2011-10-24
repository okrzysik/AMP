#ifndef included_AMP_SiloIO
#define included_AMP_SiloIO

#ifdef USE_SILO

#include <string.h>
#include "vectors/Vector.h"
#include "MeshVariable.h"
#include "DOFMap.h"
#include <silo.h>
#include <sstream>
#include <map>

namespace AMP { 
namespace Mesh {

  template <typename MANAGER>
  class SiloIO {
    public:
      typedef          MANAGER                         MeshManager;
      typedef typename MeshManager::MeshIterator       MeshIterator;
      typedef typename MeshManager::Adapter            MeshAdapter;
      typedef typename MeshAdapter::ElementIterator    ElementIterator;
      typedef typename MeshAdapter::NodeIterator       NodeIterator;
      typedef typename MeshManager::MeshNameIterator   MeshNameIterator;
      typedef typename MeshAdapter::DataIterator       DataIterator;

    private:
      DBfile    * d_FileHandle;

      std::map<size_t,size_t>     d_NodeRenumber;
      std::vector<MPI_Request>    d_OutstandingRequests;

      void         flushMPIRequests ();


      double  d_Time;
      unsigned int  d_IterCount;
      bool   d_SyncForMasterFile;
      std::string  d_FileName , d_SubdirName;

      void   computeNodeRenumber ( MeshAdapter & );
      std::vector<char *>  d_MeshName;
      std::map<std::string,std::vector<char *> >  d_VecNames;

      void   writeMesh ( MeshAdapter & , const std::string &mesh_name );
      int    writeConnectivityList ( MeshAdapter & , const std::string & );
      void   writeNodalVector ( MeshAdapter & , AMP::LinearAlgebra::Vector::shared_ptr  p , const std::string &mesh_name , const std::string &vec_name );
      void   writeCellVector ( MeshAdapter & , AMP::LinearAlgebra::Vector::shared_ptr  p , const std::string &mesh_name , const std::string &vec_name );
      void   writeVector ( MeshAdapter & , AMP::LinearAlgebra::Vector::shared_ptr  p , const std::string &mesh_name , const std::string &vec_name );

      void   writeMasterFile ( const std::string &fname );

      void   readMesh ( MeshAdapter & , const std::string & mesh_name );
      void   readNodalVector ( MeshAdapter & , AMP::LinearAlgebra::Vector::shared_ptr  p , const std::string &mesh_name , const std::string &vec_name );
      void   readCellVector ( MeshAdapter & , AMP::LinearAlgebra::Vector::shared_ptr  p , const std::string &mesh_name , const std::string &vec_name );
      void   readVector ( MeshAdapter & , AMP::LinearAlgebra::Vector::shared_ptr  p , const std::string &mesh_name , const std::string &vec_name );

    public:

      std::string getExtension () { return "silo"; }
      void  readFile ( MeshManager & , std::string &fname );

      void  writeFile ( MeshManager & , const std::string &fname );
      void  decompWriteFile ( MeshManager & , MeshAdapter& , const std::string &fname , const std::string &mastername );

      void  readFile ( MeshManager & , const std::string &fname );
      void  decompReadFile ( MeshManager & , MeshAdapter& , const std::string &fname , const std::string &mastername );

      void  updateFile ( MeshManager & , std::string &fname , double t , unsigned int );
  };

}
}

#include "SiloIO.tmpl.h"

#else

namespace AMP { 
namespace Mesh {

  template <typename MANAGER>
  class SiloIO {
    public:
      typedef          MANAGER                         MeshManager;
      typedef typename MeshManager::MeshIterator       MeshIterator;
      typedef typename MeshManager::Adapter            MeshAdapter;
      typedef typename MeshAdapter::ElementIterator    ElementIterator;
      typedef typename MeshAdapter::NodeIterator       NodeIterator;
      typedef typename MeshManager::MeshNameIterator   MeshNameIterator;
      typedef typename MeshAdapter::DataIterator       DataIterator;

      void  readFile ( MeshManager & , const std::string & ) 
      { AMP_ERROR( "SILO not built" ); }

      void  writeFile ( MeshManager & , const std::string & )
      { AMP_ERROR( "SILO not built" ); }

      void  decompWriteFile ( MeshManager & , MeshAdapter& , const std::string &fname , const std::string &mastername )
      { AMP_ERROR( "SILO not built" ); }

      void  updateFile ( MeshManager & , std::string & , double t , unsigned int )
      { AMP_ERROR( "SILO not built" ); }
      std::string getExtension () { return "silo"; }
  };

}
}

#endif
#endif


