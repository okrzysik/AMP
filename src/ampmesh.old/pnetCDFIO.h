#ifndef included_AMP_ExodusIO
#define included_AMP_ExodusIO

#include <sstream>
#include <map>

//#include <netcdf.h>
#include <pnetcdf.h>
#ifdef USE_MPI
    #include "mpi.h"
#else
    typedef int MPI_Offset
#endif


namespace AMP { 
namespace Mesh {

  template <typename MANAGER>
  class pnetCDFIO {
    private:
      int        d_ncid;
      int        d_ScalarId;
      int        d_3VectorId;

    public:
      typedef          MANAGER                    MeshManager;
      typedef typename MANAGER::MeshIterator      MeshIterator;
      typedef typename MANAGER::Adapter           MeshAdapter;
      typedef typename MeshAdapter::DataIterator  VectorIterator;

      pnetCDFIO () {}

      void  readSingleMesh ( MeshAdapter & ) {}

      void  writeSingleMesh ( MeshAdapter & , std::string );

      void  readFile ( MeshManager & , const char *fname );

      void  writeFile ( MeshManager & , const char *fname );
  };

  template <typename MANAGER>
  void pnetCDFIO<MANAGER>::writeFile ( MeshManager &meshes , const char *fname )
  {
    int idp;
    AMP::AMP_MPI globalComm(AMP_COMM_WORLD);
    d_ncid = ncmpi_create ( globalComm.getCommunicator() , fname , NC_SHARE , MPI_INFO_NULL , &idp );
    MeshIterator  cur_mesh = meshes.beginMeshes();
    ncmpi_redef ( d_ncid );
    ncmpi_def_dim ( d_ncid , "3Vector" , 3 , &d_3VectorId );
    ncmpi_def_dim ( d_ncid , "Scalar" , 1 , &d_ScalarId );
    ncmpi_enddef ( d_ncid );
    while ( cur_mesh != meshes.endMeshes() )
    {
      writeSingleMesh ( **cur_mesh , "pellet" );
      cur_mesh++;
    }
    ncmpi_close ( d_ncid );
  }

  template <typename MANAGER>
  void pnetCDFIO<MANAGER>::writeSingleMesh ( MeshAdapter &mesh , std::string prefix )
  {
    ncmpi_redef ( d_ncid );
    int  nodal_id;
    std::stringstream  nodal_dim_name;
    nodal_dim_name << prefix << "_nodal";
    ncmpi_def_dim ( d_ncid , nodal_dim_name.str().c_str() , mesh.numTotalNodes() , &nodal_id );
    ncmpi_enddef ( d_ncid );
        
    VectorIterator  curVec = mesh.beginData();
    while ( curVec != mesh.endData() )
    {
      int dim_array[2]; 
      MPI_Offset start[2], count[2], stride[2];
      dim_array[0] = nodal_id;
      dim_array[1] = (*curVec)->getVariable()->DOFsPerObject() == 1 ? d_ScalarId : d_3VectorId;
      start[0] = 0;
      start[1] = 0;
      count[0] = mesh.numTotalNodes();
      count[1] = (*curVec)->getVariable()->DOFsPerObject(); 
      stride[0] = 1;
      stride[1] = 1;
      std::stringstream  vec_name;
      vec_name << prefix << "_" << (*curVec)->getVariable()->getName();
      int var_id;
      ncmpi_redef ( d_ncid );
      ncmpi_def_var ( d_ncid , vec_name.str().c_str() , NC_DOUBLE , 2 , dim_array , &var_id );
      ncmpi_enddef ( d_ncid );
      std::cout << "!!!" << ncmpi_strerror ( ncmpi_put_vars_double_all ( d_ncid , var_id , start , count , stride , (*curVec)->template getRawDataBlock<double>() ) ) << std::endl;
      for ( int i = 0 ; i != 100 ; i++ )
        std::cout << (*curVec)->template getRawDataBlock<double>()[i] << " ";
      std::cout << std::endl;
      curVec++;
    }
  }

}
}

#endif


