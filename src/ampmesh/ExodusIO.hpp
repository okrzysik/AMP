#ifndef included_AMP_ExodusIO
#define included_AMP_ExodusIO

namespace AMP { 
namespace Mesh {

  template <typename MANAGER>
  void ExodusIO<MANAGER>::createFile ( const char *fname )
  {
    int  num1 , num2;
    num1 = num2 = sizeof ( double );  // Magic constants to ensure doubles are 64 bits..
    ExodusII::ex_create ( fname , EX_CLOBBER , &num1 , &num2 );
    openFile ( fname );
  }

  template <typename MANAGER>
  void ExodusIO<MANAGER>::openFile ( const char *fname )
  {
    int  num1 , num2 , ex_version;
    num1 = num2 = sizeof ( double );  // Magic constants to ensure doubles are 64 bits..

    d_FileHandle = ExodusII::ex_open ( fname , EX_WRITE , & num1 , & num2 , &ex_version );
  }

  template <typename MANAGER>
  void ExodusIO<MANAGER>::writeSingleMesh ( MeshAdapter &mesh , const char *fname )
  {
    openFile ( fname );
  }

}
}

#endif


