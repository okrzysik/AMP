#include "MeshManager.h"
#include "../../vectors/test/test_VectorLoops.h"
#include "vectors/CommCollectVector.h"
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/AMP_MPI.h"

class CommCollectGen
{

  public:
    static AMP::Mesh::MeshManager::Adapter::shared_ptr  mesh;

    typedef AMP::LinearAlgebra::CommCollectVector   vector;

    static  AMP::LinearAlgebra::Variable::shared_ptr   getVariable()
    {
      return AMP::LinearAlgebra::Variable::shared_ptr ();
    }

    static  AMP::LinearAlgebra::Vector::shared_ptr     getVector()
    {
      AMP::LinearAlgebra::Vector::shared_ptr t = mesh->createVector ( AMP::LinearAlgebra::Variable::shared_ptr ( new AMP::Mesh::NodalScalarVariable ( "tt" ) ) );
      AMP::AMP_MPI globalComm(AMP_COMM_WORLD);
      return AMP::LinearAlgebra::CommCollectVector::view ( t , globalComm );
    }
};



AMP::Mesh::MeshManager::Adapter::shared_ptr  CommCollectGen::mesh;

using namespace AMP::unit_test;

int main ( int argc , char ** argv )
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

    // Limit scope to ensure proper memory destruction
    {
        AMP::AMP_MPI globalComm(AMP_COMM_WORLD);
        int rank = globalComm.getRank();
        int size = globalComm.getSize();
        int key = (rank < size/2)?0:1;
        AMP::AMP_MPI newComm = globalComm.split(key,rank);
        AMP::Mesh::LibMeshAdapter::initAdapter ( argc , argv , newComm );   // We need to initialize the new communicator
        CommCollectGen::mesh = AMP::Mesh::MeshManager::Adapter::shared_ptr ( new AMP::Mesh::LibMeshAdapter () );
        CommCollectGen::mesh->readExodusIIFile ( "pellet_1x.e" );

        test_managed_vectors_bottom<CommCollectGen> ( &ut );
        CommCollectGen::mesh = AMP::Mesh::MeshManager::Adapter::shared_ptr ();

    } // local scope

    ut.report ();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;

}
