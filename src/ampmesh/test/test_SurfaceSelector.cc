#include "string.h"

#include "vectors/Variable.h"
#include "ampmesh/DOFMap.h"
#include "ampmesh/MeshManager.h"
#include "ampmesh/MeshAdapter.h"
#include "ampmesh/LibMeshAdapter.h"
#include "ampmesh/MeshManager.h"

#include "test_MeshAdapterLoop.h"
#include "test_SurfaceSelectorTests.h"

#include "utils/AMPManager.h"
#include "utils/UnitTest.h"


template <typename GENERATOR, int BID , int SIZE , bool NODAL , bool RUNTIME>
class  SurfaceVectorGenerator
{
public:
    enum
    {
        BoundaryID = BID
    };

    typedef  AMP::LinearAlgebra::Vector               vector;
    typedef  AMP::Mesh::VS_ByMeshIteratorTmpl < AMP::Mesh::MeshManager::Adapter , AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator >                          selector;

    static AMP::LinearAlgebra::Variable::shared_ptr  getVariable ()
    {
        GENERATOR meshGenerator;
        AMP::Mesh::MeshAdapter::shared_ptr mesh = meshGenerator.getMesh();
        if ( NODAL ) {
            return AMP::LinearAlgebra::Variable::shared_ptr ( new AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable,SIZE>( "test vector" , mesh ) );
        } else {
            if ( RUNTIME )
                return AMP::LinearAlgebra::Variable::shared_ptr ( new AMP::Mesh::RunTimeIntegrationPointVariable ( "test vector" , SIZE , mesh ) );
            else
                return AMP::LinearAlgebra::Variable::shared_ptr ( new AMP::LinearAlgebra::VectorVariable<AMP::Mesh::IntegrationPointVariable,SIZE>( "test vector" , mesh ) );
        }
    }

    static  AMP::LinearAlgebra::Vector::shared_ptr  getVolumeVector()
    {
        GENERATOR meshGenerator;
        AMP::Mesh::MeshAdapter::shared_ptr mesh = meshGenerator.getMesh();
        return mesh->createVector( getVariable() );
    }

    static  AMP::LinearAlgebra::Vector::shared_ptr  subsetVolumeVector ( AMP::LinearAlgebra::Vector::shared_ptr v )
    {
        GENERATOR meshGenerator;
        AMP::Mesh::MeshAdapter::shared_ptr mesh = meshGenerator.getMesh();
        AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator  b = mesh->beginOwnedBoundary ( BID );
        AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator  e = mesh->endOwnedBoundary ( BID );
        //AMP_ASSERT ( b != e );
        return v->select ( selector ( mesh , b , e , SIZE ) , "ans" );
    }

    static  AMP::LinearAlgebra::Vector::shared_ptr  getVector()
    {
        return subsetVolumeVector ( getVolumeVector() );
    }

    static  AMP::Mesh::DOFMap::shared_ptr  getDOFMap()
    {
        GENERATOR meshGenerator;
        AMP::Mesh::MeshAdapter::shared_ptr mesh = meshGenerator.getMesh();
        return mesh->getDOFMap ( getVariable() );
    }

};



template <typename GENERATOR, typename FACTORY>
void  surfaceLoops ( AMP::UnitTest *ut )
{
    test_managed_vectors_loop<FACTORY> ( ut );
    VerifySurfaceSetScalar<GENERATOR,FACTORY>::run_test ( ut );
    VerifySurfaceSetValue<GENERATOR,FACTORY>::run_test ( ut );
}

template <typename GENERATOR, int BID , bool NODAL , bool RUNTIME>
void  varySize ( AMP::UnitTest *ut )
{
    // Run tests
    surfaceLoops < GENERATOR, SurfaceVectorGenerator<GENERATOR,BID,1,NODAL,RUNTIME> > ( ut );
    // surfaceLoops < SurfaceVectorGenerator < ExodusReaderGenerator , BID , 3 , NODAL , RUNTIME> > ( ut );
}

template <typename GENERATOR, int BID , bool NODAL>
void  varyRuntime ( AMP::UnitTest *ut )
{
    varySize < GENERATOR, BID , NODAL , false > ( ut );
    // if ( !NODAL )
    //    varySize < GEN , BID , NODAL , true > ( ut );
}

template <typename GENERATOR, int BID>
void  varyNodal ( AMP::UnitTest *ut )
{
    varyRuntime <GENERATOR,BID,true> ( ut );
    // varyRuntime <GEN,BID,false> ( ut );
}

template <typename GENERATOR>
void  varyBID ( AMP::UnitTest *ut )
{
    varyNodal <GENERATOR,1> ( ut );
    // varyNodal <GEN , 2> ( ut );
    // varyNodal <GEN , 4> ( ut );
}



int main ( int argc , char **argv )
{

    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

    varyBID < ExodusReaderGenerator > ( &ut );

    ut.report ();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;

}
