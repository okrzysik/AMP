#ifndef included_test_SurfaceSelectorTests
#define included_test_SurfaceSelectorTests

namespace AMP {
namespace unit_test {


template <typename GENERATOR, typename FACTORY>
class  VerifySurfaceSetValue
{
    public:
      static const char * get_test_name () { return "setValue on a Surface"; }

      static  void  run_test ( AMP::UnitTest *utils )
      {
        GENERATOR meshGenerator;
        AMP::Mesh::MeshAdapter::shared_ptr mesh = meshGenerator.getMesh();
        AMP::LinearAlgebra::Vector::shared_ptr  volVec = FACTORY::getVolumeVector();
        AMP::LinearAlgebra::Vector::shared_ptr  surfVec = FACTORY::subsetVolumeVector( volVec );

        surfVec->setToScalar ( 1.0 );

        AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator  b = mesh->beginOwnedBoundary( FACTORY::BoundaryID );
        AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator  e = mesh->endOwnedBoundary( FACTORY::BoundaryID );
        AMP::Mesh::DOFMap::shared_ptr dofMap = mesh->getDOFMap ( volVec->getVariable() );
        unsigned int dofsPerObj = dofMap->numDOFsPerObject();
        while ( b != e )
        {
          for ( unsigned int i = 0 ; i != dofsPerObj ; i++ )
          {
            size_t  curEntry = dofMap->getGlobalID ( b->globalID () , i );
            surfVec->setValueByGlobalID ( curEntry , 1.0 );
          }
          b++;
        }

        b = mesh->beginOwnedBoundary( FACTORY::BoundaryID );
        while ( b != e )
        {
          for ( unsigned int i = 0 ; i != dofsPerObj ; i++ )
          {
            size_t  curEntry = dofMap->getGlobalID ( b->globalID () , i );
            double  curVal = volVec->getValueByGlobalID ( curEntry );
            if ( curVal != 1.0 )
            {
              utils->failure ( "Surface not set correctly" );
              return;
            }
            volVec->setValueByGlobalID ( curEntry , -1.0 );
          }
          b++;
        }

        double volMax = volVec->max();
        if ( volMax > 0.0 )
          utils->failure ( "Volume vec not set correctly" );
        else
          utils->passes ( "Volume vec set correctly" );

        double  l1norm = volVec->L1Norm();
        double  l1guess = (double) surfVec->getGlobalSize();
        if ( fabs ( l1norm - l1guess ) < 0.00001 )
          utils->passes ( "Volume vec has correct L1Norm" );
        else
          utils->failure ( "Volume vec has incorrect L1Norm" );
      }
};


template <typename GENERATOR, typename FACTORY>
class  VerifySurfaceSetScalar
{
    public:
      static const char * get_test_name () { return "setScalar on a Surface"; }

      static  void  run_test ( AMP::UnitTest *utils )
      {
        GENERATOR meshGenerator;
        AMP::Mesh::MeshAdapter::shared_ptr mesh = meshGenerator.getMesh();
        AMP::LinearAlgebra::Vector::shared_ptr  volVec = FACTORY::getVolumeVector();
        AMP::LinearAlgebra::Vector::shared_ptr  surfVec = FACTORY::subsetVolumeVector( volVec );

        surfVec->setToScalar ( 1.0 );

        AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator  b = mesh->beginOwnedBoundary ( FACTORY::BoundaryID );
        AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator  e = mesh->endOwnedBoundary ( FACTORY::BoundaryID );
        AMP::Mesh::DOFMap::shared_ptr dofMap = mesh->getDOFMap ( volVec->getVariable() );
        unsigned int dofsPerObj = dofMap->numDOFsPerObject();
        while ( b != e )
        {
          for ( unsigned int i = 0 ; i != dofsPerObj ; i++ )
          {
            size_t  curEntry = dofMap->getGlobalID ( b->globalID () , i );
            double  curVal = volVec->getValueByGlobalID ( curEntry );
            if ( curVal != 1.0 )
            {
              utils->failure ( "Surface not set correctly" );
              return;
            }
            volVec->setValueByGlobalID ( curEntry , -1.0 );
          }
          b++;
        }

        double volMax = volVec->max();
        if ( volMax > 0.0 )
          utils->failure ( "Volume vec not set correctly" );
        else
          utils->passes ( "Volume vec set correctly" );

        double  l1norm = volVec->L1Norm();
        double  l1guess = (double) surfVec->getGlobalSize();
        if ( fabs ( l1norm - l1guess ) < 0.00001 )
          utils->passes ( "Volume vec has correct L1Norm" );
        else
          utils->failure ( "Volume vec has incorrect L1Norm" );
      }
};


}
}


#endif
