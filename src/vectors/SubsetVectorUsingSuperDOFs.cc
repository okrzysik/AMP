

#include "SubsetVectorUsingSuperDOFs.h"
#include "SubsetVariable.h"



namespace AMP {
namespace LinearAlgebra {

  void  SubsetVectorUsingSuperDOFs::setValuesByLocalID ( int cnt , int *ndx , const double *vals )
  {
    if ( d_ViewVector )
    {
      d_ViewVector->setValuesByLocalID ( cnt , ndx , vals );
    }
    else
    {
      std::runtime_error ( "What?" );
    }
  }

  void  SubsetVectorUsingSuperDOFs::setLocalValuesByGlobalID ( int cnt , int *ndx , const double *vals )
  {
    if ( d_ViewVector )
    {
      d_ViewVector->setLocalValuesByGlobalID ( cnt , ndx , vals );
    }
    else
    {
      for ( int i = 0 ; i != cnt ; i++ )
      {
        d_Space[d_Indexer->getSubID ( ndx[i] )] = vals[ndx[i]];
      }
    }
  }

  void  SubsetVectorUsingSuperDOFs::addValuesByLocalID ( int cnt , int *ndx , const double *vals )
  {
    if ( d_ViewVector )
    {
      d_ViewVector->addValuesByLocalID ( cnt , ndx , vals );
    }
    else
    {
      std::runtime_error ( "What?" );
    }
  }

  void  SubsetVectorUsingSuperDOFs::addLocalValuesByGlobalID ( int cnt , int *ndx , const double *vals )
  {
    if ( d_ViewVector )
    {
      d_ViewVector->addLocalValuesByGlobalID ( cnt , ndx , vals );
    }
    else
    {
      for ( int i = 0 ; i != cnt ; i++ )
      {
        d_Space[d_Indexer->getSubID ( ndx[i] )] += vals[ndx[i]];
      }
    }
  }

  void  SubsetVectorUsingSuperDOFs::getLocalValuesByGlobalID ( int cnt , int *ndx , double *vals ) const
  {
    if ( d_ViewVector )
    {
      d_ViewVector->getLocalValuesByGlobalID ( cnt , ndx , vals );
    }
    else
    {
      for ( int i = 0 ; i != cnt ; i++ )
      {
        vals[ndx[i]] = d_Space[d_Indexer->getSubID ( ndx[i] )];
      }
    }
  }

  Vector::shared_ptr   SubsetVectorUsingSuperDOFs::view ( Vector::shared_ptr  v , Variable::shared_ptr var )
  {
    SubsetVectorUsingSuperDOFs *retVal = new SubsetVectorUsingSuperDOFs();
    AMP_ASSERT ( var->isA<SubsetVariable> () );
    retVal->setVariable ( var );
    retVal->d_ViewVector = v;
    retVal->d_Indexer = var->castTo<SubsetVariable>().getIndexer();
    retVal->d_GlobalSize = v->getGlobalSize();
    retVal->d_LocalSize = v->getLocalSize();
    retVal->d_LocalStartID = v->getCommunicationList()->getStartGID();
    retVal->setCommunicationList ( v->getCommunicationList()->subset ( var->castTo<SubsetVariable>().getIndexer() ) );
    retVal->computeIDMap ();
    return Vector::shared_ptr ( retVal );
  }

}
}

