
#include "DualVector.h"
#include <stdexcept>
#include <algorithm>
#include <math.h>


 /// \cond UNDOCUMENTED
namespace AMP {
namespace LinearAlgebra {

  DualVector::DualVector ( Vector::shared_ptr v1 , Vector::shared_ptr v2 , Variable::shared_ptr names )
    : d_pVector1 ( v1 )
    , d_pVector2 ( v2 )
  {
    setVariable ( names );
    if ( v1->isA<DataChangeFirer> () )
      v1->castTo<DataChangeFirer>().registerListener ( this );
    if ( v2->isA<DataChangeFirer> () )
      v2->castTo<DataChangeFirer>().registerListener ( this );
  }

  Vector::shared_ptr  DualVector::subsetVectorForVariable ( const Variable::shared_ptr  &name )
  {
    Vector::shared_ptr retVal;

    if ( *name == *getVariable() ) retVal = Vector::subsetVectorForVariable ( name );
    else if ( d_pVariable->castTo<DualVariable>().first() == name ) retVal = d_pVector1;
    else if ( d_pVariable->castTo<DualVariable>().second() == name ) retVal = d_pVector2;
    else 
    {
      if ( d_pVector1->isA<DualVector>() ) retVal = d_pVector1->castTo<DualVector>().subsetVectorForVariable ( name );
      if ( !retVal && d_pVector2->isA<DualVector>() ) retVal = d_pVector1->castTo<DualVector>().subsetVectorForVariable ( name );
    }

    return retVal;
  }

  void  DualVector::partitionValues ( size_t num , size_t *indices , const double *vals , size_t start2 ,
                    std::vector<size_t> &Ndx1 , std::vector<size_t> &Ndx2 , 
                    std::vector<double> &Vals1 , std::vector<double> &Vals2 ) const
  {
    size_t num1Vals = 0;
    for (size_t i=0; i<num; i++)
    {
      if ( indices[i] < start2 )
        num1Vals++;
    }
    Ndx1.resize ( num1Vals );
    Ndx2.resize ( num - num1Vals );
    Vals1.resize ( num1Vals );
    Vals2.resize ( num - num1Vals );
    int num1Off = 0;
    int num2Off = 0;
    for (size_t i=0; i<num; i++)
    {
      if ( indices[i] < start2 )
      {
        Ndx1[num1Off] = indices[i];
        Vals1[num1Off] = vals[i];
        num1Off++;
      }
      else
      {
        Ndx2[num2Off] = indices[i] - start2;
        Vals2[num2Off] = vals[i];
        num2Off++;
      }
    }
  }

  void DualVector::setValuesByLocalID ( int num , size_t *indices , const double *vals )
  {
    std::vector<size_t>  ndx1 , ndx2;
    std::vector<double>  val1 , val2;
    partitionValues ( num , indices , vals , d_pVector1->getLocalSize() , ndx1 , ndx2 , val1 , val2 );
    if ( ndx1.size() )
      d_pVector1->setValuesByLocalID ( ndx1.size() , &(ndx1[0]) , &(val1[0]) );
    if ( ndx2.size() )
      d_pVector2->setValuesByLocalID ( ndx2.size() , &(ndx2[0]) , &(val2[0]) );
  }

  void DualVector::setLocalValuesByGlobalID ( int num , size_t *indices , const double *vals )
  {
    std::vector<size_t>  ndx1 , ndx2;
    std::vector<double>  val1 , val2;
    partitionValues ( num , indices , vals , d_pVector1->getGlobalSize() , ndx1 , ndx2 , val1 , val2 );
    if ( ndx1.size() )
      d_pVector1->setLocalValuesByGlobalID ( ndx1.size() , &(ndx1[0]) , &(val1[0]) );
    if ( ndx2.size() )
      d_pVector2->setLocalValuesByGlobalID ( ndx2.size() , &(ndx2[0]) , &(val2[0]) );
  }

  void DualVector::setValuesByGlobalID ( int num , size_t *indices , const double *vals )
  {
    std::vector<size_t>  ndx1 , ndx2;
    std::vector<double>  val1 , val2;
    partitionValues ( num , indices , vals , d_pVector1->getGlobalSize() , ndx1 , ndx2 , val1 , val2 );
    if ( ndx1.size() )
      d_pVector1->setValuesByGlobalID ( ndx1.size() , &(ndx1[0]) , &(val1[0]) );
    if ( ndx2.size() )
      d_pVector2->setValuesByGlobalID ( ndx2.size() , &(ndx2[0]) , &(val2[0]) );
  }

  void DualVector::addValuesByLocalID ( int num , size_t *indices , const double *vals )
  {
    std::vector<size_t>  ndx1 , ndx2;
    std::vector<double>  val1 , val2;
    partitionValues ( num , indices , vals , d_pVector1->getLocalSize() , ndx1 , ndx2 , val1 , val2 );
    if ( ndx1.size() )
      d_pVector1->addValuesByLocalID ( ndx1.size() , &(ndx1[0]) , &(val1[0]) );
    if ( ndx2.size() )
      d_pVector2->addValuesByLocalID ( ndx2.size() , &(ndx2[0]) , &(val2[0]) );
  }

  void DualVector::addLocalValuesByGlobalID ( int num , size_t *indices , const double *vals )
  {
    std::vector<size_t>  ndx1 , ndx2;
    std::vector<double>  val1 , val2;
    partitionValues ( num , indices , vals , d_pVector1->getGlobalSize() , ndx1 , ndx2 , val1 , val2 );
    if ( ndx1.size() )
      d_pVector1->addLocalValuesByGlobalID ( ndx1.size() , &(ndx1[0]) , &(val1[0]) );
    if ( ndx2.size() )
      d_pVector2->addLocalValuesByGlobalID ( ndx2.size() , &(ndx2[0]) , &(val2[0]) );
  }

  void DualVector::addValuesByGlobalID ( int num , size_t *indices , const double *vals )
  {
    std::vector<size_t>  ndx1 , ndx2;
    std::vector<double>  val1 , val2;
    partitionValues ( num , indices , vals , d_pVector1->getGlobalSize() , ndx1 , ndx2 , val1 , val2 );
    if ( ndx1.size() )
      d_pVector1->addValuesByGlobalID ( ndx1.size() , &(ndx1[0]) , &(val1[0]) );
    if ( ndx2.size() )
      d_pVector2->addValuesByGlobalID ( ndx2.size() , &(ndx2[0]) , &(val2[0]) );
  }

  void DualVector::getValuesByGlobalID ( int num , size_t *indices , double *vals ) const
  {
    std::vector<size_t>  ndx1 , ndx2;
    std::vector<double>  val1 , val2;
    partitionValues ( num , indices , vals , d_pVector1->getGlobalSize() , ndx1 , ndx2 , val1 , val2 );

    if ( ndx1.size() )
      d_pVector1->getValuesByGlobalID ( ndx1.size() , &(ndx1[0]) , &(val1[0]) );
    if ( ndx2.size() )
      d_pVector2->getValuesByGlobalID ( ndx2.size() , &(ndx2[0]) , &(val2[0]) );

    std::vector<size_t>::iterator  ndx1_iter , ndx2_iter;
    std::vector<double>::iterator  val1_iter , val2_iter;

    if ( ndx1.size() == 0 )
      ndx1.push_back ( -1 );

    ndx1_iter = ndx1.begin();
    val1_iter = val1.begin();
    val2_iter = val2.begin();

    for ( int i = 0 ; i != num ; i++ )
    {
      if ( indices[i] == *ndx1_iter )
      {
        vals[i] = *val1_iter;
        val1_iter++;
        ndx1_iter++;
      }
      else
      {
        vals[i] = *val2_iter;
        val2_iter++;
      }
    }
  }

  void DualVector::getLocalValuesByGlobalID ( int num , size_t *indices , double *vals ) const
  {
    std::vector<size_t>  ndx1 , ndx2;
    std::vector<double>  val1 , val2;
    partitionValues ( num , indices , vals , d_pVector1->getGlobalSize() , ndx1 , ndx2 , val1 , val2 );

    if ( ndx1.size() )
      d_pVector1->getLocalValuesByGlobalID ( ndx1.size() , &(ndx1[0]) , &(val1[0]) );
    if ( ndx2.size() )
      d_pVector2->getLocalValuesByGlobalID ( ndx2.size() , &(ndx2[0]) , &(val2[0]) );

    std::vector<size_t>::iterator  ndx1_iter , ndx2_iter;
    std::vector<double>::iterator  val1_iter , val2_iter;

    if ( ndx1.size() == 0 )
      ndx1.push_back ( -1 );

    ndx1_iter = ndx1.begin();
    val1_iter = val1.begin();
    val2_iter = val2.begin();

    for ( int i = 0 ; i != num ; i++ )
    {
      if ( indices[i] == *ndx1_iter )
      {
        vals[i] = *val1_iter;
        val1_iter++;
        ndx1_iter++;
      }
      else
      {
        vals[i] = *val2_iter;
        val2_iter++;
      }
    }
  }
  void DualVector::getValuesByLocalID ( int num , size_t *indices , double *vals ) const
  {
    std::vector<size_t>  ndx1 , ndx2;
    std::vector<double>  val1 , val2;
    partitionValues ( num , indices , vals , d_pVector1->getLocalSize() , ndx1 , ndx2 , val1 , val2 );

    if ( ndx1.size() )
      d_pVector1->getValuesByLocalID ( ndx1.size() , &(ndx1[0]) , &(val1[0]) );
    if ( ndx2.size() )
      d_pVector2->getValuesByLocalID ( ndx2.size() , &(ndx2[0]) , &(val2[0]) );

    std::vector<size_t>::iterator  ndx1_iter , ndx2_iter;
    std::vector<double>::iterator  val1_iter , val2_iter;

    if ( ndx1.size() == 0 )
      ndx1.push_back ( -1 );

    ndx1_iter = ndx1.begin();
    val1_iter = val1.begin();
    val2_iter = val2.begin();

    for ( int i = 0 ; i != num ; i++ )
    {
      if ( indices[i] == *ndx1_iter )
      {
        vals[i] = *val1_iter;
        val1_iter++;
        ndx1_iter++;
      }
      else
      {
        vals[i] = *val2_iter;
        val2_iter++;
      }
    }
  }

}
}

/// \endcond
