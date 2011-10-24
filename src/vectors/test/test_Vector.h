#include <vectors/NullVariable.h>
#include <vectors/NativePetscVector.h>
#include "vectors/SimpleVector.h"
#include "test_VectorTests.h"
#include "utils/AMP_MPI.h"

/// \cond UNDOCUMENTED

namespace AMP {
namespace unit_test {


template <typename T>
class  CloneFactory
{
public:
    typedef typename T::vector                  vector;

    static AMP::LinearAlgebra::Variable::shared_ptr  getVariable() {
        return T::getVariable()->cloneVariable ();
    }

    static AMP::LinearAlgebra::Vector::shared_ptr getVector() {
        return T::getVector()->cloneVector();
    }
};


template <int I>
class  SimpleVectorFactory
{
public:
    typedef AMP::LinearAlgebra::SimpleVector                  vector;

    static AMP::LinearAlgebra::Variable::shared_ptr  getVariable() {
        return AMP::LinearAlgebra::Variable::shared_ptr ();
    }

    static AMP::LinearAlgebra::Vector::shared_ptr getVector() {
        return AMP::LinearAlgebra::SimpleVector::create ( I , getVariable() );
    }
};


template <typename T>
class  SimpleManagedVectorFactory
{
  public:
    typedef  T                            vector;

    static AMP::LinearAlgebra::Variable::shared_ptr  getVariable()
    {
      return AMP::LinearAlgebra::Variable::shared_ptr ( new AMP::LinearAlgebra::NullVariable ( "..." ));   // no variable.....
    }

    static AMP::LinearAlgebra::Vector::shared_ptr  getVector()
    {
      const int num_local = 210;
      AMP::AMP_MPI globalComm(AMP_COMM_WORLD);
      AMP::LinearAlgebra::EpetraVectorEngineParameters *p = new AMP::LinearAlgebra::EpetraVectorEngineParameters ( num_local , num_local*globalComm.getSize(), globalComm );
      for ( int i = 0; i != num_local ; i++ )
          p->addMapping ( i , i + globalComm.getRank()*num_local );
      AMP::LinearAlgebra::ManagedVectorParameters *p1 = new AMP::LinearAlgebra::ManagedVectorParameters;
      p1->d_Engine = AMP::LinearAlgebra::VectorEngine::shared_ptr ( new AMP::LinearAlgebra::EpetraVectorEngine ( AMP::LinearAlgebra::VectorEngineParameters::shared_ptr ( p ) , AMP::LinearAlgebra::VectorEngine::BufferPtr ( new AMP::LinearAlgebra::VectorEngine::Buffer ( 210 ) ) ) );
      p1->d_CommList = AMP::LinearAlgebra::CommunicationList::createEmpty ( 210 );
      AMP::LinearAlgebra::Vector::shared_ptr retval ( new T ( AMP::LinearAlgebra::VectorParameters::shared_ptr ( p1 ) ) );
      retval->setVariable ( AMP::LinearAlgebra::Variable::shared_ptr ( new AMP::LinearAlgebra::NullVariable ( "Test Vector" ) ) );
      return retval;
    }
};

template <typename T>
class  PetscManagedVectorFactory
{
  public:
    typedef  T                            vector;

    static AMP::LinearAlgebra::Variable::shared_ptr  getVariable ()
    {
      return AMP::LinearAlgebra::Variable::shared_ptr ();   // no variable.....
    }

    static AMP::LinearAlgebra::Vector::shared_ptr  getVector ()
    {
      Vec v;
      AMP::AMP_MPI globalComm(AMP_COMM_WORLD);
      VecCreate ( globalComm.getCommunicator(), &v );
      VecSetSizes ( v , 15 , PETSC_DECIDE );
      boost::shared_ptr<AMP::LinearAlgebra::NativePetscVectorParameters> npvParams ( new AMP::LinearAlgebra::NativePetscVectorParameters ( v ) );
      npvParams->d_Deleteable = true;
      AMP::LinearAlgebra::NativePetscVector *newVec = new AMP::LinearAlgebra::NativePetscVector ( npvParams );
      VecSetFromOptions ( v );
      newVec->assemble();
      AMP::LinearAlgebra::ManagedVectorParameters *p1 = new AMP::LinearAlgebra::ManagedVectorParameters;
      p1->d_Engine = AMP::LinearAlgebra::VectorEngine::shared_ptr ( newVec );
      p1->d_CommList = AMP::LinearAlgebra::CommunicationList::createEmpty ( 210 );
      AMP::LinearAlgebra::Vector::shared_ptr retval ( new T ( AMP::LinearAlgebra::VectorParameters::shared_ptr ( p1 ) ) );
      retval->setVariable ( AMP::LinearAlgebra::Variable::shared_ptr ( new AMP::LinearAlgebra::NullVariable ( "Test Vector" ) ) );
      return retval;
    }
};

template <typename TYPE1 , typename FACTORY2>
class  ViewFactory
{
  public:
    typedef  TYPE1                           vector;

    static  AMP::LinearAlgebra::Variable::shared_ptr  getVariable()
    {
      return FACTORY2::getVariable();
    }

    static  AMP::LinearAlgebra::Vector::shared_ptr getVector()
    {
      return vector::view ( FACTORY2::getVector() );
    }
};

template <typename FACTORY1,typename FACTORY2>
class  DualVectorFactory
{
  public:
    typedef  AMP::LinearAlgebra::DualVector                vector;

    static AMP::LinearAlgebra::Variable::shared_ptr  getVariable()
    {
      return AMP::LinearAlgebra::Variable::shared_ptr ( new AMP::LinearAlgebra::DualVariable ( FACTORY1::getVariable() , FACTORY2::getVariable() ) );
    }

    static AMP::LinearAlgebra::Vector::shared_ptr getVector ()
    {
      AMP::LinearAlgebra::Vector::shared_ptr p1 ( FACTORY1::getVector() );
      AMP::LinearAlgebra::Vector::shared_ptr p2 ( FACTORY2::getVector() );
      return AMP::LinearAlgebra::DualVector::create ( p1 , p2 , getVariable() );
    }

};

template <typename FACTORY1,int NUM1 ,typename FACTORY2,int NUM2>
class  MultiVectorFactory
{
  public:
    typedef  AMP::LinearAlgebra::MultiVector              vector;

    static AMP::LinearAlgebra::Variable::shared_ptr  getVariable()
    {
      AMP::LinearAlgebra::MultiVariable  *newVar = new AMP::LinearAlgebra::MultiVariable ( "var1" );
      for ( int i = 0 ; i != NUM1 ; i++ )
        newVar->add ( FACTORY1::getVariable() );
      for ( int i = 0 ; i != NUM2 ; i++ )
        newVar->add ( FACTORY2::getVariable() );
      return AMP::LinearAlgebra::Variable::shared_ptr ( newVar );
    }

    static AMP::LinearAlgebra::Vector::shared_ptr  getVector()
    {
      AMP::AMP_MPI globalComm(AMP_COMM_WORLD);
      AMP::LinearAlgebra::Vector::shared_ptr  retVal = AMP::LinearAlgebra::MultiVector::create ( getVariable () , globalComm );
      for ( int i = 0 ; i != NUM1 ; i++ )
        retVal->castTo<AMP::LinearAlgebra::MultiVector>().addVector ( FACTORY1::getVector() );
      for ( int i = 0 ; i != NUM2 ; i++ )
        retVal->castTo<AMP::LinearAlgebra::MultiVector>().addVector ( FACTORY2::getVector() );
      return retVal;
    }
};


}
}

/// \endcond
