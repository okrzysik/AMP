#include "test_Vector.h"
#include "vectors/petsc/NativePetscVector.h"
#include "vectors/petsc/ManagedPetscVector.h"
#include "test_PetscVectorTests.h"
#include "utils/AMP_MPI.h"


/// \cond UNDOCUMENTED

namespace AMP {
namespace unit_test {

template <typename FACTORY>
class  PetscCloneFactory
{
  public:
    typedef  AMP::LinearAlgebra::Vector                                 vector;

    static AMP::LinearAlgebra::Vector::shared_ptr getNativeVector()
    {
      return FACTORY::getNativeVector();
    }

    static void  destroyNativeVector ( AMP::LinearAlgebra::NativePetscVector &rhs )
    { FACTORY::destroyNativeVector ( rhs ); }

    static void  destroyNativeVector ( AMP::LinearAlgebra::Vector::shared_ptr rhs )
    { FACTORY::destroyNativeVector ( rhs ); }

    static AMP::LinearAlgebra::Vector::shared_ptr getManagedVector()
    {
      return FACTORY::getManagedVector()->cloneVector();
    }
};

template <typename FACTORY>
class  PetscViewFactory
{
  public:
    typedef  AMP::LinearAlgebra::Vector                                 vector;

    static AMP::LinearAlgebra::Vector::shared_ptr getNativeVector()
    {
      return FACTORY::getNativeVector();
    }

    static void  destroyNativeVector ( AMP::LinearAlgebra::NativePetscVector &rhs )
    { FACTORY::destroyNativeVector ( rhs ); }

    static void  destroyNativeVector ( AMP::LinearAlgebra::Vector::shared_ptr rhs )
    { FACTORY::destroyNativeVector ( rhs ); }

    static AMP::LinearAlgebra::Vector::shared_ptr getManagedVector()
    {
      return AMP::LinearAlgebra::PetscVector::view ( FACTORY::getManagedVector() );
    }
};

template <typename MANAGED_FACTORY>
class  SimplePetscVectorFactory
{
  public:
    typedef  AMP::LinearAlgebra::Vector                               vector;

    static AMP::LinearAlgebra::Vector::shared_ptr getNativeVector()
    {
      AMP::LinearAlgebra::Vector::shared_ptr  t = getManagedVector();
      Vec  ans;
      AMP::AMP_MPI globalComm(AMP_COMM_WORLD);
      VecCreate ( globalComm.getCommunicator(), &ans );
      VecSetSizes ( ans , t->getLocalSize() , PETSC_DECIDE );
      VecSetFromOptions ( ans );
      boost::shared_ptr<AMP::LinearAlgebra::NativePetscVectorParameters> npvParams ( new AMP::LinearAlgebra::NativePetscVectorParameters ( ans ) );
      npvParams->d_Deleteable = true;
      AMP::LinearAlgebra::Vector::shared_ptr retVal ( new AMP::LinearAlgebra::NativePetscVector ( npvParams ) );
      retVal->setVariable ( AMP::LinearAlgebra::Variable::shared_ptr ( new AMP::LinearAlgebra::Variable ( "petsc vector" ) ) );
      return retVal;
    }

    static void  destroyNativeVector ( AMP::LinearAlgebra::NativePetscVector &rhs )
    { VecDestroy ( rhs.getVec() ); }

    static void  destroyNativeVector ( AMP::LinearAlgebra::Vector::shared_ptr rhs )
    {  destroyNativeVector ( *boost::dynamic_pointer_cast<AMP::LinearAlgebra::NativePetscVector> ( rhs ) ); }

    static AMP::LinearAlgebra::Vector::shared_ptr getManagedVector()
    {
      return MANAGED_FACTORY::getVector();
    }
};

template <typename T>
class  SimplePetscNativeFactory : public SimplePetscVectorFactory<SimpleManagedVectorFactory<AMP::LinearAlgebra::ManagedPetscVector> >
{
  public:
    typedef T                  vector;

    static AMP::LinearAlgebra::Variable::shared_ptr  getVariable()
    {
      return AMP::LinearAlgebra::Variable::shared_ptr ( new AMP::LinearAlgebra::Variable ( "dummy" ) ); //No associated variable
    }

    static AMP::LinearAlgebra::Vector::shared_ptr getVector()
    {
      return getNativeVector();
    }
};



}
}

/// \endcond
