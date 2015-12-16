#ifndef included_test_PetscVector
#define included_test_PetscVector

#include "test_PetscVectorTests.h"
#include "test_Vector.h"
#include "utils/AMP_MPI.h"
#include "vectors/petsc/ManagedPetscVector.h"
#include "vectors/petsc/NativePetscVector.h"


/// \cond UNDOCUMENTED

namespace AMP {
namespace unit_test {

template <typename FACTORY>
class PetscCloneFactory
{
public:
    typedef AMP::LinearAlgebra::Vector vector;

    static AMP::LinearAlgebra::Vector::shared_ptr getNativeVector()
    {
        return FACTORY::getNativeVector();
    }

    static void destroyNativeVector( AMP::LinearAlgebra::NativePetscVector &rhs )
    {
        FACTORY::destroyNativeVector( rhs );
    }

    static void destroyNativeVector( AMP::LinearAlgebra::Vector::shared_ptr rhs )
    {
        FACTORY::destroyNativeVector( rhs );
    }

    static AMP::LinearAlgebra::Vector::shared_ptr getManagedVector()
    {
        return FACTORY::getManagedVector()->cloneVector();
    }
};

template <typename FACTORY>
class PetscViewFactory
{
public:
    typedef AMP::LinearAlgebra::Vector vector;

    static AMP::LinearAlgebra::Vector::shared_ptr getNativeVector()
    {
        return FACTORY::getNativeVector();
    }

    static void destroyNativeVector( AMP::LinearAlgebra::NativePetscVector &rhs )
    {
        FACTORY::destroyNativeVector( rhs );
    }

    static void destroyNativeVector( AMP::LinearAlgebra::Vector::shared_ptr rhs )
    {
        FACTORY::destroyNativeVector( rhs );
    }

    static AMP::LinearAlgebra::Vector::shared_ptr getManagedVector()
    {
        return AMP::LinearAlgebra::PetscVector::view( FACTORY::getManagedVector() );
    }
};

template <typename MANAGED_FACTORY>
class SimplePetscVectorFactory
{
public:
    typedef AMP::LinearAlgebra::Vector vector;

    static AMP::LinearAlgebra::Vector::shared_ptr getNativeVector()
    {
        AMP::LinearAlgebra::Vector::shared_ptr t = getManagedVector();
        size_t localSize                         = t->getLocalSize();
        Vec ans;
        AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
        VecCreate( globalComm.getCommunicator(), &ans );
        VecSetSizes( ans, localSize, PETSC_DECIDE );
        VecSetFromOptions( ans );
        PetscInt N;
        VecGetSize( ans, &N );
        AMP_ASSERT( N == (int) ( t->getGlobalSize() ) );
        int a, b;
        VecGetOwnershipRange( ans, &a, &b );
        AMP_ASSERT( b - a == (int) localSize );
        AMP::shared_ptr<AMP::LinearAlgebra::NativePetscVectorParameters> npvParams(
            new AMP::LinearAlgebra::NativePetscVectorParameters( ans, true ) );
        npvParams->d_Deleteable = true;
        AMP::shared_ptr<AMP::LinearAlgebra::NativePetscVector> retVal(
            new AMP::LinearAlgebra::NativePetscVector( npvParams ) );
        retVal->setVariable( AMP::LinearAlgebra::Variable::shared_ptr(
            new AMP::LinearAlgebra::Variable( "petsc vector" ) ) );
        return retVal;
    }

    static void destroyNativeVector( AMP::LinearAlgebra::NativePetscVector &rhs )
    {
#if ( PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 0 )
        VecDestroy( rhs.getVec() );
#elif ( PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 2 )
        VecDestroy( &rhs.getVec() );
#else
#error Not programmed for this version yet
#endif
    }

    static void destroyNativeVector( AMP::LinearAlgebra::Vector::shared_ptr rhs )
    {
        destroyNativeVector(
            *AMP::dynamic_pointer_cast<AMP::LinearAlgebra::NativePetscVector>( rhs ) );
    }

    static AMP::LinearAlgebra::Vector::shared_ptr getManagedVector()
    {
        return MANAGED_FACTORY::getVector();
    }
};


#ifdef USE_EXT_TRILINOS
template <typename T>
class SimplePetscNativeFactory
    : public SimplePetscVectorFactory<
          SimpleManagedVectorFactory<AMP::LinearAlgebra::ManagedPetscVector>>
{
public:
    typedef T vector;

    static AMP::LinearAlgebra::Variable::shared_ptr getVariable()
    {
        return AMP::LinearAlgebra::Variable::shared_ptr(
            new AMP::LinearAlgebra::Variable( "dummy" ) ); // No associated variable
    }

    static AMP::LinearAlgebra::Vector::shared_ptr getVector() { return getNativeVector(); }
};
#endif
}
}

/// \endcond
#endif
