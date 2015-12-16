#include "utils/AMP_MPI.h"
#include "vectors/sundials/ManagedSundialsVector.h"

#include "VectorUnitTest.h"
#include "test_SundialsVectorTests.h"

/// \cond UNDOCUMENTED

namespace AMP {
namespace unit_test {

template <typename FACTORY>
class PetscCloneFactory {
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
class PetscViewFactory {
public:
    typedef AMP::LinearAlgebra::Vector vector;

    static AMP::LinearAlgebra::Vector::shared_ptr getNativeVecto()
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
class SimplePetscVectorFactory {
public:
    typedef AMP::LinearAlgebra::Vector vector;

    static AMP::LinearAlgebra::Vector::shared_ptr getNativeVector()
    {
        AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
        AMP::LinearAlgebra::Vector::shared_ptr t = getManagedVector();
        Vec ans;
        VecCreate( globalComm, &ans );
        VecSetSizes( ans, t->getLocalSize(), PETSC_DECIDE );
        VecSetFromOptions( ans );
        return AMP::LinearAlgebra::Vector::shared_ptr(
            new AMP::LinearAlgebra::NativePetscVector( ans, globalComm, true ) );
    }

    static void destroyNativeVector( AMP::LinearAlgebra::NativePetscVector &rhs )
    {
        VecDestroy( rhs.getVec() );
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

template <typename T>
class SimplePetscNativeFactory
    : public SimplePetscVectorFactory<
          SimpleManagedVectorFactory<AMP::LinearAlgebra::ManagedPetscVector>> {
public:
    typedef T vector;

    static AMP::LinearAlgebra::Variable::shared_ptr getVariable( UnitTest & )
    {
        return AMP::LinearAlgebra::Variable::shared_ptr(); // No associated variable
    }

    static AMP::LinearAlgebra::Vector::shared_ptr getVector( UnitTest &u )
    {
        return getNativeVector( u );
    }
};
};
}
}

/// \endcond
