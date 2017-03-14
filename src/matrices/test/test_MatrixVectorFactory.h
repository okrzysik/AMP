#ifndef included_test_MatrixVectorFactory
#define included_test_MatrixVectorFactory

#include "discretization/DOF_Manager.h"
#include "matrices/Matrix.h"
#include "utils/UnitTest.h"
#include "vectors/Vector.h"

#if defined(USE_EXT_PETSC) && defined(USE_EXT_TRILINOS)
#include "matrices/petsc/PetscMatrix.h"
#include "vectors/petsc/ManagedPetscVector.h"
#include "vectors/petsc/NativePetscVector.h"
#endif


namespace AMP {
namespace LinearAlgebra {


AMP::LinearAlgebra::Matrix::shared_ptr global_cached_matrix =
    AMP::LinearAlgebra::Matrix::shared_ptr();


// Classes to serve as the vector factories
class AmpInterfaceLeftVectorFactory
{
public:
    typedef AMP::LinearAlgebra::Vector vector;

    static AMP::LinearAlgebra::Variable::shared_ptr getVariable()
    {
        return AMP::LinearAlgebra::Variable::shared_ptr(
            new AMP::LinearAlgebra::Variable( "left" ) );
    }

    static AMP::LinearAlgebra::Vector::shared_ptr getVector()
    {
        PROFILE_START( "AmpInterfaceLeftVectorFactory::getVector" );
        AMP::LinearAlgebra::Matrix::shared_ptr matrix = global_cached_matrix;
        AMP_ASSERT( global_cached_matrix != nullptr );
        AMP::LinearAlgebra::Vector::shared_ptr vector = matrix->getLeftVector();
        vector->setVariable( getVariable() );
        PROFILE_STOP( "AmpInterfaceLeftVectorFactory::getVector" );
        return vector;
    }
    static std::string name() { return "AmpInterfaceLeftVectorFactory"; }
};


class AmpInterfaceRightVectorFactory
{
public:
    typedef AMP::LinearAlgebra::Vector vector;

    static AMP::LinearAlgebra::Variable::shared_ptr getVariable()
    {
        return AMP::LinearAlgebra::Variable::shared_ptr(
            new AMP::LinearAlgebra::Variable( "right" ) );
    }

    static AMP::LinearAlgebra::Vector::shared_ptr getVector()
    {
        PROFILE_START( "AmpInterfaceRightVectorFactory::getVector" );
        AMP::LinearAlgebra::Matrix::shared_ptr matrix = global_cached_matrix;
        AMP_ASSERT( global_cached_matrix != nullptr );
        AMP::LinearAlgebra::Vector::shared_ptr vector = matrix->getRightVector();
        vector->setVariable( getVariable() );
        PROFILE_STOP( "AmpInterfaceRightVectorFactory::getVector" );
        return vector;
    }
    static std::string name() { return "AmpInterfaceRightVectorFactory"; }
};


#if defined(USE_EXT_PETSC) && defined(USE_EXT_TRILINOS)

class PETScInterfaceLeftVectorFactory
{
public:
    typedef AMP::LinearAlgebra::ManagedPetscVector vector;

    static AMP::LinearAlgebra::Variable::shared_ptr getVariable()
    {
        return AMP::LinearAlgebra::Variable::shared_ptr(
            new AMP::LinearAlgebra::Variable( "petsc_left" ) );
    }

    static AMP::LinearAlgebra::Vector::shared_ptr getVector()
    {
        PROFILE_START( "PETScInterfaceLeftVectorFactory::getVector" );
        AMP_ASSERT( global_cached_matrix != nullptr );
        AMP::LinearAlgebra::Matrix::shared_ptr matrix =
            AMP::LinearAlgebra::PetscMatrix::createView( global_cached_matrix );
        ::Mat m = matrix->castTo<AMP::LinearAlgebra::PetscMatrix>().getMat();
        ::Vec v;
        DISABLE_WARNINGS
        MatGetVecs( m, &v, nullptr );
        ENABLE_WARNINGS
        AMP::shared_ptr<AMP::LinearAlgebra::NativePetscVectorParameters> p(
            new AMP::LinearAlgebra::NativePetscVectorParameters( v, true ) );
        AMP::LinearAlgebra::Vector::shared_ptr vector(
            new AMP::LinearAlgebra::NativePetscVector( p ) );
        vector->setVariable( getVariable() );
        PROFILE_STOP( "PETScInterfaceLeftVectorFactory::getVector" );
        return vector;
    }

    static AMP::LinearAlgebra::Vector::shared_ptr getNativeVector() { return getVector(); }

    static AMP::LinearAlgebra::Vector::shared_ptr getManagedVector()
    {
        AMP::LinearAlgebra::Matrix::shared_ptr matrix = global_cached_matrix;
        AMP_ASSERT( global_cached_matrix != nullptr );
        return AMP::LinearAlgebra::PetscVector::view( matrix->getLeftVector() );
    }
    static std::string name() { return "PETScInterfaceLeftVectorFactory"; }
};


class PETScInterfaceRightVectorFactory
{
public:
    typedef AMP::LinearAlgebra::ManagedPetscVector vector;

    static AMP::LinearAlgebra::Variable::shared_ptr getVariable()
    {
        return AMP::LinearAlgebra::Variable::shared_ptr(
            new AMP::LinearAlgebra::Variable( "petsc_right" ) );
    }

    static AMP::LinearAlgebra::Vector::shared_ptr getVector()
    {
        PROFILE_START( "PETScInterfaceRightVectorFactory::getVector" );
        AMP_ASSERT( global_cached_matrix != nullptr );
        AMP::LinearAlgebra::Matrix::shared_ptr matrix =
            AMP::LinearAlgebra::PetscMatrix::createView( global_cached_matrix );
        ::Mat m = matrix->castTo<AMP::LinearAlgebra::PetscMatrix>().getMat();
        ::Vec v;
        DISABLE_WARNINGS
        MatGetVecs( m, &v, nullptr );
        ENABLE_WARNINGS
        AMP::shared_ptr<AMP::LinearAlgebra::NativePetscVectorParameters> p(
            new AMP::LinearAlgebra::NativePetscVectorParameters( v, true ) );
        AMP::LinearAlgebra::Vector::shared_ptr vector(
            new AMP::LinearAlgebra::NativePetscVector( p ) );
        vector->setVariable( getVariable() );
        PROFILE_STOP( "PETScInterfaceRightVectorFactory::getVector" );
        return vector;
    }

    static AMP::LinearAlgebra::Vector::shared_ptr getNativeVector() { return getVector(); }

    static AMP::LinearAlgebra::Vector::shared_ptr getManagedVector()
    {
        AMP_ASSERT( global_cached_matrix != nullptr );
        AMP::LinearAlgebra::Matrix::shared_ptr matrix = global_cached_matrix;
        return AMP::LinearAlgebra::PetscVector::view( matrix->getRightVector() );
    }
    static std::string name() { return "PETScInterfaceRightVectorFactory"; }
};

#endif
}
}
#endif
