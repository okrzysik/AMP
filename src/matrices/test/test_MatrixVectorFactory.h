#ifndef included_test_MatrixVectorFactory
#define included_test_MatrixVectorFactory

#include "discretization/DOF_Manager.h"
#include "matrices/Matrix.h"
#include "utils/UnitTest.h"
#include "vectors/Vector.h"

#ifdef USE_EXT_PETSC
#include "matrices/petsc/PetscMatrix.h"
#include "vectors/petsc/ManagedPetscVector.h"
#include "vectors/petsc/NativePetscVector.h"
#endif


namespace AMP {
namespace unit_test {


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
};


#if defined( USE_EXT_PETSC ) && defined( USE_EXT_PETSC )

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
        MatGetVecs( m, &v, nullptr );
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
        MatGetVecs( m, &v, nullptr );
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
};

#endif
}
}
#endif
