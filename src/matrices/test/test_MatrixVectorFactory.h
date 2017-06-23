#ifndef included_test_MatrixVectorFactory
#define included_test_MatrixVectorFactory

#include "discretization/DOF_Manager.h"
#include "matrices/Matrix.h"
#include "utils/UnitTest.h"
#include "vectors/Vector.h"

#include <vectors/testHelpers/VectorFactory.h>

#if defined(USE_EXT_PETSC) && defined(USE_EXT_TRILINOS)
#include "matrices/petsc/PetscMatrix.h"
#include "vectors/petsc/ManagedPetscVector.h"
#include "vectors/petsc/NativePetscVector.h"
#include <vectors/testHelpers/petsc/PetscVectorFactory.h>
#endif




namespace AMP {
namespace LinearAlgebra {


AMP::LinearAlgebra::Matrix::shared_ptr global_cached_matrix =
    AMP::LinearAlgebra::Matrix::shared_ptr();


// Classes to serve as the vector factories
class AmpInterfaceLeftVectorFactory: public VectorFactory
{
public:

    virtual AMP::LinearAlgebra::Variable::shared_ptr getVariable() const override
    {
        return AMP::LinearAlgebra::Variable::shared_ptr(
            new AMP::LinearAlgebra::Variable( "left" ) );
    }

    virtual AMP::LinearAlgebra::Vector::shared_ptr getVector() const override
    {
        PROFILE_START( "AmpInterfaceLeftVectorFactory::getVector" );
        AMP::LinearAlgebra::Matrix::shared_ptr matrix = global_cached_matrix;
        AMP_ASSERT( global_cached_matrix != nullptr );
        AMP::LinearAlgebra::Vector::shared_ptr vector = matrix->getLeftVector();
        vector->setVariable( getVariable() );
        PROFILE_STOP( "AmpInterfaceLeftVectorFactory::getVector" );
        return vector;
    }
    virtual std::string name() const override { return "AmpInterfaceLeftVectorFactory"; }

    virtual AMP::Discretization::DOFManager::shared_ptr getDOFMap() const override
    {
        return getVector()->getDOFManager();
    }
};


class AmpInterfaceRightVectorFactory: public VectorFactory
{
public:

    virtual AMP::LinearAlgebra::Variable::shared_ptr getVariable() const override
    {
        return AMP::LinearAlgebra::Variable::shared_ptr(
            new AMP::LinearAlgebra::Variable( "right" ) );
    }

    virtual AMP::LinearAlgebra::Vector::shared_ptr getVector() const override
    {
        PROFILE_START( "AmpInterfaceRightVectorFactory::getVector" );
        AMP::LinearAlgebra::Matrix::shared_ptr matrix = global_cached_matrix;
        AMP_ASSERT( global_cached_matrix != nullptr );
        AMP::LinearAlgebra::Vector::shared_ptr vector = matrix->getRightVector();
        vector->setVariable( getVariable() );
        PROFILE_STOP( "AmpInterfaceRightVectorFactory::getVector" );
        return vector;
    }
    virtual std::string name() const override { return "AmpInterfaceRightVectorFactory"; }

    virtual AMP::Discretization::DOFManager::shared_ptr getDOFMap() const override
    {
        return getVector()->getDOFManager();
    }
};


#if defined(USE_EXT_PETSC) && defined(USE_EXT_TRILINOS)

class PETScInterfaceLeftVectorFactory: public VectorFactory, PetscVectorFactory
{
public:

    virtual AMP::LinearAlgebra::Variable::shared_ptr getVariable() const override
    {
        return AMP::LinearAlgebra::Variable::shared_ptr(
            new AMP::LinearAlgebra::Variable( "petsc_left" ) );
    }

    virtual AMP::LinearAlgebra::Vector::shared_ptr getVector() const override
    {
        PROFILE_START( "PETScInterfaceLeftVectorFactory::getVector" );
        AMP_ASSERT( global_cached_matrix != nullptr );
        auto matrix = AMP::dynamic_pointer_cast<AMP::LinearAlgebra::PetscMatrix>(
            AMP::LinearAlgebra::PetscMatrix::createView( global_cached_matrix ) );
        ::Mat m = matrix->getMat();
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

    virtual AMP::LinearAlgebra::Vector::shared_ptr getNativeVector() const override { return getVector(); }

    virtual AMP::LinearAlgebra::Vector::shared_ptr getManagedVector() const override
    {
        AMP::LinearAlgebra::Matrix::shared_ptr matrix = global_cached_matrix;
        AMP_ASSERT( global_cached_matrix != nullptr );
        return AMP::LinearAlgebra::PetscVector::view( matrix->getLeftVector() );
    }

    virtual std::string name() const override { return "PETScInterfaceLeftVectorFactory"; }

    virtual AMP::Discretization::DOFManager::shared_ptr getDOFMap() const override
    {
        return getVector()->getDOFManager();
    }
};


class PETScInterfaceRightVectorFactory: public VectorFactory,  PetscVectorFactory
{
public:

    virtual AMP::LinearAlgebra::Variable::shared_ptr getVariable() const override
    {
        return AMP::LinearAlgebra::Variable::shared_ptr(
            new AMP::LinearAlgebra::Variable( "petsc_right" ) );
    }

    virtual AMP::LinearAlgebra::Vector::shared_ptr getVector() const override
    {
        PROFILE_START( "PETScInterfaceRightVectorFactory::getVector" );
        AMP_ASSERT( global_cached_matrix != nullptr );
        auto matrix = AMP::dynamic_pointer_cast<AMP::LinearAlgebra::PetscMatrix>(
            AMP::LinearAlgebra::PetscMatrix::createView( global_cached_matrix ) );
        ::Mat m = matrix->getMat();
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

    virtual AMP::LinearAlgebra::Vector::shared_ptr getNativeVector() const override { return getVector(); }

    virtual AMP::LinearAlgebra::Vector::shared_ptr getManagedVector() const override
    {
        AMP_ASSERT( global_cached_matrix != nullptr );
        AMP::LinearAlgebra::Matrix::shared_ptr matrix = global_cached_matrix;
        return AMP::LinearAlgebra::PetscVector::view( matrix->getRightVector() );
    }

    virtual std::string name() const override { return "PETScInterfaceRightVectorFactory"; }

    virtual AMP::Discretization::DOFManager::shared_ptr getDOFMap() const override
    {
        return getVector()->getDOFManager();
    }
};

#endif

} // LinearAlgebra namespace
} // AMP namespace

#endif
