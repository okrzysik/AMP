#ifndef included_AMP_test_MatrixVectorFactory
#define included_AMP_test_MatrixVectorFactory

#include "AMP/discretization/DOF_Manager.h"
#include "AMP/matrices/Matrix.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"
#include "AMP/vectors/testHelpers/VectorFactory.h"

#if defined( USE_EXT_PETSC )
#include "AMP/matrices/petsc/PetscMatrix.h"
#include "AMP/vectors/petsc/PetscHelpers.h"
#include "AMP/vectors/petsc/PetscVector.h"
#include "AMP/vectors/testHelpers/petsc/PetscVectorFactory.h"
#include "petscmat.h"
#endif

#include "ProfilerApp.h"

namespace AMP {
namespace LinearAlgebra {


// Classes to serve as the vector factories
class AmpInterfaceLeftVectorFactory : public VectorFactory
{
public:
    explicit AmpInterfaceLeftVectorFactory( AMP::LinearAlgebra::Matrix::shared_ptr matrix )
        : d_matrix( matrix )
    {
    }
    AMP::LinearAlgebra::Vector::shared_ptr getVector() const override
    {
        PROFILE_START( "AmpInterfaceLeftVectorFactory::getVector" );
        auto vector = d_matrix->getLeftVector();
        vector->setVariable( std::make_shared<AMP::LinearAlgebra::Variable>( "left" ) );
        PROFILE_STOP( "AmpInterfaceLeftVectorFactory::getVector" );
        return vector;
    }
    std::string name() const override
    {
        return "AmpInterfaceLeftVectorFactory<" + getVector()->type() + ">";
    }

private:
    AMP::LinearAlgebra::Matrix::shared_ptr d_matrix;
};


class AmpInterfaceRightVectorFactory : public VectorFactory
{
public:
    explicit AmpInterfaceRightVectorFactory( AMP::LinearAlgebra::Matrix::shared_ptr matrix )
        : d_matrix( matrix )
    {
    }
    AMP::LinearAlgebra::Vector::shared_ptr getVector() const override
    {
        PROFILE_START( "AmpInterfaceRightVectorFactory::getVector" );
        auto vector = d_matrix->getRightVector();
        vector->setVariable( std::make_shared<AMP::LinearAlgebra::Variable>( "right" ) );
        PROFILE_STOP( "AmpInterfaceRightVectorFactory::getVector" );
        return vector;
    }
    std::string name() const override
    {
        return "AmpInterfaceRightVectorFactory<" + getVector()->type() + ">";
    }

private:
    AMP::LinearAlgebra::Matrix::shared_ptr d_matrix;
};


#if defined( USE_EXT_PETSC )

class PETScInterfaceLeftVectorFactory : public PetscVectorFactory
{
public:
    explicit PETScInterfaceLeftVectorFactory( AMP::LinearAlgebra::Matrix::shared_ptr matrix )
        : d_matrix( matrix )
    {
    }
    AMP::LinearAlgebra::Vector::shared_ptr getVector() const override
    {
        PROFILE_START( "PETScInterfaceLeftVectorFactory::getVector" );
        auto view = AMP::LinearAlgebra::PetscMatrix::view( d_matrix );
        ::Mat m   = view->getMat();
        ::Vec v;
        DISABLE_WARNINGS
        MatGetVecs( m, &v, nullptr );
        ENABLE_WARNINGS
        auto vector = createVector( v, true );
        vector->setVariable( std::make_shared<AMP::LinearAlgebra::Variable>( "petsc_left" ) );
        PROFILE_STOP( "PETScInterfaceLeftVectorFactory::getVector" );
        return vector;
    }
    std::shared_ptr<Vec> getVec( AMP::LinearAlgebra::Vector::shared_ptr vec ) const override
    {
        auto data = std::dynamic_pointer_cast<NativePetscVectorData>( vec->getVectorData() );
        auto ptr  = std::make_shared<Vec>( data->getVec() );
        return ptr;
    }
    std::string name() const override { return "PETScInterfaceLeftVectorFactory"; };

private:
    AMP::LinearAlgebra::Matrix::shared_ptr d_matrix;
};


class PETScInterfaceRightVectorFactory : public PetscVectorFactory
{
public:
    explicit PETScInterfaceRightVectorFactory( AMP::LinearAlgebra::Matrix::shared_ptr matrix )
        : d_matrix( matrix )
    {
    }
    AMP::LinearAlgebra::Vector::shared_ptr getVector() const override
    {
        PROFILE_START( "PETScInterfaceRightVectorFactory::getVector" );
        auto view = AMP::LinearAlgebra::PetscMatrix::view( d_matrix );
        ::Mat m   = view->getMat();
        ::Vec v;
        DISABLE_WARNINGS
        MatGetVecs( m, &v, nullptr );
        ENABLE_WARNINGS
        auto vector = createVector( v, true );
        vector->setVariable( std::make_shared<AMP::LinearAlgebra::Variable>( "petsc_right" ) );
        PROFILE_STOP( "PETScInterfaceRightVectorFactory::getVector" );
        return vector;
    }
    std::shared_ptr<Vec> getVec( AMP::LinearAlgebra::Vector::shared_ptr vec ) const override
    {
        auto data = std::dynamic_pointer_cast<NativePetscVectorData>( vec->getVectorData() );
        auto ptr  = std::make_shared<Vec>( data->getVec() );
        return ptr;
    }
    std::string name() const override { return "PETScInterfaceRightVectorFactory"; }

private:
    AMP::LinearAlgebra::Matrix::shared_ptr d_matrix;
};

#endif

} // namespace LinearAlgebra
} // namespace AMP

#endif
