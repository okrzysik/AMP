#ifndef included_test_PetscVectorFactory
#define included_test_PetscVectorFactory

#include "vectors/testHelpers/VectorFactory.h"
#include "vectors/testHelpers/petsc/PetscVectorTests.h"
#include "vectors/testHelpers/generateVectorFactories.h"

#include "utils/AMP_MPI.h"
#include "vectors/petsc/PetscHelpers.h"
#include "vectors/petsc/ManagedPetscVector.h"
#include "vectors/petsc/NativePetscVector.h"


/// \cond UNDOCUMENTED

namespace AMP {
namespace LinearAlgebra {


class PetscCloneFactory: public PetscVectorFactory
{
public:
    PetscCloneFactory( AMP::shared_ptr<const PetscVectorFactory> factory ): d_factory(factory) {}

    virtual AMP::LinearAlgebra::Vector::shared_ptr getNativeVector() const override
    {
        return d_factory->getNativeVector();
    }

    virtual void destroyNativeVector( AMP::LinearAlgebra::NativePetscVector &rhs ) const override
    {
        d_factory->destroyNativeVector( rhs );
    }

    virtual void destroyNativeVector( AMP::LinearAlgebra::Vector::shared_ptr rhs ) const override
    {
        d_factory->destroyNativeVector( rhs );
    }

    virtual AMP::LinearAlgebra::Vector::shared_ptr getManagedVector() const override
    {
        return d_factory->getManagedVector()->cloneVector();
    }
    virtual std::string name() const override { return "PetscCloneFactory"; }

private:
    PetscCloneFactory();
    AMP::shared_ptr<const PetscVectorFactory> d_factory;
};


class PetscViewFactory: public PetscVectorFactory
{
public:
    PetscViewFactory( AMP::shared_ptr<const PetscVectorFactory> factory ): d_factory(factory) {}

    virtual AMP::LinearAlgebra::Vector::shared_ptr getNativeVector() const override
    {
        return d_factory->getNativeVector();
    }

    virtual void destroyNativeVector( AMP::LinearAlgebra::NativePetscVector &rhs ) const override
    {
        d_factory->destroyNativeVector( rhs );
    }

    virtual void destroyNativeVector( AMP::LinearAlgebra::Vector::shared_ptr rhs ) const override
    {
        d_factory->destroyNativeVector( rhs );
    }

    virtual AMP::LinearAlgebra::Vector::shared_ptr getManagedVector() const override
    {
        return AMP::LinearAlgebra::PetscVector::view( d_factory->getManagedVector() );
    }

    virtual std::string name() const override { return "PetscViewFactory"; }

private:
    AMP::shared_ptr<const PetscVectorFactory> d_factory;
};


class SimplePetscVectorFactory: public PetscVectorFactory
{
public:

    SimplePetscVectorFactory( AMP::shared_ptr<const VectorFactory> factory ): d_factory(factory) {}

    virtual AMP::LinearAlgebra::Vector::shared_ptr getNativeVector() const override
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

    virtual void destroyNativeVector( AMP::LinearAlgebra::NativePetscVector &rhs ) const override
    {
        PETSC::vecDestroy( &rhs.getVec() );
    }

    virtual void destroyNativeVector( AMP::LinearAlgebra::Vector::shared_ptr rhs ) const override
    {
        destroyNativeVector(
            *AMP::dynamic_pointer_cast<AMP::LinearAlgebra::NativePetscVector>( rhs ) );
    }

    virtual AMP::LinearAlgebra::Vector::shared_ptr getManagedVector() const override
    {
        return d_factory->getVector();
    }
    virtual std::string name() const override { return "SimplePetscVectorFactory"; }

protected:
    AMP::shared_ptr<const VectorFactory> d_factory;
};


#ifdef USE_EXT_TRILINOS
class SimplePetscNativeFactory : public VectorFactory, SimplePetscVectorFactory
{
public:
    SimplePetscNativeFactory( ):
        SimplePetscVectorFactory( generateVectorFactory( "SimpleManagedVectorFactory<ManagedPetscVector>" ) )
    {
    }

    virtual AMP::LinearAlgebra::Variable::shared_ptr getVariable() const override
    {
        return AMP::LinearAlgebra::Variable::shared_ptr(
            new AMP::LinearAlgebra::Variable( "dummy" ) ); // No associated variable
    }

    virtual AMP::LinearAlgebra::Vector::shared_ptr getVector() const override { return getNativeVector(); }

    virtual std::string name() const override { return "SimplePetscNativeFactory"; }

    virtual AMP::Discretization::DOFManager::shared_ptr getDOFMap() const override { return d_factory->getDOFMap(); }

};


#endif
}
}

/// \endcond
#endif
