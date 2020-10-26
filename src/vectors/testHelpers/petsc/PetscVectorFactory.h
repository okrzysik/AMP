#ifndef included_AMP_test_PetscVectorFactory
#define included_AMP_test_PetscVectorFactory

#include "AMP/utils/AMP_MPI.h"
#include "AMP/vectors/VectorBuilder.h"
#include "AMP/vectors/petsc/ManagedPetscVector.h"
#include "AMP/vectors/petsc/NativePetscVectorData.h"
#include "AMP/vectors/petsc/PetscHelpers.h"
#include "AMP/vectors/testHelpers/VectorFactory.h"
#include "AMP/vectors/testHelpers/generateVectorFactories.h"

#include "petscvec.h"


/// \cond UNDOCUMENTED

namespace AMP {
namespace LinearAlgebra {


/**
 * \class PetscVectorFactory
 * \brief A helper class to generate vectors
 */
class PetscVectorFactory
{
public:
    virtual ~PetscVectorFactory() {}
    virtual AMP::LinearAlgebra::Vector::shared_ptr getNativeVector() const = 0;
    virtual AMP::LinearAlgebra::Vector::shared_ptr getManagedVector() const = 0;
    virtual std::string name() const = 0;
protected:
    PetscVectorFactory() {}
    PetscVectorFactory( const PetscVectorFactory & );
};


class PetscCloneFactory : public PetscVectorFactory
{
public:
    PetscCloneFactory() = delete;
    explicit PetscCloneFactory( std::shared_ptr<const PetscVectorFactory> factory )
        : d_factory( factory )
    {
    }
    AMP::LinearAlgebra::Vector::shared_ptr getNativeVector() const override
    {
        return d_factory->getNativeVector();
    }
    AMP::LinearAlgebra::Vector::shared_ptr getManagedVector() const override
    {
        return d_factory->getManagedVector()->cloneVector();
    }
    std::string name() const override { return "PetscCloneFactory"; }
private:
    std::shared_ptr<const PetscVectorFactory> d_factory;
};


class PetscViewFactory : public PetscVectorFactory
{
public:
    explicit PetscViewFactory( std::shared_ptr<const PetscVectorFactory> factory )
        : d_factory( factory )
    {
    }
    AMP::LinearAlgebra::Vector::shared_ptr getNativeVector() const override
    {
        return d_factory->getNativeVector();
    }
    AMP::LinearAlgebra::Vector::shared_ptr getManagedVector() const override
    {
        return AMP::LinearAlgebra::PetscVector::view( d_factory->getManagedVector() )
            ->getManagedVec();
    }
    std::string name() const override { return "PetscViewFactory"; }
private:
    std::shared_ptr<const PetscVectorFactory> d_factory;
};


class SimplePetscVectorFactory : public PetscVectorFactory
{
public:
    explicit SimplePetscVectorFactory( std::shared_ptr<const VectorFactory> factory )
        : d_factory( factory )
    {
    }
    AMP::LinearAlgebra::Vector::shared_ptr getNativeVector() const override
    {
        auto t           = getManagedVector();
        size_t localSize = t->getLocalSize();
        Vec ans;
        AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
        VecCreate( t->getComm().getCommunicator(), &ans );
        VecSetSizes( ans, localSize, PETSC_DECIDE );
        VecSetFromOptions( ans );
        PetscInt N;
        VecGetSize( ans, &N );
        AMP_ASSERT( N == (int) ( t->getGlobalSize() ) );
        int a, b;
        VecGetOwnershipRange( ans, &a, &b );
        AMP_ASSERT( b - a == (int) localSize );
        auto retVal = createVector( ans, true );
        retVal->setVariable( std::make_shared<AMP::LinearAlgebra::Variable>( "petsc vector" ) );
        return retVal;
    }
    AMP::LinearAlgebra::Vector::shared_ptr getManagedVector() const override
    {
        return d_factory->getVector();
    }
    std::string name() const override { return "SimplePetscVectorFactory"; }
protected:
    std::shared_ptr<const VectorFactory> d_factory;
};


class SimplePetscNativeFactory : public VectorFactory, SimplePetscVectorFactory
{
public:
    SimplePetscNativeFactory()
        : SimplePetscVectorFactory(
              generateVectorFactory( "ManagedPetscVectorFactory<SimpleVectorFactory<45,true>>" ) )
    {
    }
    AMP::LinearAlgebra::Vector::shared_ptr getVector() const override
    {
        return getNativeVector();
    }
    std::string name() const override { return "SimplePetscNativeFactory"; }
};


class ManagedPetscVectorFactory : public VectorFactory
{
public:
    ManagedPetscVectorFactory( std::shared_ptr<const VectorFactory> factory ) : d_factory( factory )
    {
    }
    AMP::LinearAlgebra::Vector::shared_ptr getVector() const override
    {
        auto engine = d_factory->getVector();
        auto retval = std::make_shared<ManagedPetscVector>( engine );
        retval->setVariable( std::make_shared<AMP::LinearAlgebra::Variable>( "Test Vector" ) );
        return retval;
    }
    std::string name() const override
    {
        return "ManagedPetscVectorFactory<" + d_factory->name() + ">";
    }
private:
    std::shared_ptr<const VectorFactory> d_factory;
};


class NativePetscVectorFactory : public VectorFactory
{
public:
    AMP::LinearAlgebra::Vector::shared_ptr getVector() const override
    {
        Vec v;
        AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
        VecCreate( globalComm.getCommunicator(), &v );
        PetscInt local_size = 15;
        VecSetSizes( v, local_size, PETSC_DECIDE );
        VecSetType( v, VECMPI ); // this line will have to be modified for the no mpi and cuda cases
        auto newVec = createVector( v, true );
        VecSetFromOptions( v );
        newVec->getVectorData()->assemble();
        newVec->setVariable(
            std::make_shared<AMP::LinearAlgebra::Variable>( "Test NativePetscVector" ) );
        return newVec;
    }
    std::string name() const override { return "NativePetscVectorFactory"; }
};

template<typename T>
class PetscManagedVectorFactory : public VectorFactory
{
public:
    AMP::LinearAlgebra::Vector::shared_ptr getVector() const override
    {
        Vec v;
        AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
        VecCreate( globalComm.getCommunicator(), &v );
        VecSetSizes( v, 15, PETSC_DECIDE );
        auto newVec = createVector( v, true );
        VecSetFromOptions( v );
        newVec->getVectorData()->assemble();
        auto retval = std::make_shared<T>( newVec );
        retval->setVariable( std::make_shared<AMP::LinearAlgebra::Variable>( "Test Vector" ) );
        return retval;
    }
    std::string name() const override { return "PetscManagedVectorFactory"; }
};


} // namespace LinearAlgebra
} // namespace AMP

/// \endcond
#endif
