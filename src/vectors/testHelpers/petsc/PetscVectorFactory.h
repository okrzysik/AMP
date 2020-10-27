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
class PetscVectorFactory : public VectorFactory
{
public:
    virtual ~PetscVectorFactory() {}
    virtual Vec getVec( AMP::LinearAlgebra::Vector::shared_ptr ) const = 0;

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
    Vec getVec( AMP::LinearAlgebra::Vector::shared_ptr vec ) const override
    {
        return d_factory->getVec( vec );
    }
    AMP::LinearAlgebra::Vector::shared_ptr getVector() const override
    {
        return d_factory->getVector()->cloneVector();
    }
    std::string name() const override { return "PetscCloneFactory<" + d_factory->name() + ">"; }

private:
    std::shared_ptr<const PetscVectorFactory> d_factory;
};


class PetscViewFactory : public PetscVectorFactory
{
public:
    explicit PetscViewFactory( std::shared_ptr<const VectorFactory> factory ) : d_factory( factory )
    {
    }
    Vec getVec( AMP::LinearAlgebra::Vector::shared_ptr vec ) const override
    {
        return std::dynamic_pointer_cast<AMP::LinearAlgebra::ManagedPetscVector>( vec )->getVec();
    }
    AMP::LinearAlgebra::Vector::shared_ptr getVector() const override
    {
        return std::dynamic_pointer_cast<AMP::LinearAlgebra::ManagedPetscVector>(
            AMP::LinearAlgebra::PetscVector::view( d_factory->getVector() ) );
    }
    std::string name() const override { return "PetscViewFactory<" + d_factory->name() + ">"; }

private:
    std::shared_ptr<const VectorFactory> d_factory;
};


class ManagedPetscVectorFactory : public PetscVectorFactory
{
public:
    ManagedPetscVectorFactory( std::shared_ptr<const VectorFactory> factory ) : d_factory( factory )
    {
    }
    Vec getVec( AMP::LinearAlgebra::Vector::shared_ptr vec ) const override
    {
        return std::dynamic_pointer_cast<AMP::LinearAlgebra::ManagedPetscVector>( vec )->getVec();
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


class NativePetscVectorFactory : public PetscVectorFactory
{
public:
    Vec getVec( AMP::LinearAlgebra::Vector::shared_ptr vec ) const override
    {
        return std::dynamic_pointer_cast<AMP::LinearAlgebra::NativePetscVectorData>(
                   vec->getVectorData() )
            ->getVec();
    }
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


} // namespace LinearAlgebra
} // namespace AMP

/// \endcond
#endif
