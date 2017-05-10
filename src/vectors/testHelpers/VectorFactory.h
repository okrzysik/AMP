#ifndef included_test_VectorFactory
#define included_test_VectorFactory

#include "utils/AMP_MPI.h"
#include "vectors/ArrayVector.h"
#include "vectors/ManagedVector.h"
#include "vectors/MultiVariable.h"
#include "vectors/MultiVector.h"
#include "vectors/NullVector.h"
#include "vectors/SimpleVector.h"
#include "vectors/Variable.h"
#include "vectors/testHelpers/VectorTests.h"
#ifdef USE_EXT_PETSC
#include <vectors/petsc/NativePetscVector.h>
#endif
#ifdef USE_EXT_TRILINOS
#include <vectors/trilinos/epetra/EpetraVectorEngine.h>
#endif


namespace AMP {
namespace LinearAlgebra {



class CloneFactory: public VectorFactory
{
public:
    explicit CloneFactory( AMP::shared_ptr<const VectorFactory> factory ): d_factory(factory) {}

    virtual AMP::LinearAlgebra::Variable::shared_ptr getVariable() const override
    {
        return d_factory->getVariable()->cloneVariable( "noname" );
    }

    virtual AMP::LinearAlgebra::Vector::shared_ptr getVector() const override
    {
        AMP::LinearAlgebra::Vector::shared_ptr vec = d_factory->getVector();
        return vec->cloneVector();
    }

    virtual std::string name() const override { return "CloneFactory<" + d_factory->name() + ">"; }

    virtual AMP::Discretization::DOFManager::shared_ptr getDOFMap() const override
    {
        return d_factory->getDOFMap();
    }

private:
    CloneFactory();
    AMP::shared_ptr<const VectorFactory> d_factory;
};


class NullVectorFactory: public VectorFactory
{
public:
    typedef AMP::LinearAlgebra::NullVector vector;

    virtual AMP::LinearAlgebra::Variable::shared_ptr getVariable() const override
    {
        return AMP::LinearAlgebra::Variable::shared_ptr(
            new AMP::LinearAlgebra::Variable( "null" ) );
    }

    virtual AMP::LinearAlgebra::Vector::shared_ptr getVector() const override
    {
        return AMP::LinearAlgebra::NullVector::create( "null" );
    }

    virtual std::string name() const override { return "NullVectorFactory"; }

    virtual AMP::Discretization::DOFManager::shared_ptr getDOFMap() const override
    {
        return AMP::Discretization::DOFManager::shared_ptr();
    }
};


template <class TYPE=double>
class SimpleVectorFactory: public VectorFactory
{
public:

    SimpleVectorFactory( int i, bool global ): I(i), GLOBAL(global) {}

    virtual AMP::LinearAlgebra::Variable::shared_ptr getVariable() const override
    {
        return AMP::LinearAlgebra::Variable::shared_ptr(
            new AMP::LinearAlgebra::Variable( "simple" ) );
    }

    virtual AMP::LinearAlgebra::Vector::shared_ptr getVector() const override
    {
        AMP::LinearAlgebra::Vector::shared_ptr vec;
        if ( GLOBAL )
            vec = AMP::LinearAlgebra::SimpleVector<TYPE>::create(
                I, getVariable(), AMP_MPI( AMP_COMM_WORLD ) );
        else
            vec = AMP::LinearAlgebra::SimpleVector<TYPE>::create( I, getVariable() );
        return vec;
    }

    virtual std::string name() const override { return "SimpleVectorFactory"; }

    virtual AMP::Discretization::DOFManager::shared_ptr getDOFMap() const override
    {
        return getVector()->getDOFManager();
    }

private:
    SimpleVectorFactory();
    int I;
    bool GLOBAL;
};


template <class TYPE=double>
class ArrayVectorFactory: public VectorFactory
{
public:

    ArrayVectorFactory( int d, int i, bool global ): D(d), I(i), GLOBAL(global) {}

    virtual AMP::LinearAlgebra::Variable::shared_ptr getVariable() const override
    {
        return AMP::LinearAlgebra::Variable::shared_ptr(
            new AMP::LinearAlgebra::Variable( "array" ) );
    }

    virtual AMP::LinearAlgebra::Vector::shared_ptr getVector() const override
    {
        AMP::LinearAlgebra::Vector::shared_ptr vec;
        if ( GLOBAL )
            vec = AMP::LinearAlgebra::ArrayVector<TYPE>::create(
                std::vector<size_t>( D, I ), getVariable(), AMP_MPI( AMP_COMM_WORLD ) );
        else
            vec = AMP::LinearAlgebra::ArrayVector<TYPE>::create( std::vector<size_t>( D, I ),
                                                              getVariable() );
        return vec;
    }

    virtual std::string name() const override { return "ArrayVectorFactory"; }

    virtual AMP::Discretization::DOFManager::shared_ptr getDOFMap() const override
    {
        return getVector()->getDOFManager();
    }
private:
    ArrayVectorFactory();
    int D, I;
    bool GLOBAL;
};


#ifdef USE_EXT_TRILINOS
template <typename TYPE>
class SimpleManagedVectorFactory: public VectorFactory
{
public:
    SimpleManagedVectorFactory() {}

    virtual AMP::LinearAlgebra::Variable::shared_ptr getVariable() const override
    {
        return AMP::LinearAlgebra::Variable::shared_ptr(
            new AMP::LinearAlgebra::Variable( "..." ) );
    }

    virtual AMP::LinearAlgebra::Vector::shared_ptr getVector() const override
    {
        const int num_local = 210;
        AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
        AMP::shared_ptr<AMP::LinearAlgebra::EpetraVectorEngineParameters> epetraParams(
            new AMP::LinearAlgebra::EpetraVectorEngineParameters(
                num_local, num_local * globalComm.getSize(), globalComm ) );
        AMP::shared_ptr<AMP::LinearAlgebra::ManagedVectorParameters> managedParams(
            new AMP::LinearAlgebra::ManagedVectorParameters );
        AMP::LinearAlgebra::VectorEngine::BufferPtr buffer( new std::vector<double>( 120 ) );
        managedParams->d_Engine = AMP::LinearAlgebra::VectorEngine::shared_ptr(
            new AMP::LinearAlgebra::EpetraVectorEngine( epetraParams, buffer ) );
        managedParams->d_CommList =
            AMP::LinearAlgebra::CommunicationList::createEmpty( 210, globalComm );
        managedParams->d_DOFManager = AMP::Discretization::DOFManager::shared_ptr(
            new AMP::Discretization::DOFManager( 210, globalComm ) );
        AMP::shared_ptr<TYPE> retval( new TYPE( managedParams ) );
        retval->setVariable( AMP::LinearAlgebra::Variable::shared_ptr(
            new AMP::LinearAlgebra::Variable( "Test Vector" ) ) );
        return retval;
    }

    virtual std::string name() const override { return "SimpleManagedVectorFactory"; }

    virtual AMP::Discretization::DOFManager::shared_ptr getDOFMap() const override
    {
        return getVector()->getDOFManager();
    }
};
#endif


#ifdef USE_EXT_PETSC
class NativePetscVectorFactory: public VectorFactory
{
public:

    virtual AMP::LinearAlgebra::Variable::shared_ptr getVariable() const override
    {
        return AMP::LinearAlgebra::Variable::shared_ptr(); // no variable.....
    }

    virtual AMP::LinearAlgebra::Vector::shared_ptr getVector() const override
    {
        Vec v;
        AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
        VecCreate( globalComm.getCommunicator(), &v );
        PetscInt local_size = 15;
        VecSetSizes( v, local_size, PETSC_DECIDE );
        VecSetType( v, VECMPI );  // this line will have to be modified for the no mpi and cuda cases
        AMP::shared_ptr<AMP::LinearAlgebra::NativePetscVectorParameters> npvParams(
            new AMP::LinearAlgebra::NativePetscVectorParameters( v, true ) );
        AMP::shared_ptr<AMP::LinearAlgebra::NativePetscVector> newVec( new AMP::LinearAlgebra::NativePetscVector( npvParams ) );
        VecSetFromOptions( v );
        newVec->assemble();
        newVec->setVariable( AMP::LinearAlgebra::Variable::shared_ptr(
            new AMP::LinearAlgebra::Variable( "Test NativePetscVector" ) ) );
        return newVec;
    }

    virtual std::string name() const override { return "NativePetscVectorFactory"; }

    virtual AMP::Discretization::DOFManager::shared_ptr getDOFMap() const override
    {
        return getVector()->getDOFManager();
    }
};

template <typename T>
class PetscManagedVectorFactory: public VectorFactory
{
public:
    typedef T vector;

    virtual AMP::LinearAlgebra::Variable::shared_ptr getVariable() const override
    {
        return AMP::LinearAlgebra::Variable::shared_ptr(); // no variable.....
    }

    virtual AMP::LinearAlgebra::Vector::shared_ptr getVector() const override
    {
        Vec v;
        AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
        VecCreate( globalComm.getCommunicator(), &v );
        VecSetSizes( v, 15, PETSC_DECIDE );
        AMP::shared_ptr<AMP::LinearAlgebra::NativePetscVectorParameters> npvParams(
            new AMP::LinearAlgebra::NativePetscVectorParameters( v, true ) );
        auto newVec = AMP::make_shared<AMP::LinearAlgebra::NativePetscVector> ( npvParams );
        VecSetFromOptions( v );
        newVec->assemble();
        auto p1        = new AMP::LinearAlgebra::ManagedVectorParameters;
        p1->d_Engine   = newVec;
        p1->d_CommList = AMP::LinearAlgebra::CommunicationList::createEmpty( 210, globalComm );
        AMP::LinearAlgebra::Vector::shared_ptr retval(
            new T( AMP::LinearAlgebra::VectorParameters::shared_ptr( p1 ) ) );
        retval->setVariable( AMP::LinearAlgebra::Variable::shared_ptr(
            new AMP::LinearAlgebra::Variable( "Test Vector" ) ) );
        return retval;
    }

    virtual std::string name() const override { return "PetscManagedVectorFactory"; }
};
#endif

#ifdef USE_EXT_SUNDIALS
class NativeSundialsVectorFactory: public VectorFactory
{
public:

    virtual AMP::LinearAlgebra::Variable::shared_ptr getVariable() const override
    {
        return AMP::LinearAlgebra::Variable::shared_ptr(); // no variable.....
    }

    virtual AMP::LinearAlgebra::Vector::shared_ptr getVector() const override
    {
        AMP::LinearAlgebra::Vector::shared_ptr newVec;
        AMP_ERROR("Not implemented");
        return newVec;
    }

    virtual std::string name() const override { return "NativeSundialsVectorFactory"; }

    virtual AMP::Discretization::DOFManager::shared_ptr getDOFMap() const override
    {
        return getVector()->getDOFManager();
    }
};
#endif

template <typename TYPE>
class ViewFactory: public VectorFactory
{
public:
    explicit ViewFactory( AMP::shared_ptr<const VectorFactory> factory ): d_factory(factory) {}

    virtual AMP::LinearAlgebra::Variable::shared_ptr getVariable() const override
    {
        return d_factory->getVariable();
    }

    virtual AMP::LinearAlgebra::Vector::shared_ptr getVector() const override
    {
        return TYPE::view( d_factory->getVector() );
    }

    virtual std::string name() const override { return "ViewFactory"; }

    virtual AMP::Discretization::DOFManager::shared_ptr getDOFMap() const override
    {
        return d_factory->getDOFMap();
    }

private:
    ViewFactory();
    AMP::shared_ptr<const VectorFactory> d_factory;
};


class MultiVectorFactory: public VectorFactory
{
public:
    MultiVectorFactory( AMP::shared_ptr<const VectorFactory> factory1, int N1,
                        AMP::shared_ptr<const VectorFactory> factory2, int N2 ):
        NUM1(N1), 
        NUM2(N2), 
        FACTORY1(factory1),
        FACTORY2(factory2)
    {
    }

    virtual AMP::LinearAlgebra::Variable::shared_ptr getVariable() const override
    {
        AMP::LinearAlgebra::MultiVariable *newVar = new AMP::LinearAlgebra::MultiVariable( "var1" );
        for ( int i = 0; i != NUM1; i++ )
            newVar->add( FACTORY1->getVariable() );
        for ( int i = 0; i != NUM2; i++ )
            newVar->add( FACTORY2->getVariable() );
        return AMP::LinearAlgebra::Variable::shared_ptr( newVar );
    }

    virtual AMP::LinearAlgebra::Vector::shared_ptr getVector() const override
    {
        AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
        auto retVal = AMP::LinearAlgebra::MultiVector::create( getVariable(), globalComm );
        for ( int i = 0; i != NUM1; i++ )
            retVal->addVector( FACTORY1->getVector() );
        for ( int i = 0; i != NUM2; i++ )
            retVal->addVector( FACTORY2->getVector() );
        return retVal;
    }

    virtual std::string name() const override { return "MultiVectorFactory"; }

    virtual AMP::Discretization::DOFManager::shared_ptr getDOFMap() const override
    {
        return getVector()->getDOFManager();
    }

private:
    MultiVectorFactory();
    int NUM1, NUM2;
    AMP::shared_ptr<const VectorFactory> FACTORY1;
    AMP::shared_ptr<const VectorFactory> FACTORY2;
};


}
}

/// \endcond

#endif
