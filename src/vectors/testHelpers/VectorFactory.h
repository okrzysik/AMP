#ifndef included_AMP_test_VectorFactory
#define included_AMP_test_VectorFactory

#include "AMP/utils/AMP_MPI.h"
#include "AMP/vectors/ArrayVector.h"
#include "AMP/vectors/ManagedVector.h"
#include "AMP/vectors/MultiVariable.h"
#include "AMP/vectors/MultiVector.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/VectorBuilder.h"
#include "AMP/vectors/testHelpers/VectorTests.h"
#ifdef USE_PETSC
#include "AMP/vectors/trilinos/epetra/EpetraVectorEngine.h"
#include "petscvec.h"
#endif
#ifdef USE_EXT_TRILINOS
#include "AMP/vectors/trilinos/epetra/EpetraVectorEngine.h"
#endif


namespace AMP {
namespace LinearAlgebra {


class CloneFactory : public VectorFactory
{
public:
    explicit CloneFactory( std::shared_ptr<const VectorFactory> factory ) : d_factory( factory ) {}

    virtual AMP::LinearAlgebra::Variable::shared_ptr getVariable() const override
    {
        return d_factory->getVariable()->cloneVariable( "noname" );
    }

    virtual AMP::LinearAlgebra::Vector::shared_ptr getVector() const override
    {
        auto vec = d_factory->getVector();
        return vec->cloneVector();
    }

    virtual std::string name() const override { return "CloneFactory<" + d_factory->name() + ">"; }

    virtual AMP::Discretization::DOFManager::shared_ptr getDOFMap() const override
    {
        return d_factory->getDOFMap();
    }

private:
    CloneFactory();
    std::shared_ptr<const VectorFactory> d_factory;
};


class NullVectorFactory : public VectorFactory
{
public:
    typedef AMP::LinearAlgebra::Vector vector;

    virtual AMP::LinearAlgebra::Variable::shared_ptr getVariable() const override
    {
        return std::make_shared<AMP::LinearAlgebra::Variable>( "null" );
    }

    virtual AMP::LinearAlgebra::Vector::shared_ptr getVector() const override
    {
        return std::make_shared<AMP::LinearAlgebra::Vector>( "null" );
    }

    virtual std::string name() const override { return "NullVectorFactory"; }

    virtual AMP::Discretization::DOFManager::shared_ptr getDOFMap() const override
    {
        return AMP::Discretization::DOFManager::shared_ptr();
    }
};


template<class TYPE       = double,
         typename VecOps  = VectorOperationsDefault<TYPE>,
         typename VecData = VectorDataCPU<TYPE>>
class SimpleVectorFactory : public VectorFactory
{
public:
    SimpleVectorFactory( int i, bool global ) : I( i ), GLOBAL( global ) {}

    virtual AMP::LinearAlgebra::Variable::shared_ptr getVariable() const override
    {
        return std::make_shared<AMP::LinearAlgebra::Variable>( "simple" );
    }

    virtual AMP::LinearAlgebra::Vector::shared_ptr getVector() const override
    {
        AMP::LinearAlgebra::Vector::shared_ptr vec;
        if ( GLOBAL )
            vec = AMP::LinearAlgebra::createSimpleVector<TYPE, VecOps, VecData>(
                I, getVariable(), AMP_MPI( AMP_COMM_WORLD ) );
        else
            vec = AMP::LinearAlgebra::createSimpleVector<TYPE, VecOps, VecData>( I, getVariable() );
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


template<class TYPE = double>
class ArrayVectorFactory : public VectorFactory
{
public:
    ArrayVectorFactory( int d, int i, bool global ) : D( d ), I( i ), GLOBAL( global ) {}

    virtual AMP::LinearAlgebra::Variable::shared_ptr getVariable() const override
    {
        return std::make_shared<AMP::LinearAlgebra::Variable>( "array" );
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

template<typename TYPE>
class SimpleManagedVectorFactory : public VectorFactory
{
public:
    SimpleManagedVectorFactory() {}

    virtual AMP::LinearAlgebra::Variable::shared_ptr getVariable() const override
    {
        return std::make_shared<AMP::LinearAlgebra::Variable>( "..." );
    }

    virtual AMP::LinearAlgebra::Vector::shared_ptr getVector() const override
    {
        const int nLocal = 210;
        AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
        const int start   = nLocal * globalComm.getRank();
        const int nGlobal = nLocal * globalComm.getSize();
        auto commList   = AMP::LinearAlgebra::CommunicationList::createEmpty( nLocal, globalComm );
        auto dofManager = std::make_shared<AMP::Discretization::DOFManager>( nLocal, globalComm );
        auto epetraParams = std::make_shared<AMP::LinearAlgebra::EpetraVectorEngineParameters>(
            commList, dofManager );
        auto managedParams = std::make_shared<AMP::LinearAlgebra::ManagedVectorParameters>();
        managedParams->d_Buffer =
            std::make_shared<AMP::LinearAlgebra::VectorDataCPU<double>>( start, nLocal, nGlobal );
        managedParams->d_Engine = std::make_shared<AMP::LinearAlgebra::EpetraVectorEngine>(
            epetraParams, managedParams->d_Buffer );
        managedParams->d_CommList = commList;

        managedParams->d_DOFManager = dofManager;

        auto retval = std::make_shared<TYPE>( managedParams );
        retval->setVariable( std::make_shared<AMP::LinearAlgebra::Variable>( "Test Vector" ) );
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
class NativePetscVectorFactory : public VectorFactory
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
        VecSetType( v, VECMPI ); // this line will have to be modified for the no mpi and cuda cases
        auto newVec = createVector( v, true );
        VecSetFromOptions( v );
        newVec->getVectorData()->assemble();
        newVec->setVariable(
            std::make_shared<AMP::LinearAlgebra::Variable>( "Test NativePetscVector" ) );
        return newVec;
    }

    virtual std::string name() const override { return "NativePetscVectorFactory"; }

    virtual AMP::Discretization::DOFManager::shared_ptr getDOFMap() const override
    {
        return getVector()->getDOFManager();
    }
};

template<typename T>
class PetscManagedVectorFactory : public VectorFactory
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
        auto newVec = createVector( v, true );
        VecSetFromOptions( v );
        newVec->getVectorData()->assemble();
        auto p1        = std::make_shared<AMP::LinearAlgebra::ManagedVectorParameters>();
        p1->d_Engine   = newVec;
        p1->d_CommList = AMP::LinearAlgebra::CommunicationList::createEmpty( 210, globalComm );
        auto retval    = std::make_shared<T>( p1 );
        retval->setVariable( std::make_shared<AMP::LinearAlgebra::Variable>( "Test Vector" ) );
        return retval;
    }

    virtual std::string name() const override { return "PetscManagedVectorFactory"; }
};
#endif

#ifdef USE_EXT_SUNDIALS
class NativeSundialsVectorFactory : public VectorFactory
{
public:
    virtual AMP::LinearAlgebra::Variable::shared_ptr getVariable() const override
    {
        return AMP::LinearAlgebra::Variable::shared_ptr(); // no variable.....
    }

    virtual AMP::LinearAlgebra::Vector::shared_ptr getVector() const override
    {
        AMP::LinearAlgebra::Vector::shared_ptr newVec;
        AMP_ERROR( "Not implemented" );
        return newVec;
    }

    virtual std::string name() const override { return "NativeSundialsVectorFactory"; }

    virtual AMP::Discretization::DOFManager::shared_ptr getDOFMap() const override
    {
        return getVector()->getDOFManager();
    }
};
#endif

template<typename TYPE>
class ViewFactory : public VectorFactory
{
public:
    explicit ViewFactory( std::shared_ptr<const VectorFactory> factory ) : d_factory( factory ) {}

    virtual AMP::LinearAlgebra::Variable::shared_ptr getVariable() const override
    {
        return d_factory->getVariable();
    }

    virtual AMP::LinearAlgebra::Vector::shared_ptr getVector() const override
    {
        auto rtn = TYPE::view( d_factory->getVector() );
        auto vec = dynamic_cast<TYPE *>( rtn.get() );
        AMP_INSIST( vec != nullptr, "Failed to cast view to type" );
        auto native = vec->getNativeVec();
        NULL_USE( native );
        return rtn;
    }

    virtual std::string name() const override { return "ViewFactory<" + d_factory->name() + ">"; }

    virtual AMP::Discretization::DOFManager::shared_ptr getDOFMap() const override
    {
        return d_factory->getDOFMap();
    }

private:
    ViewFactory();
    std::shared_ptr<const VectorFactory> d_factory;
};


class MultiVectorFactory : public VectorFactory
{
public:
    MultiVectorFactory( std::shared_ptr<const VectorFactory> factory1,
                        int N1,
                        std::shared_ptr<const VectorFactory> factory2,
                        int N2 )
        : NUM1( N1 ), NUM2( N2 ), FACTORY1( factory1 ), FACTORY2( factory2 )
    {
    }

    virtual AMP::LinearAlgebra::Variable::shared_ptr getVariable() const override
    {
        auto newVar = std::make_shared<AMP::LinearAlgebra::MultiVariable>( "var1" );
        for ( int i = 0; i != NUM1; i++ )
            newVar->add( FACTORY1->getVariable() );
        for ( int i = 0; i != NUM2; i++ )
            newVar->add( FACTORY2->getVariable() );
        return newVar;
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
    std::shared_ptr<const VectorFactory> FACTORY1;
    std::shared_ptr<const VectorFactory> FACTORY2;
};


} // namespace LinearAlgebra
} // namespace AMP

/// \endcond

#endif
