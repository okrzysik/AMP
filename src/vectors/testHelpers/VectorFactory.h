#ifndef included_AMP_test_VectorFactory
#define included_AMP_test_VectorFactory

#include "AMP/utils/AMP_MPI.h"
#include "AMP/vectors/MultiVariable.h"
#include "AMP/vectors/MultiVector.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/VectorBuilder.h"
#include "AMP/vectors/testHelpers/VectorTests.h"


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
    SimpleVectorFactory( int i, bool global, std::string name = "SimpleVectorFactory" )
        : I( i ), GLOBAL( global ), NAME( std::move( name ) )
    {
    }

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

    virtual std::string name() const override { return NAME; }

    virtual AMP::Discretization::DOFManager::shared_ptr getDOFMap() const override
    {
        return getVector()->getDOFManager();
    }

private:
    SimpleVectorFactory();
    int I;
    bool GLOBAL;
    std::string NAME;
};


template<class TYPE = double>
class ArrayVectorFactory : public VectorFactory
{
public:
    ArrayVectorFactory( size_t d, size_t i, bool global ) : D( d ), I( i ), GLOBAL( global ) {}

    virtual AMP::LinearAlgebra::Variable::shared_ptr getVariable() const override
    {
        return std::make_shared<AMP::LinearAlgebra::Variable>( "array" );
    }

    virtual AMP::LinearAlgebra::Vector::shared_ptr getVector() const override
    {
        AMP::LinearAlgebra::Vector::shared_ptr vec;
        if ( GLOBAL )
            vec = AMP::LinearAlgebra::createArrayVector<TYPE>(
                { D, I }, getVariable(), AMP_MPI( AMP_COMM_WORLD ) );
        else
            vec = AMP::LinearAlgebra::createArrayVector<TYPE>( { D, I }, getVariable() );
        return vec;
    }

    virtual std::string name() const override { return "ArrayVectorFactory"; }

    virtual AMP::Discretization::DOFManager::shared_ptr getDOFMap() const override
    {
        return getVector()->getDOFManager();
    }

private:
    ArrayVectorFactory();
    size_t D, I;
    bool GLOBAL;
};


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
        auto vec = TYPE::view( d_factory->getVector() );
        AMP_INSIST( vec != nullptr, "Failed to cast view to type" );
        auto native = vec->getNativeVec();
        NULL_USE( native );
        return vec->getManagedVec();
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

    virtual std::string name() const override
    {
        return "MultiVectorFactory<" + FACTORY1->name() + "," + std::to_string( NUM1 ) + "," +
               FACTORY2->name() + "," + std::to_string( NUM2 ) + ">";
    }

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
