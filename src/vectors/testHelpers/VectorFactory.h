#ifndef included_AMP_test_VectorFactory
#define included_AMP_test_VectorFactory

#include "AMP/utils/AMP_MPI.h"
#include "AMP/vectors/MultiVariable.h"
#include "AMP/vectors/MultiVector.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/VectorBuilder.h"
#include "AMP/vectors/VectorBuilder.hpp"
#include "AMP/vectors/testHelpers/VectorTests.h"


namespace AMP::LinearAlgebra {


class CloneFactory : public VectorFactory
{
public:
    explicit CloneFactory( std::shared_ptr<const VectorFactory> factory ) : d_factory( factory ) {}
    AMP::LinearAlgebra::Vector::shared_ptr getVector() const override
    {
        auto vec = d_factory->getVector();
        return vec->clone();
    }
    std::string name() const override { return "CloneFactory<" + d_factory->name() + ">"; }

private:
    CloneFactory();
    std::shared_ptr<const VectorFactory> d_factory;
};


class NullVectorFactory : public VectorFactory
{
public:
    AMP::LinearAlgebra::Vector::shared_ptr getVector() const override
    {
        return std::make_shared<AMP::LinearAlgebra::Vector>( "null" );
    }
    std::string name() const override { return "NullVectorFactory"; }
};


template<class TYPE       = double,
         typename VecOps  = VectorOperationsDefault<TYPE>,
         typename VecData = VectorDataDefault<TYPE>>
class SimpleVectorFactory : public VectorFactory
{
public:
    SimpleVectorFactory( int i, bool global, std::string name = "SimpleVectorFactory" )
        : I( i ), GLOBAL( global ), NAME( std::move( name ) )
    {
    }
    AMP::LinearAlgebra::Vector::shared_ptr getVector() const override
    {
        AMP::LinearAlgebra::Vector::shared_ptr vec;
        auto var = std::make_shared<AMP::LinearAlgebra::Variable>( "simple" );
        if ( GLOBAL )
            vec = AMP::LinearAlgebra::createSimpleVector<TYPE, VecOps, VecData>(
                I, var, AMP_MPI( AMP_COMM_WORLD ) );
        else
            vec = AMP::LinearAlgebra::createSimpleVector<TYPE, VecOps, VecData>( I, var );
        return vec;
    }
    std::string name() const override { return NAME; }

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
    AMP::LinearAlgebra::Vector::shared_ptr getVector() const override
    {
        AMP::LinearAlgebra::Vector::shared_ptr vec;
        auto var = std::make_shared<AMP::LinearAlgebra::Variable>( "array" );
        if ( GLOBAL ) {
            AMP_MPI comm( AMP_COMM_WORLD );
            ArraySize index( comm.getRank(), 0, 0 );
            vec = AMP::LinearAlgebra::createArrayVector<TYPE>( { D, I }, index, comm, var );
        } else {
            vec = AMP::LinearAlgebra::createArrayVector<TYPE>( { D, I }, var );
        }
        return vec;
    }
    std::string name() const override { return "ArrayVectorFactory"; }

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
    AMP::LinearAlgebra::Vector::shared_ptr getVector() const override
    {
        auto vec = TYPE::view( d_factory->getVector() );
        AMP_INSIST( vec != nullptr, "Failed to cast view to type" );
        auto native = vec->getNativeVec();
        NULL_USE( native );
        return vec->getManagedVec();
    }
    std::string name() const override { return "ViewFactory<" + d_factory->name() + ">"; }

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
    AMP::LinearAlgebra::Vector::shared_ptr getVector() const override
    {
        AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
        auto var    = std::make_shared<AMP::LinearAlgebra::MultiVariable>( "var1" );
        auto retVal = AMP::LinearAlgebra::MultiVector::create( var, globalComm );
        for ( int i = 0; i != NUM1; i++ )
            retVal->addVector( FACTORY1->getVector() );
        for ( int i = 0; i != NUM2; i++ )
            retVal->addVector( FACTORY2->getVector() );
        return retVal;
    }
    std::string name() const override
    {
        return "MultiVectorFactory<" + FACTORY1->name() + "," + std::to_string( NUM1 ) + "," +
               FACTORY2->name() + "," + std::to_string( NUM2 ) + ">";
    }

private:
    MultiVectorFactory();
    int NUM1, NUM2;
    std::shared_ptr<const VectorFactory> FACTORY1;
    std::shared_ptr<const VectorFactory> FACTORY2;
};


} // namespace AMP::LinearAlgebra

/// \endcond

#endif
