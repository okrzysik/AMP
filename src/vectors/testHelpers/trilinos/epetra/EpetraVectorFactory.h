#ifndef included_AMP_EpetraVectorFactor
#define included_AMP_EpetraVectorFactor


#include "AMP/utils/UnitTest.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/trilinos/epetra/ManagedEpetraVector.h"
#include "AMP/vectors/testHelpers/VectorTests.h"


/// \cond UNDOCUMENTED


namespace AMP {
namespace LinearAlgebra {


class NativeEpetraFactory : public VectorFactory
{
public:
    NativeEpetraFactory() {}

    AMP::LinearAlgebra::Variable::shared_ptr getVariable() const override
    {
        return std::make_shared<AMP::LinearAlgebra::Variable>( "epetra" );
    }

    AMP::LinearAlgebra::Vector::shared_ptr getVector() const override
    {
        const int nLocal = 210;
        AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
        const int start   = nLocal * globalComm.getRank();
        const int nGlobal = nLocal * globalComm.getSize();
        auto commList   = AMP::LinearAlgebra::CommunicationList::createEmpty( nLocal, globalComm );
        auto dofManager = std::make_shared<AMP::Discretization::DOFManager>( nLocal, globalComm );
        auto buffer =
            std::make_shared<AMP::LinearAlgebra::VectorDataCPU<double>>( start, nLocal, nGlobal );
        auto vec = createEpetraVector( commList, dofManager, buffer );
        return vec;
    }

    std::string name() const override { return "NativeEpetraFactory"; }

    AMP::Discretization::DOFManager::shared_ptr getDOFMap() const override
    {
        return getVector()->getDOFManager();
    }
};


class ManagedEpetraVectorFactory : public VectorFactory
{
public:
    ManagedEpetraVectorFactory( std::shared_ptr<const VectorFactory> factory ) : d_factory( factory ) {}

    virtual AMP::LinearAlgebra::Variable::shared_ptr getVariable() const override
    {
        return std::make_shared<AMP::LinearAlgebra::Variable>( "..." );
    }

    virtual AMP::LinearAlgebra::Vector::shared_ptr getVector() const override
    {
        auto engine = d_factory->getVector();
        auto retval = std::make_shared<ManagedEpetraVector>( engine );
        retval->setVariable( std::make_shared<AMP::LinearAlgebra::Variable>( "Test Vector" ) );
        return retval;
    }

    virtual std::string name() const override { return "ManagedEpetraVectorFactory<" + d_factory->name() + ">"; }

    virtual AMP::Discretization::DOFManager::shared_ptr getDOFMap() const override
    {
        return getVector()->getDOFManager();
    }
private:
    std::shared_ptr<const VectorFactory> d_factory;
};



} // namespace LinearAlgebra
} // namespace AMP

/// \endcond

#endif
