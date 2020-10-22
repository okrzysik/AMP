#ifndef included_AMP_ThyraVectorFactor
#define included_AMP_ThyraVectorFactor


#include "AMP/utils/UnitTest.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/testHelpers/VectorTests.h"


/// \cond UNDOCUMENTED


namespace AMP {
namespace LinearAlgebra {


class NativeThyraFactory : public VectorFactory
{
public:
    NativeThyraFactory() {}

    virtual AMP::LinearAlgebra::Variable::shared_ptr getVariable() const override;

    virtual AMP::LinearAlgebra::Vector::shared_ptr getVector() const override;

    virtual std::string name() const override { return "NativeThyraFactory"; }

    virtual AMP::Discretization::DOFManager::shared_ptr getDOFMap() const override;
};


class ManagedThyraFactory : public VectorFactory
{
public:
    explicit ManagedThyraFactory( std::shared_ptr<VectorFactory> factory ) : d_factory( factory ) {}

    virtual AMP::LinearAlgebra::Variable::shared_ptr getVariable() const override;

    virtual AMP::LinearAlgebra::Vector::shared_ptr getVector() const override;

    virtual std::string name() const override
    {
        return "ManagedThyraFactory<" + d_factory->name() + ">";
    }

    virtual AMP::Discretization::DOFManager::shared_ptr getDOFMap() const override;

private:
    ManagedThyraFactory();
    std::shared_ptr<VectorFactory> d_factory;
};


class ManagedNativeThyraFactory : public VectorFactory
{
public:
    explicit ManagedNativeThyraFactory( std::shared_ptr<VectorFactory> factory )
        : d_factory( factory )
    {
    }

    virtual AMP::LinearAlgebra::Variable::shared_ptr getVariable() const override;

    virtual AMP::LinearAlgebra::Vector::shared_ptr getVector() const override;

    virtual std::string name() const override
    {
        return "ManagedNativeThyraFactory<" + d_factory->name() + ">";
    }

    virtual AMP::Discretization::DOFManager::shared_ptr getDOFMap() const override;

private:
    ManagedNativeThyraFactory();
    std::shared_ptr<VectorFactory> d_factory;
};

#ifdef USE_TRILINOS_BELOS
void testBelosThyraVector( AMP::UnitTest &utils, const VectorFactory &factory );
#endif
} // namespace LinearAlgebra
} // namespace AMP

/// \endcond

#endif
