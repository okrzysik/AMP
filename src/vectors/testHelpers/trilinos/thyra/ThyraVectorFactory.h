#ifndef included_ThyraVectorFactor
#define included_ThyraVectorFactor


#include "utils/UnitTest.h"
#include "vectors/Vector.h"
#include "vectors/testHelpers/VectorTests.h"
#include <vectors/trilinos/thyra/ManagedThyraVector.h>
#include <vectors/trilinos/thyra/NativeThyraVector.h>


/// \cond UNDOCUMENTED


namespace AMP {
namespace LinearAlgebra {


class NativeThyraFactory: public VectorFactory
{
public:
    
    NativeThyraFactory() {}

    virtual AMP::LinearAlgebra::Variable::shared_ptr getVariable() const override;

    virtual AMP::LinearAlgebra::Vector::shared_ptr getVector() const override;

    virtual std::string name() const override { return "NativeThyraFactory"; }

    virtual AMP::Discretization::DOFManager::shared_ptr getDOFMap() const override;
};


class ManagedThyraFactory: public VectorFactory
{
public:
    
    ManagedThyraFactory( AMP::shared_ptr<VectorFactory> factory ): d_factory(factory) {}

    virtual AMP::LinearAlgebra::Variable::shared_ptr getVariable() const override;

    virtual AMP::LinearAlgebra::Vector::shared_ptr getVector() const override;

    virtual std::string name() const override { return "ManagedThyraFactory"; }

    virtual AMP::Discretization::DOFManager::shared_ptr getDOFMap() const override;

private:
    ManagedThyraFactory();
    AMP::shared_ptr<VectorFactory> d_factory;
};


class ManagedNativeThyraFactory: public VectorFactory
{
public:

    ManagedNativeThyraFactory( AMP::shared_ptr<VectorFactory> factory ): d_factory(factory) {}

    virtual AMP::LinearAlgebra::Variable::shared_ptr getVariable() const override;

    virtual AMP::LinearAlgebra::Vector::shared_ptr getVector() const override;

    virtual std::string name() const override { return "ManagedNativeThyraFactory"; }

    virtual AMP::Discretization::DOFManager::shared_ptr getDOFMap() const override;

private:
    ManagedNativeThyraFactory();
    AMP::shared_ptr<VectorFactory> d_factory;
};

#ifdef USE_TRILINOS_BELOS
void testBelosThyraVector( AMP::UnitTest &utils, const VectorFactory& factory );
#endif

}
}

/// \endcond

#endif
