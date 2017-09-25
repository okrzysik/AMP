#ifndef included_test_StridedVectorFactory
#define included_test_StridedVectorFactory

#include "vectors/MultiVector.h"
#include "vectors/VectorSelector.h"
#include "vectors/testHelpers/VectorTests.h"

#include "utils/AMPManager.h"
#include "utils/UnitTest.h"


namespace AMP {
namespace LinearAlgebra {


class StridedVectorFactory : public VectorFactory
{
public:
    explicit StridedVectorFactory( AMP::shared_ptr<const VectorFactory> factory )
        : d_factory( factory )
    {
    }

    virtual AMP::LinearAlgebra::Variable::shared_ptr getVariable() const override
    {
        return AMP::LinearAlgebra::Variable::shared_ptr();
    }

    virtual AMP::LinearAlgebra::Vector::shared_ptr getVector() const override
    {
        auto vec        = d_factory->getVector();
        auto criterion  = AMP::LinearAlgebra::VS_Stride( 1, 3 );
        auto vec_select = vec->select( criterion, "thirds" );
        size_t N1       = vec->getGlobalSize();
        size_t N2       = vec_select->getGlobalSize();
        AMP_ASSERT( N1 / 3 == N2 );
        return vec_select;
    }

    virtual std::string name() const override { return "StridedVectorFactory"; }

    virtual AMP::Discretization::DOFManager::shared_ptr getDOFMap() const override
    {
        return d_factory->getDOFMap();
    }

private:
    AMP::shared_ptr<const VectorFactory> d_factory;
};
} // namespace LinearAlgebra
} // namespace AMP

#endif
