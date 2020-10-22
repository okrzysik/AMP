#ifndef included_AMP_SundialsVectorFactory
#define included_AMP_SundialsVectorFactory

#include "AMP/utils/AMP_MPI.h"
#include "AMP/vectors/testHelpers/VectorTests.h"


namespace AMP {
namespace LinearAlgebra {


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


} // namespace LinearAlgebra
} // namespace AMP

/// \endcond

#endif
