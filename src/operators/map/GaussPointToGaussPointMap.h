
#ifndef included_AMP_GaussPointToGaussPointMap
#define included_AMP_GaussPointToGaussPointMap

#include "AMP/operators/map/NodeToNodeMap.h"

#include <vector>

namespace AMP {
namespace Operator {

class GaussPointToGaussPointMap : public NodeToNodeMap
{
public:
    explicit GaussPointToGaussPointMap(
        const AMP::shared_ptr<AMP::Operator::OperatorParameters> &params );

    virtual void applyStart( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                             AMP::LinearAlgebra::Vector::shared_ptr f ) override;

    virtual void applyFinish( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                              AMP::LinearAlgebra::Vector::shared_ptr f ) override;

    static bool validMapType( const std::string &t );

    void setFrozenInputVector( AMP::LinearAlgebra::Vector::shared_ptr u ) { d_frozenInputVec = u; }

    virtual ~GaussPointToGaussPointMap() {}

protected:
    void createIdxMap( AMP::shared_ptr<AMP::Operator::OperatorParameters> params );

    void correctLocalOrdering();

    bool d_useFrozenInputVec;

    AMP::LinearAlgebra::Vector::shared_ptr d_frozenInputVec;

    std::vector<std::vector<unsigned int>> d_idxMap;
};
} // namespace Operator
} // namespace AMP

#endif
