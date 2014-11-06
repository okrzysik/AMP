
#ifndef included_AMP_GaussPointToGaussPointMap
#define included_AMP_GaussPointToGaussPointMap

#include "operators/map/NodeToNodeMap.h"

#include <vector>

namespace AMP {
namespace Operator {

class GaussPointToGaussPointMap : public NodeToNodeMap {
public:

    GaussPointToGaussPointMap(const AMP::shared_ptr<AMP::Operator::OperatorParameters> & params);

    virtual void applyStart(AMP::LinearAlgebra::Vector::const_shared_ptr f, AMP::LinearAlgebra::Vector::const_shared_ptr u,
            AMP::LinearAlgebra::Vector::shared_ptr r, const double a = -1.0, const double b = 1.0);

    virtual void applyFinish(AMP::LinearAlgebra::Vector::const_shared_ptr f, AMP::LinearAlgebra::Vector::const_shared_ptr u,
            AMP::LinearAlgebra::Vector::shared_ptr r, const double a = -1.0, const double b = 1.0);

    static bool  validMapType ( const std::string &t );

    void setFrozenInputVector(AMP::LinearAlgebra::Vector::shared_ptr u) {
          d_frozenInputVec = u;
    }

    virtual ~GaussPointToGaussPointMap() { }

protected:
    void createIdxMap(AMP::shared_ptr<AMP::Operator::OperatorParameters> params);

    void correctLocalOrdering();

    bool d_useFrozenInputVec;

    AMP::LinearAlgebra::Vector::shared_ptr d_frozenInputVec; 

    std::vector<std::vector<unsigned int> > d_idxMap;
};


}
}

#endif


