
#ifndef included_AMP_GaussPointToGaussPointMap
#define included_AMP_GaussPointToGaussPointMap

#include "operators/map/NodeToNodeMap.h"

#include <vector>

namespace AMP {
  namespace Operator {

    class GaussPointToGaussPointMap : public NodeToNodeMap {
      public:
        GaussPointToGaussPointMap(const boost::shared_ptr<AMP::Operator::OperatorParameters> & params)
          : NodeToNodeMap(params) {
            createIdxMap(params);
            d_useFrozenInputVec = params->d_db->getBoolWithDefault("FrozenInput", false);
          }

        void applyStart(const AMP::LinearAlgebra::Vector::shared_ptr &f, const AMP::LinearAlgebra::Vector::shared_ptr &u,
            AMP::LinearAlgebra::Vector::shared_ptr &r, const double a = -1.0, const double b = 1.0) {
          AMP::LinearAlgebra::Vector::shared_ptr uInternal = u;
          if(d_useFrozenInputVec) {
            uInternal = d_frozenInputVec;
          }
          AMP::Operator::NodeToNodeMap::applyStart(f, uInternal, r, a, b);
        }

        void applyFinish(const AMP::LinearAlgebra::Vector::shared_ptr &f, const AMP::LinearAlgebra::Vector::shared_ptr &u,
            AMP::LinearAlgebra::Vector::shared_ptr &r, const double a = -1.0, const double b = 1.0) {
          AMP::Operator::NodeToNodeMap::applyFinish(f, u, r, a, b);
          correctLocalOrdering();
        }

        static bool  validMapType ( const std::string &t ) {
          if ( t == "GaussPointToGaussPoint" ) {
            return true;
          }
          return false;
        }

        void setFrozenInputVector(AMP::LinearAlgebra::Vector::shared_ptr u) {
          d_frozenInputVec = u;
        }

        virtual ~GaussPointToGaussPointMap() { }

      protected:
        void createIdxMap(boost::shared_ptr<AMP::Operator::OperatorParameters> params);

        void correctLocalOrdering();

        bool d_useFrozenInputVec;

        AMP::LinearAlgebra::Vector::shared_ptr d_frozenInputVec; 

        std::vector<std::vector<unsigned int> > d_idxMap;
    };

  }
}

#endif


