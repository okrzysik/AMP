
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

        virtual ~GaussPointToGaussPointMap() { }

      protected:
        void createIdxMap(boost::shared_ptr<AMP::Operator::OperatorParameters> params);

        void correctLocalOrdering();

        std::vector<std::vector<unsigned int> > d_idxMap;
    };

  }
}

#endif


