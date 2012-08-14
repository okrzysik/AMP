
#ifndef included_AMP_MeshBasedOperator
#define included_AMP_MeshBasedOperator

#include "operators/Operator.h"
#include "operators/MeshBasedOperatorParameters.h"

namespace AMP {
  namespace Operator {

    class MeshBasedOperator : public Operator {
      public :
        MeshBasedOperator(const boost::shared_ptr<MeshBasedOperatorParameters> & params)
          : Operator(params) {
            d_Mesh = params->d_Mesh;
          }

        virtual ~MeshBasedOperator() { }

        AMP::Mesh::Mesh::shared_ptr getMesh() {
          return d_Mesh;
        }

      protected:
        AMP::Mesh::Mesh::shared_ptr d_Mesh;
    };

  }
}

#endif


