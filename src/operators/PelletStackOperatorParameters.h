
#ifndef included_AMP_PelletStackOperatorParameters
#define included_AMP_PelletStackOperatorParameters

#include "operators/map/AsyncMapColumnOperator.h"

namespace AMP {
  namespace Operator {

    class PelletStackOperatorParameters : public OperatorParameters {
      public :

        PelletStackOperatorParameters(const boost::shared_ptr<AMP::Database>& db)
          : OperatorParameters(db) {
            d_currentPellet = static_cast<unsigned int>(-1);
          }

        virtual ~PelletStackOperatorParameters() { }

        unsigned int d_currentPellet;
        AMP_MPI d_pelletStackComm;
        boost::shared_ptr<AMP::Operator::AsyncMapColumnOperator> d_n2nMaps;
        AMP::Mesh::MeshManager::shared_ptr d_meshManager;
    };

  }  
}


#endif



