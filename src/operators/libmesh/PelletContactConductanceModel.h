#ifndef included_AMP_PelletContactConductanceModel
#define included_AMP_PelletContactConductanceModel

#include <iostream>
#include <string>
#include <vector>

#include "operators/ElementPhysicsModelFactory.h"
#include "operators/diffusion/DiffusionTransportModel.h"
#include "operators/boundary/libmesh/RobinPhysicsModel.h"
#include "utils/shared_ptr.h"

namespace AMP {
namespace Operator {

class PelletContactConductanceModel  : public RobinPhysicsModel
 {
    public :
      
      explicit PelletContactConductanceModel(const AMP::shared_ptr<RobinPhysicsModelParameters>& params);

      virtual ~PelletContactConductanceModel() {}

      void getConductance(std::vector<double> & beta,
          std::vector<double> & gamma, 
          const std::vector<std::vector <double> > &arguments);

    protected:

      unsigned int d_nTransportModels; /**< Number of Transport Models. */
      std::vector< AMP::shared_ptr< DiffusionTransportModel > > d_transportModels;

    private:

 };

}
}

#endif

