#ifndef included_AMP_PelletContactConductanceModel
#define included_AMP_PelletContactConductanceModel

#include <iostream>
#include <string>
#include <vector>

#include "operators/ElementPhysicsModelFactory.h"
#include "operators/diffusion/DiffusionTransportModel.h"
#include "operators/boundary/RobinPhysicsModel.h"
#include "boost/shared_ptr.hpp"

namespace AMP {
namespace Operator {

class PelletContactConductanceModel  : public RobinPhysicsModel
 {
    public :
      
      PelletContactConductanceModel(const boost::shared_ptr<RobinPhysicsModelParameters>& params);

      virtual ~PelletContactConductanceModel() {}

      void getConductance(std::vector<double> & beta,
          std::vector<double> & gamma, 
          const std::vector<std::vector <double> > &arguments);

    protected:

      unsigned int d_nTransportModels; /**< Number of Transport Models. */
      std::vector< boost::shared_ptr< DiffusionTransportModel > > d_transportModels;

    private:

 };

}
}

#endif

