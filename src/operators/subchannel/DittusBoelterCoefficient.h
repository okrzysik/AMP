#ifndef included_AMP_DittusBoelterCoefficient
#define included_AMP_DittusBoelterCoefficient

#include <iostream>
#include <string>
#include <vector>

#include "ElementPhysicsModelFactory.h"
#include "RobinPhysicsModel.h"
#include "boost/shared_ptr.hpp"

namespace AMP {
namespace Operator {

class DittusBoelterCoefficient  : public RobinPhysicsModel
 {
    public :
      
      DittusBoelterCoefficient(const boost::shared_ptr<RobinPhysicsModelParameters>& params);

      virtual ~DittusBoelterCoefficient() {}

      void getConductance(std::vector<double> & beta,
          std::vector<double> & gamma, 
          const std::vector<std::vector <double> > &arguments);

    protected:

      double d_k , d_De , d_Re , d_Pr;



    private:

 };

}
}

#endif

