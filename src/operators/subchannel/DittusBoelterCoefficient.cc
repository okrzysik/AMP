#include "DittusBoelterCoefficient.h"
#include "utils/Utilities.h"
#include <cmath>


namespace AMP {
  namespace Operator {


    DittusBoelterCoefficient :: DittusBoelterCoefficient(const
        boost::shared_ptr<RobinPhysicsModelParameters>& params) : RobinPhysicsModel (params) {

      d_k  = params->d_db->getDouble("ThermalConductivity"); 
      d_De = params->d_db->getDouble("ChannelDiameter");
      d_Re = params->d_db->getDouble("ReynoldsNumber");
      d_Pr = params->d_db->getDouble("PrandtlNumber");

    }

    void DittusBoelterCoefficient :: getConductance(std::vector<double> & beta,
        std::vector<double> & gamma, const std::vector<std::vector <double> > & inputVectors)
    {

      for (size_t l=1; l<inputVectors[0].size(); l++)
      {
        beta[l]  =  (0.023*d_k/d_De) * pow(d_Re, 0.8) * pow(d_Pr, 0.4);
        gamma[l] =  (0.023*d_k/d_De) * pow(d_Re, 0.8) * pow(d_Pr, 0.4);
      }

    }

  }
}


