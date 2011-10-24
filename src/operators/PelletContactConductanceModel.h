#ifndef included_AMP_PelletContactConductanceModel
#define included_AMP_PelletContactConductanceModel

#include <iostream>
#include <string>
#include <vector>

#include "ElementPhysicsModelFactory.h"
#include "DiffusionTransportModel.h"
#include "RobinPhysicsModel.h"
#include "boost/shared_ptr.hpp"

namespace AMP {
namespace Operator {

class PelletContactConductanceModel  : public RobinPhysicsModel
 {
    public :
      
      PelletContactConductanceModel(const
          boost::shared_ptr<RobinPhysicsModelParameters>& params) : RobinPhysicsModel (params) {
      
          d_nTransportModels   =  (params->d_db)->getInteger("Number_TransportModels");
          d_transportModels.resize( d_nTransportModels );
          boost::shared_ptr<ElementPhysicsModel> elementPhysicsModel;
          
          for (unsigned int i = 0; i < d_nTransportModels ; i++)
          {
            char key[100];
            sprintf(key, "DiffusionTransportModel_%d", (int)i);
            boost::shared_ptr<Database> transportModel_db = (params->d_db)->getDatabase(key);
            elementPhysicsModel = ElementPhysicsModelFactory::createElementPhysicsModel(transportModel_db);
            d_transportModels[i] = boost::dynamic_pointer_cast<DiffusionTransportModel>(elementPhysicsModel) ;
          }

      }

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

