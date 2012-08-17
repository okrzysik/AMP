#ifndef included_AMP_ConvctiveHeatCoefficient
#define included_AMP_ConvctiveHeatCoefficient

#include <iostream>
#include <string>
#include <vector>

#include "ElementPhysicsModelFactory.h"
#include "RobinPhysicsModel.h"
#include "boost/shared_ptr.hpp"

namespace AMP {
namespace Operator {

class ConvctiveHeatCoefficient  : public RobinPhysicsModel
 {
    public :
      
      ConvctiveHeatCoefficient(const boost::shared_ptr<RobinPhysicsModelParameters>& params);

      virtual ~ConvctiveHeatCoefficient() {}

      void getConductance(std::vector<double> & beta,
          std::vector<double> & gamma, 
          const std::vector<std::vector <double> > &arguments);

      AMP::Materials::Material::shared_ptr getMaterial(){return d_material;}
      AMP::Materials::PropertyPtr getProperty(){return d_property;}

    protected:

      AMP::Materials::Material::shared_ptr d_material;

      AMP::Materials::PropertyPtr d_property;

      double d_k , d_De , d_Re , d_Pr;

      std::vector<std::string> d_argNames;
    
      std::vector<double> d_defaults;

    private:

 };

}
}

#endif

