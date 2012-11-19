#ifndef included_AMP_ConvectiveHeatCoefficient
#define included_AMP_ConvectiveHeatCoefficient

#include <iostream>
#include <string>
#include <vector>

#include "ElementPhysicsModelFactory.h"
#include "RobinPhysicsModel.h"
#include "boost/shared_ptr.hpp"

namespace AMP {
namespace Operator {

class ConvectiveHeatCoefficient  : public RobinPhysicsModel
 {
    public :
      
      ConvectiveHeatCoefficient(const boost::shared_ptr<RobinPhysicsModelParameters>& params);

      virtual ~ConvectiveHeatCoefficient() {}

      void getConductance(std::vector<double> & beta,
          std::vector<double> & gamma, 
          const std::vector<std::vector <double> > &arguments);

      AMP::Materials::Material::shared_ptr getMaterial(){return d_material;}
      AMP::Materials::PropertyPtr getProperty(){return d_property;}

    protected:

      AMP::Materials::Material::shared_ptr d_material;

      AMP::Materials::PropertyPtr d_property;

      std::vector<std::string> d_argNames;
    
      std::vector<double> d_defaults;

    private:

 };

}
}

#endif

