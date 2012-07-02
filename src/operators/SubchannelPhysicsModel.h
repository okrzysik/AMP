#ifndef included_AMP_SubchannelPhysicsModel
#define included_AMP_SubchannelPhysicsModel

#include <cstring>
#include "boost/shared_ptr.hpp"
#include <map>
#include <vector>

#include "ElementPhysicsModel.h"
#include "materials/Material.h"
#include "ElementPhysicsModelParameters.h"


namespace AMP {
namespace Operator {


typedef ElementPhysicsModelParameters SubchannelPhysicsModelParameters;


/**
  This class acts a interface to the water material library for the subchannel
  flow operators.
*/
class SubchannelPhysicsModel : public ElementPhysicsModel
{
public :

      SubchannelPhysicsModel (const boost::shared_ptr<SubchannelPhysicsModelParameters>& params );

     virtual ~SubchannelPhysicsModel() { }

     void getProperty(    std::string property,
                          std::vector<double>                & result, 
                          std::map<std::string, boost::shared_ptr<std::vector<double> > > & args);

protected :

     bool d_useMaterialsLibrary;

     std::string d_physicsName;

     AMP::Materials::Material::shared_ptr d_material;

     std::map<std::string,AMP::Materials::PropertyPtr> d_properties;

private :

     boost::shared_ptr<ElementPhysicsModel> elementPhysicsModel; 

     boost::shared_ptr<ElementPhysicsModelParameters> elementPhysicsParams; 

  };

}
}

#endif
