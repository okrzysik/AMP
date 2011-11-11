#ifndef included_AMP_SourcePhysicsModel
#define included_AMP_SourcePhysicsModel

#include <cstring>
#include "boost/shared_ptr.hpp"

#include "ElementPhysicsModel.h"
#include "materials/Material.h"
#include "ElementPhysicsModelParameters.h"
#include "MassDensityModel.h"


namespace AMP {
namespace Operator {


typedef ElementPhysicsModelParameters SourcePhysicsModelParameters;


/**
  This class acts a interface to various ElementPhysicsModel and also 
  for source fuction evaluation for the VolumeIntegral Operator. The 
  Element Operation of this operator passes the input vectors and coordinates 
  through the getConstitutiveProperty function for calculating the source terms.
*/
class SourcePhysicsModel : public ElementPhysicsModel
{
public :

      /**
        This Constructor reads the value for the key USE_MATERIALS_LIBRARY (false by default)
        and also USE_ELEMENT_PHYSICS from the database. These specify if material library or the 
        element physics model are used for calculating the source terms.
        */
      SourcePhysicsModel (const boost::shared_ptr<SourcePhysicsModelParameters>& params );

     /**
      * Destructor.
      */ 
     virtual ~SourcePhysicsModel() { }

     void getConstitutiveProperty( std::vector<double>                & result, 
                                   std::vector< std::vector<double> > & InputVec, 
                                   std::vector< std::vector<double> > & ,
                                   const std::vector<Point> & Coordinates  );

protected :

     bool d_useMaterialsLibrary;

     std::string d_physicsName;

     AMP::Materials::Material::shared_ptr d_material;

     AMP::Materials::PropertyPtr d_property;

private :

     double d_DefaultTemperature;
     double d_DefaultConcentration;
     double d_DefaultBurnup;

     std::vector<double> d_defaults;

     boost::shared_ptr<ElementPhysicsModel> elementPhysicsModel; 

     boost::shared_ptr<ElementPhysicsModelParameters> elementPhysicsParams; 

  };

}
}

#endif
