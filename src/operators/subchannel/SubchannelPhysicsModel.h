#ifndef included_AMP_SubchannelPhysicsModel
#define included_AMP_SubchannelPhysicsModel

#include "AMP/materials/Material.h"
#include "AMP/operators/ElementPhysicsModel.h"
#include "AMP/operators/ElementPhysicsModelParameters.h"
#include <memory>

#include <cstring>
#include <map>
#include <vector>


namespace AMP {
namespace Operator {


/**
  This class acts a interface to the water material library for the subchannel
  flow operators.
*/
class SubchannelPhysicsModel : public ElementPhysicsModel
{
public:
    /**
      Constructor
      */
    explicit SubchannelPhysicsModel( const std::shared_ptr<ElementPhysicsModelParameters> &params );

    /**
      Destructor
      */
    virtual ~SubchannelPhysicsModel() {}

    /**
      Function to evaluate property functions
      @param [in]  property property identifier string
      @param [out] result   output vector
      @param [in]  args of input argument strings and input vectors
      */
    void getProperty( std::string property,
                      std::vector<double> &result,
                      std::map<std::string, std::shared_ptr<std::vector<double>>> &args );

    /**
      Function to return a pointer to the material
      */
    std::shared_ptr<AMP::Materials::Material> getMaterial() { return d_material; }

protected:
    // pointer to the material
    std::shared_ptr<AMP::Materials::Material> d_material;

    // map of property identifier strings and property pointers
    std::map<std::string, std::shared_ptr<AMP::Materials::Property<double>>> d_properties;
};
} // namespace Operator
} // namespace AMP

#endif
