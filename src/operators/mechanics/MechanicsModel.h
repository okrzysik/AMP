
#ifndef included_AMP_MechanicsModel
#define included_AMP_MechanicsModel

#include <cstring>

#include "operators/ElementPhysicsModel.h"
#include "operators/MechanicsModelParameters.h"

#include "materials/Material.h"
#include "vectors/Vector.h"

#include "boost/shared_ptr.hpp"

namespace AMP {
  namespace Operator {

    class MechanicsModel : public ElementPhysicsModel
    {
      public :

        /** Constructor */
        MechanicsModel(const boost::shared_ptr<MechanicsModelParameters>& params)
          : ElementPhysicsModel(params) { 
            bool useMaterialsLibrary = (params->d_db)->getBoolWithDefault("USE_MATERIALS_LIBRARY", false);
            if(useMaterialsLibrary == true) {
              AMP_INSIST( (params->d_db->keyExists("Material")), "Key ''Material'' is missing!" );
              std::string matname = params->d_db->getString("Material");
              d_material = AMP::voodoo::Factory<AMP::Materials::Material>::instance().create(matname);
            }
          }

        /** Destructor */
        virtual ~MechanicsModel() { }

        /** This function is used by the Mechanics operator to pass relevant
         * values to this material model. These values will be used to compute
         * the stress and/or tangent. */
        virtual void reset(boost::shared_ptr<MechanicsModelParameters> params) {
          d_deformationGradient = params->d_deformationGradient;
        }          

        /** This function should return the consistent/continuum tangent values in the vector
         * provided.   */
        virtual void getTangent(boost::shared_ptr<AMP::LinearAlgebra::Vector> tangent) { }

        /** This function should return the stress values in the vector
         * provided.   */
        virtual void getStress(boost::shared_ptr<AMP::LinearAlgebra::Vector> stress) { }

      protected :

        /** The material object that can be used to get material properties */
        boost::shared_ptr<AMP::Materials::Material> d_material; 

        /** Pointer to the deformation gradient that is required to compute the
         * stress and/or tangent. */
        boost::shared_ptr<AMP::LinearAlgebra::Vector> d_deformationGradient;
    };

  }
}

#endif


