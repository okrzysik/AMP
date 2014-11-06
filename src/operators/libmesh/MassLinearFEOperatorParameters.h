#ifndef included_AMP_MassLinearFEOperatorParameters
#define included_AMP_MassLinearFEOperatorParameters

#include "operators/libmesh/LinearFEOperatorParameters.h"

#include "operators/libmesh/MassDensityModel.h"

namespace AMP {
namespace Operator {
    
    class MassLinearFEOperatorParameters : public LinearFEOperatorParameters {
        public :
        
        MassLinearFEOperatorParameters(const AMP::shared_ptr<AMP::Database> &db)
        : LinearFEOperatorParameters(db) {  }
        
        virtual ~MassLinearFEOperatorParameters() {}

        AMP::shared_ptr<MassDensityModel> d_densityModel;
        
        protected :
        
        private :
        
    };
    
}
}

#endif


