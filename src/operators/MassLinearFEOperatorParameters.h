#ifndef included_AMP_MassLinearFEOperatorParameters
#define included_AMP_MassLinearFEOperatorParameters

#include "LinearFEOperatorParameters.h"

#include "MassDensityModel.h"

namespace AMP {
namespace Operator {
    
    class MassLinearFEOperatorParameters : public LinearFEOperatorParameters {
        public :
        
        MassLinearFEOperatorParameters(const boost::shared_ptr<AMP::Database> &db)
        : LinearFEOperatorParameters(db) {  }
        
        virtual ~MassLinearFEOperatorParameters() {}

        boost::shared_ptr<MassDensityModel> d_densityModel;
        
        protected :
        
        private :
        
    };
    
}
}

#endif


