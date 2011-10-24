#ifndef included_AMP_MassLinearFEOperatorParameters
#define included_AMP_MassLinearFEOperatorParameters

#include "FEOperatorParameters.h"

#include "MassDensityModel.h"

namespace AMP {
namespace Operator {
    
    class MassLinearFEOperatorParameters : public FEOperatorParameters {
        public :
        
        MassLinearFEOperatorParameters(const boost::shared_ptr<AMP::Database> &db)
        : FEOperatorParameters(db) {  }
        
        virtual ~MassLinearFEOperatorParameters() {}

        boost::shared_ptr<MassDensityModel> d_densityModel;
        
        protected :
        
        private :
        
    };
    
}
}

#endif


