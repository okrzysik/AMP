#ifndef included_AMP_MassElement
#define included_AMP_MassElement

#include <vector>

#include "utils/shared_ptr.h"

/* AMP files */
#include "operators/ElementOperation.h"
#include "operators/libmesh/MassDensityModel.h"

/* Libmesh files */
#include "libmesh/fe_type.h"
#include "libmesh/fe_base.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/elem.h"

namespace AMP {
namespace Operator {
    
    class MassElement : public ElementOperation 
    {
        public :
        
        explicit MassElement(const AMP::shared_ptr<ElementOperationParameters>& params);
        
        virtual ~MassElement() {  }
        
        void initializeForCurrentElement( const ::Elem *elem,
                  const AMP::shared_ptr<MassDensityModel> & densityModel );
        
        protected :
        
        AMP::shared_ptr < ::FEType > d_feType;
        
        AMP::shared_ptr < ::FEBase > d_fe;
        
        AMP::shared_ptr < ::QBase > d_qrule;
        
        const std::vector<Real> *d_JxW;
        
        const std::vector<std::vector<Real> > *d_phi;
        
        const ::Elem *d_elem;
        
        AMP::shared_ptr<MassDensityModel> d_densityModel;
        
        private :
        
    };
    
}
}

#endif




