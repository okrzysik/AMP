#ifndef included_AMP_MassElement
#define included_AMP_MassElement

#include <vector>

#include "boost/shared_ptr.hpp"

/* AMP files */
#include "ElementOperation.h"
#include "MassDensityModel.h"

/* Libmesh files */
#include "fe_type.h"
#include "fe_base.h"
#include "quadrature_gauss.h"
#include "elem.h"

namespace AMP {
namespace Operator {
    
    class MassElement : public ElementOperation 
    {
        public :
        
        MassElement(const boost::shared_ptr<ElementOperationParameters>& params);
        
        virtual ~MassElement() {  }
        
        void initializeForCurrentElement( const ::Elem *elem,
                  const boost::shared_ptr<MassDensityModel> & densityModel );
        
        protected :
        
        boost::shared_ptr < ::FEType > d_feType;
        
        boost::shared_ptr < ::FEBase > d_fe;
        
        boost::shared_ptr < ::QBase > d_qrule;
        
        const std::vector<Real> *d_JxW;
        
        const std::vector<std::vector<Real> > *d_phi;
        
        const ::Elem *d_elem;
        
        boost::shared_ptr<MassDensityModel> d_densityModel;
        
        private :
        
    };
    
}
}

#endif




