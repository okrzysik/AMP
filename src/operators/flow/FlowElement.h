
#ifndef included_AMP_FlowElement
#define included_AMP_FlowElement

#include <vector>

#include "boost/shared_ptr.hpp"

/* AMP files */
#include "operators/ElementOperation.h"
#include "operators/flow/FlowTransportModel.h"

/* Libmesh files */
#include "elem.h"
#include "fe_type.h"
#include "fe_base.h"
#include "quadrature.h"

namespace AMP {
namespace Operator {

  class FlowElement : public ElementOperation 
  {
    public :


      FlowElement(const boost::shared_ptr<ElementOperationParameters>& params);

      virtual ~FlowElement() {  }

      void initializeForCurrentElement( const ::Elem *elem, 
          const boost::shared_ptr<FlowTransportModel> & transportModel )
      {
        d_elem = elem;
        d_transportModel = transportModel;
      }

    protected :

      boost::shared_ptr < ::FEType > d_feType; /**< Type of polynomial used for the
                                                               finite element shape functions. This includes
                                                               both the polynomial order: 
                                                               First order/Second order etc. and polynomial family:
                                                               Lagrange/Hierarchic/Hermite etc.  */

      boost::shared_ptr < ::FEBase > d_fe; /**< Finite element shape functions. */

      boost::shared_ptr < ::QBase > d_qrule; /**< Quadtrature rule used for numerical integration. */

      const ::Elem *d_elem; /**< Pointer to the current element within the finite element assembly. */

      boost::shared_ptr<FlowTransportModel> d_transportModel; 

    private :

  };

}
}

#endif


