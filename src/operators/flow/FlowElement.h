
#ifndef included_AMP_FlowElement
#define included_AMP_FlowElement

#include <vector>

#include "utils/shared_ptr.h"

/* AMP files */
#include "operators/ElementOperation.h"
#include "operators/flow/FlowTransportModel.h"

/* Libmesh files */
#include "libmesh/elem.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_type.h"
#include "libmesh/quadrature.h"

namespace AMP {
namespace Operator {

class FlowElement : public ElementOperation
{
public:
    explicit FlowElement( const AMP::shared_ptr<ElementOperationParameters> &params );

    virtual ~FlowElement() {}

    void initializeForCurrentElement( const ::Elem *elem,
                                      const AMP::shared_ptr<FlowTransportModel> &transportModel )
    {
        d_elem           = elem;
        d_transportModel = transportModel;
    }

protected:
    AMP::shared_ptr<::FEType>
        d_feType; /**< Type of polynomial used for the
                                    finite element shape functions. This includes
                                    both the polynomial order:
                                    First order/Second order etc. and polynomial family:
                                    Lagrange/Hierarchic/Hermite etc.  */

    AMP::shared_ptr<::FEBase> d_fe; /**< Finite element shape functions. */

    AMP::shared_ptr<::QBase> d_qrule; /**< Quadtrature rule used for numerical integration. */

    const ::Elem *d_elem; /**< Pointer to the current element within the finite element assembly. */

    AMP::shared_ptr<FlowTransportModel> d_transportModel;

private:
};
}
}

#endif
