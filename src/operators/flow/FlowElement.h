
#ifndef included_AMP_FlowElement
#define included_AMP_FlowElement

#include <vector>

#include "AMP/operators/ElementOperation.h"
#include "AMP/operators/flow/FlowTransportModel.h"
#include <memory>

// Libmesh headers
DISABLE_WARNINGS
#include "libmesh/auto_ptr.h"
#include "libmesh/cell_hex8.h"
#include "libmesh/elem.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_type.h"
#include "libmesh/quadrature.h"
#include "libmesh/string_to_enum.h"
ENABLE_WARNINGS

namespace AMP {
namespace Operator {

class FlowElement : public ElementOperation
{
public:
    explicit FlowElement( const std::shared_ptr<ElementOperationParameters> &params );

    virtual ~FlowElement() {}

    void initializeForCurrentElement( const libMesh::Elem *elem,
                                      const std::shared_ptr<FlowTransportModel> &transportModel )
    {
        d_elem           = elem;
        d_transportModel = transportModel;
    }

protected:
    std::shared_ptr<libMesh::FEType>
        d_feType; /**< Type of polynomial used for the
                                    finite element shape functions. This includes
                                    both the polynomial order:
                                    First order/Second order etc. and polynomial family:
                                    Lagrange/Hierarchic/Hermite etc.  */

    std::shared_ptr<libMesh::FEBase> d_fe; /**< Finite element shape functions. */

    std::shared_ptr<libMesh::QBase>
        d_qrule; /**< Quadtrature rule used for numerical integration. */

    const libMesh::Elem
        *d_elem; /**< Pointer to the current element within the finite element assembly. */

    std::shared_ptr<FlowTransportModel> d_transportModel;

private:
};
} // namespace Operator
} // namespace AMP

#endif
