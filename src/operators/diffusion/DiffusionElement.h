#ifndef included_AMP_DiffusionElement
#define included_AMP_DiffusionElement

#include <vector>

#include "AMP/operators/ElementOperation.h"
#include "AMP/operators/diffusion/DiffusionTransportModel.h"
#include "AMP/operators/diffusion/DiffusionTransportTensorModel.h"
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

class DiffusionElement : public ElementOperation
{
public:
    explicit DiffusionElement( const std::shared_ptr<ElementOperationParameters> &params );

    virtual ~DiffusionElement() {}

    void
    initializeForCurrentElement( const ::Elem *elem,
                                 const std::shared_ptr<DiffusionTransportModel> &transportModel );

protected:
    std::shared_ptr<::FEType> d_feType;

    std::shared_ptr<::FEBase> d_fe;

    std::shared_ptr<::QBase> d_qrule;

    const std::vector<Real> *d_JxW;

    const std::vector<std::vector<Real>> *d_phi;

    const std::vector<std::vector<RealGradient>> *d_dphi;

    const ::Elem *d_elem;

    std::shared_ptr<DiffusionTransportModel> d_transportModel;
    std::shared_ptr<DiffusionTransportTensorModel> d_transportTensorModel;

private:
};
} // namespace Operator
} // namespace AMP

#endif
