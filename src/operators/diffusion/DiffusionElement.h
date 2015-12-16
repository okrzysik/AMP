#ifndef included_AMP_DiffusionElement
#define included_AMP_DiffusionElement

#include <vector>

#include "utils/shared_ptr.h"

/* AMP files */
#include "operators/ElementOperation.h"
#include "operators/diffusion/DiffusionTransportModel.h"
#include "operators/diffusion/DiffusionTransportTensorModel.h"

/* Libmesh files */
#include "libmesh/fe_type.h"
#include "libmesh/fe_base.h"
#include "libmesh/elem.h"
#include "libmesh/quadrature_gauss.h"

namespace AMP {
namespace Operator {

class DiffusionElement: public ElementOperation {
public:

    explicit DiffusionElement(
            const AMP::shared_ptr<ElementOperationParameters>& params);

    virtual ~DiffusionElement() {
    }

    void initializeForCurrentElement(const ::Elem *elem,
            const AMP::shared_ptr<DiffusionTransportModel> & transportModel);

protected:

    AMP::shared_ptr< ::FEType> d_feType;

    AMP::shared_ptr< ::FEBase> d_fe;

    AMP::shared_ptr< ::QBase> d_qrule;

    const std::vector<Real> *d_JxW;

    const std::vector<std::vector<Real> > *d_phi;

    const std::vector<std::vector<RealGradient> > *d_dphi;

    const ::Elem *d_elem;

    AMP::shared_ptr<DiffusionTransportModel> d_transportModel;
    AMP::shared_ptr<DiffusionTransportTensorModel> d_transportTensorModel;

private:

};

}
}

#endif

