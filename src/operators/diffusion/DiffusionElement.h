#ifndef included_AMP_DiffusionElement
#define included_AMP_DiffusionElement

#include <vector>

#include "boost/shared_ptr.hpp"

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

    DiffusionElement(
            const boost::shared_ptr<ElementOperationParameters>& params);

    virtual ~DiffusionElement() {
    }

    void initializeForCurrentElement(const ::Elem *elem,
            const boost::shared_ptr<DiffusionTransportModel> & transportModel);

protected:

    boost::shared_ptr< ::FEType> d_feType;

    boost::shared_ptr< ::FEBase> d_fe;

    boost::shared_ptr< ::QBase> d_qrule;

    const std::vector<Real> *d_JxW;

    const std::vector<std::vector<Real> > *d_phi;

    const std::vector<std::vector<RealGradient> > *d_dphi;

    const ::Elem *d_elem;

    boost::shared_ptr<DiffusionTransportModel> d_transportModel;
    boost::shared_ptr<DiffusionTransportTensorModel> d_transportTensorModel;

private:

};

}
}

#endif

