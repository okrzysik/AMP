
#ifndef included_AMP_DiffusionLinearFEOperator
#define included_AMP_DiffusionLinearFEOperator

/* AMP files */
#include "operators/LinearFEOperator.h"
#include "operators/diffusion/DiffusionLinearFEOperatorParameters.h"
#include "operators/diffusion/DiffusionLinearElement.h"
#include "utils/Utilities.h"
#include "ampmesh/MeshElement.h"

/* Boost files */
#include "boost/shared_ptr.hpp"

#include <vector>


namespace AMP {
namespace Operator {


class DiffusionLinearFEOperator: public LinearFEOperator {
public:

        DiffusionLinearFEOperator(const boost::shared_ptr<
            DiffusionLinearFEOperatorParameters>& params);

        ~DiffusionLinearFEOperator() {
        }

        void preAssembly(const boost::shared_ptr<OperatorParameters>& params);

        void postAssembly();

        void preElementOperation(const AMP::Mesh::MeshElement &);

        void postElementOperation();

        AMP::LinearAlgebra::Variable::shared_ptr getInputVariable();

        AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable();

protected:

        bool d_useConstantTemperature;

        bool d_useConstantConcentration;

        bool d_useConstantBurnup;

        AMP::LinearAlgebra::Vector::shared_ptr d_temperature;

        AMP::LinearAlgebra::Vector::shared_ptr d_concentration;

        AMP::LinearAlgebra::Vector::shared_ptr d_burnup;

        std::vector<std::vector<double> > d_elementStiffnessMatrix;

        boost::shared_ptr<DiffusionLinearElement> d_diffLinElem;

        boost::shared_ptr<DiffusionTransportModel> d_transportModel;

        boost::shared_ptr<AMP::LinearAlgebra::Variable> d_inpVariable;

        boost::shared_ptr<AMP::LinearAlgebra::Variable> d_outVariable;

};


}
}

#endif

