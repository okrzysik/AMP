
#ifndef included_AMP_MassLinearFEOperator
#define included_AMP_MassLinearFEOperator

/* AMP files */
#include "operators/libmesh/LinearFEOperator.h"
#include "operators/libmesh/MassLinearFEOperatorParameters.h"
#include "operators/libmesh/MassLinearElement.h"
#include "utils/Utilities.h"

/* Boost files */
#include "boost/shared_ptr.hpp"

#include <vector>

namespace AMP {
  namespace Operator {

    class MassLinearFEOperator : public LinearFEOperator 
    {
      public :

        MassLinearFEOperator(const boost::shared_ptr<MassLinearFEOperatorParameters>& params);

        virtual ~MassLinearFEOperator() { }

        void preAssembly(const boost::shared_ptr<AMP::Operator::OperatorParameters>&);

        void postAssembly();

        void preElementOperation(const AMP::Mesh::MeshElement &);

        void postElementOperation();

        AMP::LinearAlgebra::Variable::shared_ptr getInputVariable();

        AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() ;

        boost::shared_ptr<MassDensityModel> getDensityModel() { return d_densityModel; };

      protected :

        bool d_useConstantTemperature;

        bool d_useConstantConcentration;

        bool d_useConstantBurnup;

        double d_constantTemperatureValue;

        double d_constantConcentrationValue;

        double d_constantBurnupValue;

        AMP::LinearAlgebra::Vector::shared_ptr d_temperature;

        AMP::LinearAlgebra::Vector::shared_ptr d_concentration;

        AMP::LinearAlgebra::Vector::shared_ptr d_burnup;

        std::vector<std::vector<double> > d_elementMassMatrix;

        boost::shared_ptr<MassLinearElement> d_massLinElem;

        boost::shared_ptr<MassDensityModel> d_densityModel;

        boost::shared_ptr<AMP::LinearAlgebra::Variable> d_inpVariable;

        boost::shared_ptr<AMP::LinearAlgebra::Variable> d_outVariable;

      private :

    };

  }
}

#endif


