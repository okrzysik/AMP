
#ifndef included_AMP_MassLinearFEOperator
#define included_AMP_MassLinearFEOperator

/* AMP files */
#include "LinearFEOperator.h"
#include "MassLinearFEOperatorParameters.h"
#include "MassLinearElement.h"
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

      ~MassLinearFEOperator() { }

      void preAssembly(const boost::shared_ptr<AMP::Operator::OperatorParameters>&);

      void postAssembly();

      void preElementOperation(const AMP::Mesh::MeshElement &, const std::vector<AMP::Discretization::DOFManager::shared_ptr> &);

      void postElementOperation();

      AMP::LinearAlgebra::Variable::shared_ptr createInputVariable (const std::string & name, int varId = -1);

      AMP::LinearAlgebra::Variable::shared_ptr createOutputVariable (const std::string & name, int varId = -1);

      AMP::LinearAlgebra::Variable::shared_ptr getInputVariable(int varId = -1);

      AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() ;

      void setInputVariableName(const std::string & name, int varId = -1);

      void setOutputVariableName(const std::string & name, int varId = -1);

      unsigned int numberOfDOFMaps();

      AMP::LinearAlgebra::Variable::shared_ptr getVariableForDOFMap(unsigned int );

      boost::shared_ptr<MassDensityModel> getDensityModel();

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

      std::vector<unsigned int> d_dofIndices;

      std::vector<std::vector<double> > d_elementMassMatrix;

      boost::shared_ptr<MassLinearElement> d_massLinElem;

      boost::shared_ptr<MassDensityModel> d_densityModel;

      //boost::shared_ptr<AMP::Mesh::NodalScalarVariable> d_inpVariable;

      //boost::shared_ptr<AMP::Mesh::NodalScalarVariable> d_outVariable;

    private :

  };

}
}

#endif


