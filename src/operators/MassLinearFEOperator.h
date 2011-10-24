
#ifndef included_AMP_MassLinearFEOperator
#define included_AMP_MassLinearFEOperator

/* AMP files */
#include "LinearFEOperator.h"
#include "MassLinearFEOperatorParameters.h"
#include "MassLinearElement.h"
#include "ampmesh/MeshVariable.h"
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

      void preElementOperation(const AMP::Mesh::MeshManager::Adapter::Element &, const std::vector<AMP::Mesh::DOFMap::shared_ptr> &);

      void postElementOperation();

      AMP::LinearAlgebra::Variable::shared_ptr createInputVariable (const std::string & name, int varId = -1) {
        (void) varId;      
        return d_inpVariable->cloneVariable(name);
      }

      AMP::LinearAlgebra::Variable::shared_ptr createOutputVariable (const std::string & name, int varId = -1) {
        (void) varId;      
        return d_outVariable->cloneVariable(name);
      }

      AMP::LinearAlgebra::Variable::shared_ptr getInputVariable(int varId = -1) {
        return d_inpVariable;
      }

      AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() {
        return d_outVariable;
      }

      void setInputVariableName(const std::string & name, int varId = -1) {
        (void) varId;      
        d_inpVariable->setName(name);
      }

      void setOutputVariableName(const std::string & name, int varId = -1) {
        (void) varId;      
        d_outVariable->setName(name);
      }

      unsigned int numberOfDOFMaps() {
        return 1;
      }

      AMP::LinearAlgebra::Variable::shared_ptr getVariableForDOFMap(unsigned int ) {
        return d_inpVariable;
      }

      boost::shared_ptr<MassDensityModel> getDensityModel(){return d_densityModel;}

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

      boost::shared_ptr<AMP::Mesh::NodalScalarVariable> d_inpVariable;

      boost::shared_ptr<AMP::Mesh::NodalScalarVariable> d_outVariable;

    private :

  };

}
}

#endif


