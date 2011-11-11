
#ifndef included_AMP_VolumeIntegralOperator
#define included_AMP_VolumeIntegralOperator

#include "NonlinearFEOperator.h"
#include "VolumeIntegralOperatorParameters.h"
#include "SourceNonlinearElement.h"
#include "vectors/MultiVariable.h"
#include "matrices/Matrix.h"

#include <vector>

namespace AMP {
namespace Operator {

  /**
    This class is used to integrate any arbitrary functional over the entire domain. The 
    evaluation of the functional at quadrature points and volume integration are its key 
    basic operations.
    */
  class VolumeIntegralOperator : public NonlinearFEOperator
  {
    public :

      /**
        Constructor. This reads the values for the following keys from the database object contained in
        the parameter object, params:
        1) Primary(Active)InputVariables - List of active input variables names. The supported variable types are:
        NodalScalar, Nodalvector, IntegrationPointScalar, IntegrationPointVector.
        2) AuxillaryInputVariables - List of auxillary input variables names. These are frozen variables and are 
        temporarily not used in any formulation.
        3) OutputVariable - Name of the output variable
        */
      VolumeIntegralOperator(const boost::shared_ptr<VolumeIntegralOperatorParameters> & params);

      /**
        Destructor.
        */
      ~VolumeIntegralOperator() { }

      /**
        This function is called at the beginning of the FE assembly
        */
      void preAssembly(const boost::shared_ptr<AMP::LinearAlgebra::Vector> &u, boost::shared_ptr<AMP::LinearAlgebra::Vector> &r);

      /**
        This function is called at the end of the FE assembly
        */
      void postAssembly();

      /**
        This function is called at the beginning of the element computation
        */
      void preElementOperation(const AMP::Mesh::MeshElement &, const std::vector<AMP::Discretization::DOFManager::shared_ptr> &);

      /**
        This function is called at the end of the element computation
        */
      void postElementOperation();

      /**
        This is used to update the operator between successive solves with the operator. 
        */
      void reset(const boost::shared_ptr<OperatorParameters>&);

      /**
        This is used to compute the information required to reset the corresponding Linear (Jacobian) operator
        */
      boost::shared_ptr<OperatorParameters>
          getJacobianParameters(const boost::shared_ptr<AMP::LinearAlgebra::Vector>&);

      void setInputVariableName(const std::string & name, int varId = -1) {
          if(varId == -1) {
              d_inpVariables->setName(name);
          } else {
              (d_inpVariables->getVariable(varId))->setName(name);
          }
      }

      void setOutputVariableName(const std::string & name, int varId = -1) {
          (void) varId;
          d_outVariable->setName(name);
      }

      void setAuxVariable(const boost::shared_ptr<AMP::LinearAlgebra::MultiVariable>& var) {
        d_auxVariables = var;
      }

      AMP::LinearAlgebra::Variable::shared_ptr createInputVariable(const std::string & name, int varId = -1);

      AMP::LinearAlgebra::Variable::shared_ptr createOutputVariable(const std::string & name, int varId = -1) {
    (void) varId;
          return d_outVariable->cloneVariable(name);
      }

      AMP::LinearAlgebra::Variable::shared_ptr getInputVariable(int varId) {
        if(varId == -1) {
          return d_inpVariables; 
        } else {
          return d_inpVariables->getVariable(varId);
        }
      }

      AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() {
          return d_outVariable;
      }

      /**
        @return The number of different DOFMaps required to assemble this operator. Returns 2 if
        at least one of Temperature, Burnup and Oxygen concentration is active and 1 otherwise.
        */
      unsigned int numberOfDOFMaps() {
          return 1;
      }

      /**
        @param [in] id Identifier for the type of DOFMap required. It is
        a number between 0 (included) and numberOfDOFMaps (excluded)
        @return The variable corresponding to the DOFMap specified by id. 
        @see numberOfDOFMaps
        */
      AMP::LinearAlgebra::Variable::shared_ptr getVariableForDOFMap(unsigned int id){
          return (d_inpVariables->getVariable(id));
      }

      boost::shared_ptr<SourcePhysicsModel> getSourcePhysicsModel(){
          return d_sourcePhysicsModel;
      }

      boost::shared_ptr<SourceNonlinearElement> getSourceElement() { 
          return d_srcNonlinElem;
      }

      boost::shared_ptr<AMP::LinearAlgebra::Matrix> getLinearizedVolumeIntegralOperator(const boost::shared_ptr<OperatorParameters>& params);

    protected :

      double d_sourceValues;

      bool d_csource;

      std::string d_isInputType; 

      void init(const boost::shared_ptr<VolumeIntegralOperatorParameters>& params);

      std::vector<unsigned int> d_inpDofIndices; /**< DOF indices for the input variable. */

      std::vector<unsigned int> d_outDofIndices; /**< DOF indices for the output variable. */

      unsigned int d_numNodesForCurrentElement;/**< Number of nodes in the current element. */

      unsigned int d_numDofsForCurrentElement;/**< Number of Dofs in the current element. */

      std::vector<std::string> d_activeVariableNames;/**< A list of strings to store Active Variable names. */

      unsigned int d_numPrimaryVariables;/**< Number of Active Variables. */

      unsigned int d_numAuxillaryVariables;/**< Number of Auxillary(Frozen) Variables. */

      std::vector<double> d_elementOutputVector;/**< Output vector for the Element Operation. */

      boost::shared_ptr<SourceNonlinearElement> d_srcNonlinElem;/**< Element Operation. */

      boost::shared_ptr<SourcePhysicsModel> d_sourcePhysicsModel;/**< Source Physics Model. */

      std::vector<AMP::LinearAlgebra::Vector::shared_ptr> d_inVec;/**< Input vector for active variables. */

      std::vector<AMP::LinearAlgebra::Vector::shared_ptr> d_auxVec;/**< Input vector for auxillary variables. */

      AMP::LinearAlgebra::Vector::shared_ptr d_multiAuxPtr;

      AMP::LinearAlgebra::Vector::shared_ptr d_outVec; /**< Output vector. */

    private :

      boost::shared_ptr<AMP::LinearAlgebra::MultiVariable> d_inpVariables;/**< Input Active variables. */

      boost::shared_ptr<AMP::LinearAlgebra::MultiVariable> d_auxVariables;/**< Input Auxillary variables. */
     
      boost::shared_ptr<AMP::LinearAlgebra::Variable> d_outVariable;/**< Output variables. */

      boost::shared_ptr<AMP::LinearAlgebra::Matrix> d_pDiagonalMatrix;
      boost::shared_ptr<AMP::LinearAlgebra::Vector> d_pDiagonalVector;
      boost::shared_ptr<AMP::LinearAlgebra::Vector> d_pNullVector;
 
      bool d_bMatrixAndVectorsCloned;
  };

}
}

#endif

