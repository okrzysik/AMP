
#ifndef included_AMP_MechanicsLinearFEOperator
#define included_AMP_MechanicsLinearFEOperator

/* AMP files */
#include "operators/LinearFEOperator.h"
#include "MechanicsConstants.h"
#include "MechanicsLinearFEOperatorParameters.h"
#include "MechanicsLinearElement.h"
#include "MechanicsLinearUpdatedLagrangianElement.h"
#include "ampmesh/MeshVariable.h"
#include "vectors/MultiVariable.h"

#include <vector>

namespace AMP {
namespace Operator {

  /**
    A class for representing the linear operator for linear/nonlinear mechanics.
    In the case of nonlinear mechanics, this operator will result from the
    linearization (or some approximate linearization) of the nonlinear operator.
    This class can be used to compute the finite element (FE) stiffness 
    matrix corresponding to the mechanical equilibrium equations for a 
    solid body. This class only deals with the volume integration, 
    the boundary conditions are handled separately by the boundary operators. 
    */
  class MechanicsLinearFEOperator : public LinearFEOperator 
  {
    public :

      /**
        Constructor. This allocates memory for the stiffness matrix. This also computes the entries of the stiffness matrix unless
        (a) this operator is the jacobian of the nonlinear mechanics operator and (b) the nonlinear mechanics operator is not already
        initialized at the time of construction of this operator. This reads the values for the following keys from the database object contained in
        the parameter object, params:
        1) isAttachedToNonlinearOperator (false by default) - Is this a jacobian of the nonlinear mechanics operator?
        2) isNonlinearOperatorInitialized (false by default) - If this is a jacobian of the nonlinear mechanics 
        operator, is the nonlinear mechanics operator already initialized at the time of construction of this operator?
        3) InputVariable (No default value) - Name of the input variable
        4) OutputVariable (No default value) - Name of the output variable
        */
      MechanicsLinearFEOperator(const boost::shared_ptr<MechanicsLinearFEOperatorParameters>& params);

      /**
        Destructor
        */
      ~MechanicsLinearFEOperator() { }

      /**
        This is called at the start of the FE assembly. The matrix is set to 0.
        */
      void preAssembly(const boost::shared_ptr<OperatorParameters>& params);

      /**
        This is called at the end of the FE assembly. The entries of the matrix corresponding
        to nodes that are shared between two or more processors are made consistent.
        */
      void postAssembly();

      /**
        This function will be called once for each element, just before performing
        the element operation. This function extracts the local information from
        the global mesh objects (DOFMap, global vectors and matrices) and
        passes them to the element operation.
        */
      void preElementOperation(const AMP::Mesh::MeshManager::Adapter::Element &, const std::vector<AMP::Mesh::DOFMap::shared_ptr> &);

      /**
        This function will be called once for each element, just after performing the element operation.
        The element stiffness matrix is added to the global stiffness matrix in this function.
        */
      void postElementOperation();

      /**
        Sets the name for the input variable. 
        @param [in] name Name of the input variable.
        @param [in] varId This parameter is not used.
        */
      void setInputVariableName(const std::string & name, int varId = -1) {
        (void) varId;	  
        d_inpVariable->setName(name);
      }

      /**
        Sets the name for the output variable.
        @param [in] name Name of the output variable
        @param [in] varId This parameter is not used.
        */
      void setOutputVariableName(const std::string & name, int varId = -1) {
        (void) varId;	  
        d_outVariable->setName(name);
      }

      /**
        Creates a variable that has the same type as the input variable and sets its name.
        @param [in] name Name of the variable.
        @param [in] varId This parameter is not used.
        */
      static AMP::LinearAlgebra::Variable::shared_ptr createInputVariable (const std::string & name, int varId = -1) {
        (void) varId;	  
        AMP::LinearAlgebra::Variable::shared_ptr inpVar(new AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable, 3>(name) );
        return inpVar;
      }

      /**
        Creates a variable that has the same type as the output variable and sets its name.
        @param [in] name Name of the variable. 
        @param [in] varId This parameter is not used.
        */
      static AMP::LinearAlgebra::Variable::shared_ptr createOutputVariable (const std::string & name, int varId = -1) {
        (void) varId;	  
        AMP::LinearAlgebra::Variable::shared_ptr outVar(new AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable, 3>(name) );
        return outVar;
      }

      /**
        @return The variable for the input vector. 
        */
      AMP::LinearAlgebra::Variable::shared_ptr getInputVariable(int varId = -1) {
        return d_inpVariable;
      }

      /**
        @return The variable for the output vector
        */
      AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() {
        return d_outVariable;
      }

      /**
        @return The number of different DOFMaps required to assemble this operator = 1. 
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
      AMP::LinearAlgebra::Variable::shared_ptr getVariableForDOFMap(unsigned int ) {
        return d_inpVariable;
      }

      /**
        Writes the stress and strain at each Gauss point to a file.
        The 6 components of stress and strain at each Gauss point are arranged in the order:
        xx, yy, zz, yz, xz and  xy.
        @param [in] disp Displacement vector
        @param [in] fname Name of the output file
        */
      void printStressAndStrain(AMP::LinearAlgebra::Vector::shared_ptr disp, const std::string & fname);

      /**
        Creates vectors containing the stress and strain values at each Gauss point.
        The 6 components of stress and strain at each Gauss point are arranged in the order:
        xx, yy, zz, yz, xz and  xy.
        @param [in] disp Displacements 
        @param [out] stress Stresses
        @param [out] strain Strains
        */
      void computeStressesAndStrains(AMP::LinearAlgebra::Vector::shared_ptr disp,
          AMP::LinearAlgebra::Vector::shared_ptr & stress, AMP::LinearAlgebra::Vector::shared_ptr & strain);

    protected :

      std::vector<unsigned int> d_dofIndices[3]; /**< DOF indices */

      std::vector<std::vector<double> > d_elementStiffnessMatrix; /**< Element stiffness matrix. */

      boost::shared_ptr< MechanicsLinearElement > d_mechLinElem; /**< Element operation. */

      boost::shared_ptr< MechanicsLinearUpdatedLagrangianElement > d_mechLinULElem; /**< Linear Updated Lagrangian Element operation. */

      boost::shared_ptr< MechanicsMaterialModel > d_materialModel; /**< Material model. */

      AMP::LinearAlgebra::Vector::shared_ptr d_refXYZ; /**< Reference x, y and z coordinates. */

      AMP::LinearAlgebra::Vector::shared_ptr d_dispVec;

      bool d_useUpdatedLagrangian;

      unsigned int d_numNodesForCurrentElement; /**< Number of nodes in the current element. */

    private :

      boost::shared_ptr<AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable, 3> > d_inpVariable; /**< Input variable */

      boost::shared_ptr<AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable, 3> > d_outVariable; /**< Output variable */

  };

}
}

#endif

