
#ifndef included_AMP_MechanicsNonlinearFEOperator
#define included_AMP_MechanicsNonlinearFEOperator

/* AMP files */
#include "operators/NonlinearFEOperator.h"
#include "MechanicsConstants.h"
#include "MechanicsNonlinearFEOperatorParameters.h"
#include "MechanicsNonlinearElement.h"
#include "MechanicsNonlinearUpdatedLagrangianElement.h"
#include "ampmesh/MeshVariable.h"
#include "vectors/MultiVariable.h"

#include "elem.h"
#include <vector>

namespace AMP {
  namespace Operator {

    /**
      A class used for representing the nonlinear mechanics operator.
      This class can be used to compute the finite element (FE) residual
      vector corresponding to the mechanical equilibrium equations for a 
      solid body. This class only deals with the volume integration, 
      the boundary conditions are handled separately by the boundary operators. 
      */
    class MechanicsNonlinearFEOperator : public NonlinearFEOperator 
    {
      public :

        /**
          Constructor. This reads the values for the following keys from the database object contained in
          the parameter object, params:
          1) RESET_REUSES_RADIAL_RETURN (true by default) - Can we assume that the apply() function is called before a call to 
          the reset() function and can we reuse the results of the radial return algorithm computed in the apply()
          function in the reset() function? Note, this is typically true unless the reset() function
          is evaluated at a different state from the last call to apply().
          2) JACOBIAN_REUSES_RADIAL_RETURN (true by default) - Can we assume that the apply() function is called before a call to 
          the getJacobianParameters() function and can we reuse the results of the radial return algorithm computed in the apply()
          function in the getJacobianParameters() function? Note, this is typically true unless the getJacobianParameters() function
          is evaluated at a different state from the last call to apply().
          3) ActiveInputVariables (No default value) - List of active input variables names. The supported variable types are:
          DISPLACEMENT, TEMPERATURE, BURNUP, OXYGEN_CONCENTRATION and LHGR. DISPLACEMENT must be active.
          4) FREEZE_TEMPERATURE/FREEZE_BURNUP/FREEZE_OXYGEN_CONCENTRATION/FREEZE_LHGR (true by default) - Are these variables frozen? This 
          will be ignored if the corresponding variable is not active.
          5) OutputVariable (No default value) - Name of the output variable
          */
        MechanicsNonlinearFEOperator(const boost::shared_ptr<MechanicsNonlinearFEOperatorParameters>& params);

        /**
          Destructor.
          */
        ~MechanicsNonlinearFEOperator() { }

        /**
          This function is called at the beginning of the FE assembly. The output vector, r, is set to 0.
          The values of the input vector, u, on the nodes shared between two or more processors are made consistent.
          @param [in] u  input vector
          @param [out] r output vector
          */
        void preAssembly(const boost::shared_ptr<AMP::LinearAlgebra::Vector>  &u, boost::shared_ptr<AMP::LinearAlgebra::Vector>  &r);

        /**
          This function is called at the end of the FE assembly.
          The values of the output vector on the nodes shared between two or more processors are made consistent.
          */
        void postAssembly();

        /**
          This function is called at the beginning of the element computation. The part of the
          input vector that is relevant for the computation in the current element is extracted 
          and passed to MechanicsNonlinearElement.
          */
        void preElementOperation(const AMP::Mesh::MeshManager::Adapter::Element &, const std::vector<AMP::Mesh::DOFMap::shared_ptr> &);

        /**
          This function is called at the end of the element computation. The entries of the 
          element output vector are added to the corresponding entries of the global output vector.
          */
        void postElementOperation();

        /**
          This is used to update the operator between successive solves with the operator. 
          */
        void reset(const boost::shared_ptr<OperatorParameters>& );

        /**
          This is used to compute the information required to reset the corresponding Linear (Jacobian) operator
          */
        boost::shared_ptr<OperatorParameters> 
          getJacobianParameters(const boost::shared_ptr<AMP::LinearAlgebra::Vector>& );

        /**
          This performs a dummy loop over the elements and gauss points so that the mechanics material model classes can 
          allocate memory and/or initialize their data as required.
          */
        void init();

        /**
          This function is used to set the reference temperature when using temperature dependent material models.
          @param [in] refTemp Reference temperature
          */
        void setReferenceTemperature(AMP::LinearAlgebra::Vector::shared_ptr refTemp) {
          d_referenceTemperature = refTemp->subsetVectorForVariable(d_inpVariables->getVariable(Mechanics::TEMPERATURE));
          d_referenceTemperature->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
          if(d_useUpdatedLagrangian) {
            d_inVec_pre[Mechanics::TEMPERATURE]->copyVector(refTemp);
            d_inVec_pre[Mechanics::TEMPERATURE]->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
          }
        }

        /**
          This function is used to set frozen vectors in this operator. This is used when some of the 
          variables are solved for in an uncoupled manner.
          @param [in] id Variable Identifier - One of AMP::Mechanics::DISPLACEMENT/TEMPERATURE/BURNUP/OXYGEN_CONCENTRATION/LHGR
          @param [in] frozenVec Frozen vector
          @see MechanicsConstants.h
          */
        void setVector(unsigned int id, AMP::LinearAlgebra::Vector::shared_ptr &frozenVec) {
          d_inVec[id] = frozenVec->subsetVectorForVariable(d_inpVariables->getVariable(id));
          (d_inVec[id])->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
        }

        /**
          Sets the name for each component of the input vector. If varId is equal to -1, it sets
          the name for the multivariable corresponding to the entire vector.
          @param [in] name Name of the component
          @param [in] varId Identifier for the component - One of AMP::Mechanics::DISPLACEMENT/TEMPERATURE/BURNUP/OXYGEN_CONCENTRATION/LHGR
          @see MechanicsConstants.h
          */
        void setInputVariableName(const std::string & name, int varId = -1) {
          if(varId == -1) {
            d_inpVariables->setName(name);
          } else {
            (d_inpVariables->getVariable(varId))->setName(name);
          }
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
          Creates a variable that has the same type as the specified component of the input vector and sets its name.
          @param [in] name Name of the component
          @param [in] varId Identifier for the component - One of AMP::Mechanics::DISPLACEMENT/TEMPERATURE/BURNUP/OXYGEN_CONCENTRATION/LHGR
          @see MechanicsConstants.h
          */
        static AMP::LinearAlgebra::Variable::shared_ptr createInputVariable(const std::string & name, int varId = -1);

        /**
          Creates a variable that has the same type as the output variable and sets its name.
          @param [in] name Name of the variable. 
          @param [in] varId This parameter is not used.
          */
        static AMP::LinearAlgebra::Variable::shared_ptr createOutputVariable(const std::string & name, int varId = -1) {
          (void) varId;	  
          AMP::LinearAlgebra::Variable::shared_ptr outVar(new AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable, 3>(name) );
          return outVar;
        }

        /**
          @param [in] varId Identifier for the component - One of AMP::Mechanics::DISPLACEMENT/TEMPERATURE/BURNUP/OXYGEN_CONCENTRATION/LHGR
          @return The variable for the specified component of the input vector. If varId is equal to -1, it
          returns the multivariable for the entire vector.
          */
        AMP::LinearAlgebra::Variable::shared_ptr getInputVariable(int varId = -1) {
          if(varId == -1) {
            return d_inpVariables; 
          } else {
            return d_inpVariables->getVariable(varId);
          }
        }

        /**
          @return The variable for the output vector
          */
        AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() {
          return d_outVariable;
        }

        /**
          @return The number of different DOFMaps required to assemble this operator. Returns 2 if
          at least one of Temperature, Burnup, Oxygen concentration and LHGR is active and 1 otherwise.
          */
        unsigned int numberOfDOFMaps();

        /**
          @param [in] id Identifier for the type of DOFMap required. It is
          a number between 0 (included) and numberOfDOFMaps (excluded)
          @return The variable corresponding to the DOFMap specified by id. 
          @see numberOfDOFMaps
          */
        AMP::LinearAlgebra::Variable::shared_ptr getVariableForDOFMap(unsigned int id);

        /**
          Writes the stress and strain at each Gauss point to a file.
          The 6 components of stress and strain at each Gauss point are arranged in the order:
          xx, yy, zz, yz, xz and  xy.
          @param [in] u Input vector
          @param [in] fname Name of the output file
          */
        void printStressAndStrain(AMP::LinearAlgebra::Vector::shared_ptr u, const std::string & fname);

        /**
          Creates vectors containing the stress and strain values at each Gauss point.
          The 6 components of stress and strain at each Gauss point are arranged in the order:
          xx, yy, zz, yz, xz and  xy.
          @param [in] u Input vector 
          @param [out] stress Stresses
          @param [out] strain Strains
          */
        void computeStressesAndStrains(AMP::LinearAlgebra::Vector::shared_ptr u,
            AMP::LinearAlgebra::Vector::shared_ptr & stress, AMP::LinearAlgebra::Vector::shared_ptr & strain);

        boost::shared_ptr<MechanicsMaterialModel> getMaterialModel() { return d_materialModel; }

      protected :

        template <MechanicsNonlinearElement::MaterialUpdateType updateType>
          void updateMaterialForElement(const AMP::Mesh::MeshManager::Adapter::Element & , 
              const std::vector<AMP::Mesh::DOFMap::shared_ptr> & );

        template <MechanicsNonlinearUpdatedLagrangianElement::MaterialUpdateType updateType>
          void updateMaterialForUpdatedLagrangianElement(const AMP::Mesh::MeshManager::Adapter::Element & ,
              const std::vector<AMP::Mesh::DOFMap::shared_ptr> & );

        void updateMaterialForElementCommonFunction(const AMP::Mesh::MeshManager::Adapter::Element & , const std::vector<AMP::Mesh::DOFMap::shared_ptr> &, 
            std::vector<std::vector<double> >&, std::vector<std::vector<double> >& );

        std::vector<unsigned int> d_type0DofIndices[3]; /**< DOF indices for the DISPLACEMENT variable type. */
        std::vector<unsigned int> d_type1DofIndices; /**< DOF indices for the TEMPERATURE/BURNUP/OXYGEN_CONCENTRATION/LHGR variable types. */

        unsigned int d_numNodesForCurrentElement; /**< Number of nodes in the current element. */

        std::vector<double> d_elementOutputVector; /**< Element output vector. */

        boost::shared_ptr<MechanicsNonlinearElement> d_mechNonlinElem; /**< Element operation. */

        boost::shared_ptr<MechanicsNonlinearUpdatedLagrangianElement> d_mechNULElem; /**< Nonlinear Updated Lagrangian Element operation. */

        boost::shared_ptr<MechanicsMaterialModel> d_materialModel; /**< Material model. */

        std::vector<AMP::LinearAlgebra::Vector::shared_ptr> d_inVec; /**< Input vector. */

        std::vector<AMP::LinearAlgebra::Vector::shared_ptr> d_inVec_pre; /**< Input vector. at nth (previous) configuration */

        AMP::LinearAlgebra::Vector::shared_ptr d_refXYZ; /**< Reference x, y and z coordinates. */

        AMP::LinearAlgebra::Vector::shared_ptr d_referenceTemperature; /**< Reference temperature. */

        AMP::LinearAlgebra::Vector::shared_ptr d_outVec; /**< Output vector. */

        bool d_resetReusesRadialReturn; /**< A flag that is true if the reset() function can reuse 
                                          the results from the radial return computation
                                          in the last call to the apply() function and false otherwise. */

        bool d_jacobianReusesRadialReturn; /**< A flag that is true if the getJacobianParameters() function can reuse 
                                             the results from the radial return computation
                                             in the last call to the apply() function and false otherwise. */

        std::vector<bool> d_isActive; /**< A list of flags to determine which variables are active. */

        std::vector<bool> d_isFrozen; /**< A list of flags to determine which variables are frozen. */

        bool d_useUpdatedLagrangian; /**< A flag that checks whether to use Updated Lagrangian formulation or not. */

      private :

        bool d_isInitialized; /**< A flag that is true if init() has been called and false otherwsie. */

        boost::shared_ptr<AMP::LinearAlgebra::MultiVariable> d_inpVariables; /**< Input variables. */

        boost::shared_ptr<AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable, 3> > d_outVariable; /**< Output variable. */

    };

    template <MechanicsNonlinearElement::MaterialUpdateType updateType>
      void MechanicsNonlinearFEOperator :: updateMaterialForElement(const
          AMP::Mesh::MeshManager::Adapter::Element & elem, const std::vector<AMP::Mesh::DOFMap::shared_ptr> & dof_maps)
      {
        std::vector<std::vector<double> > elementInputVectors1(Mechanics::TOTAL_NUMBER_OF_VARIABLES);
        std::vector<std::vector<double> > elementInputVectors_pre1(Mechanics::TOTAL_NUMBER_OF_VARIABLES);

        updateMaterialForElementCommonFunction(elem, dof_maps, elementInputVectors1, elementInputVectors_pre1);

        d_mechNonlinElem->updateMaterialModel<updateType>(elementInputVectors1);
      }

    template <MechanicsNonlinearUpdatedLagrangianElement::MaterialUpdateType updateType>
      void MechanicsNonlinearFEOperator :: updateMaterialForUpdatedLagrangianElement(const
          AMP::Mesh::MeshManager::Adapter::Element & elem, const std::vector<AMP::Mesh::DOFMap::shared_ptr> & dof_maps)
      {
        std::vector<std::vector<double> > elementInputVectors2(Mechanics::TOTAL_NUMBER_OF_VARIABLES);
        std::vector<std::vector<double> > elementInputVectors_pre2(Mechanics::TOTAL_NUMBER_OF_VARIABLES);

        updateMaterialForElementCommonFunction(elem, dof_maps, elementInputVectors2, elementInputVectors_pre2);

        d_mechNULElem->updateMaterialModel<updateType>(elementInputVectors2, elementInputVectors_pre2);
      }
  }
}

#endif



