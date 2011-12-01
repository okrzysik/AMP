
#include "MechanicsNonlinearFEOperator.h"
#include "utils/Utilities.h"
#include "utils/InputDatabase.h"

namespace AMP {
  namespace Operator {

    MechanicsNonlinearFEOperator :: MechanicsNonlinearFEOperator (
        const boost::shared_ptr<MechanicsNonlinearFEOperatorParameters> & params)
      : NonlinearFEOperator (params) {
        AMP_INSIST( ((params.get()) != NULL), "NULL parameter!" );
        AMP_INSIST( (((params->d_db).get()) != NULL), "NULL database!" );

        d_resetReusesRadialReturn = (params->d_db)->getBoolWithDefault("RESET_REUSES_RADIAL_RETURN", true);
        d_jacobianReusesRadialReturn = (params->d_db)->getBoolWithDefault("JACOBIAN_REUSES_RADIAL_RETURN", true);

        d_useUpdatedLagrangian = (params->d_db)->getBoolWithDefault("USE_UPDATED_LAGRANGIAN", false);

        if(d_useUpdatedLagrangian) {
          d_mechNULElem = boost::dynamic_pointer_cast<MechanicsNonlinearUpdatedLagrangianElement>(d_elemOp);
        } else {
          d_mechNonlinElem = boost::dynamic_pointer_cast<MechanicsNonlinearElement>(d_elemOp);
        }

        if(d_useUpdatedLagrangian) {
          AMP_INSIST( ((d_mechNULElem.get()) != NULL), "d_elemOp is not of type MechanicsNonlinearUpdatedLagrangianElement" );
        } else {
          AMP_INSIST( ((d_mechNonlinElem.get()) != NULL), "d_elemOp is not of type MechanicsNonlinearElement" );
        }

        d_materialModel = params->d_materialModel;

        d_isActive.resize(Mechanics::TOTAL_NUMBER_OF_VARIABLES);
        d_isFrozen.resize(Mechanics::TOTAL_NUMBER_OF_VARIABLES);
        d_inVec.resize(Mechanics::TOTAL_NUMBER_OF_VARIABLES);
        if(d_useUpdatedLagrangian) {
          d_inVec_pre.resize(Mechanics::TOTAL_NUMBER_OF_VARIABLES); // storage for variables at the previous config
        }

        AMP_INSIST( params->d_db->keyExists("ActiveInputVariables"), "key not found" );
        boost::shared_ptr<AMP::Database> activeInpVar_db = params->d_db->getDatabase("ActiveInputVariables");

        AMP_INSIST(activeInpVar_db->keyExists("DISPLACEMENT"), "DISPLACEMENT must be active");

        d_isActive[Mechanics::DISPLACEMENT] = true;
        d_isActive[Mechanics::TEMPERATURE] = activeInpVar_db->keyExists("TEMPERATURE");
        d_isActive[Mechanics::BURNUP] = activeInpVar_db->keyExists("BURNUP");
        d_isActive[Mechanics::OXYGEN_CONCENTRATION] = activeInpVar_db->keyExists("OXYGEN_CONCENTRATION");
        d_isActive[Mechanics::LHGR] = activeInpVar_db->keyExists("LHGR");

        d_isFrozen[Mechanics::DISPLACEMENT] = false;
        d_isFrozen[Mechanics::TEMPERATURE] = (params->d_db)->getBoolWithDefault("FREEZE_TEMPERATURE", true);
        d_isFrozen[Mechanics::BURNUP] = (params->d_db)->getBoolWithDefault("FREEZE_BURNUP", true);
        d_isFrozen[Mechanics::OXYGEN_CONCENTRATION] = (params->d_db)->getBoolWithDefault("FREEZE_OXYGEN_CONCENTRATION", true);
        d_isFrozen[Mechanics::LHGR] = (params->d_db)->getBoolWithDefault("FREEZE_LHGR", true);

        d_inpVariables.reset(new AMP::LinearAlgebra::MultiVariable("myInpVar"));
        for(unsigned int i = 0; i < Mechanics::TOTAL_NUMBER_OF_VARIABLES; i++) {
          AMP::LinearAlgebra::Variable::shared_ptr dummyVar;
          d_inpVariables->add(dummyVar);
        }//end for i

        std::string dispVarName = activeInpVar_db->getString("DISPLACEMENT");
        AMP::LinearAlgebra::Variable::shared_ptr dispVar(new AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable,
            3>(dispVarName, d_MeshAdapter) ); 
        d_inpVariables->setVariable(Mechanics::DISPLACEMENT, dispVar);

        if(d_isActive[Mechanics::TEMPERATURE]) {
          std::string varName = activeInpVar_db->getString("TEMPERATURE");
          AMP::LinearAlgebra::Variable::shared_ptr dummyVar(new AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable, 1>(varName, d_MeshAdapter) );
          d_inpVariables->setVariable(Mechanics::TEMPERATURE, dummyVar);
          if(d_isFrozen[Mechanics::TEMPERATURE]) {
            if( params->d_FrozenTemperature != NULL ) {
              setVector(Mechanics::TEMPERATURE, params->d_FrozenTemperature);
            }
          }
        }

        if(d_isActive[Mechanics::BURNUP]) {
          std::string varName = activeInpVar_db->getString("BURNUP");
          AMP::LinearAlgebra::Variable::shared_ptr dummyVar(new AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable, 1>(varName, d_MeshAdapter) );
          d_inpVariables->setVariable(Mechanics::BURNUP, dummyVar);
          if(d_isFrozen[Mechanics::BURNUP]) {
            if( params->d_FrozenBurnup != NULL ) {
              setVector(Mechanics::BURNUP, params->d_FrozenBurnup);
            }
          }
        }

        if(d_isActive[Mechanics::OXYGEN_CONCENTRATION]) {
          std::string varName = activeInpVar_db->getString("OXYGEN_CONCENTRATION");
          AMP::LinearAlgebra::Variable::shared_ptr dummyVar(new AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable, 1>(varName, d_MeshAdapter) );
          d_inpVariables->setVariable(Mechanics::OXYGEN_CONCENTRATION, dummyVar);
          if(d_isFrozen[Mechanics::OXYGEN_CONCENTRATION]) {
            if( params->d_FrozenOxygenConcentration != NULL ) {
              setVector(Mechanics::OXYGEN_CONCENTRATION, params->d_FrozenOxygenConcentration);
            }
          }
        }

        if(d_isActive[Mechanics::LHGR]) {
          std::string varName = activeInpVar_db->getString("LHGR");
          AMP::LinearAlgebra::Variable::shared_ptr dummyVar(new AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable, 1>(varName, d_MeshAdapter) );
          d_inpVariables->setVariable(Mechanics::LHGR, dummyVar);
          if(d_isFrozen[Mechanics::LHGR]) {
            if( params->d_FrozenLHGR != NULL ) {
              setVector(Mechanics::LHGR, params->d_FrozenLHGR);
            }
          }
        }

        //memory for reference coordinates (used in UL formulation)
        if(d_useUpdatedLagrangian) {
          d_refXYZ = d_MeshAdapter->createVector(d_inpVariables->getVariable(Mechanics::DISPLACEMENT));
          d_refXYZ->zero();
        }

        // memory for variables in the previous config
        if(d_useUpdatedLagrangian) {
          for(unsigned int i=0; i < Mechanics::TOTAL_NUMBER_OF_VARIABLES; i++)
          {
            if(d_isActive[i]) {
              d_inVec_pre[i] = d_MeshAdapter->createVector(d_inpVariables->getVariable(i)); 
              d_inVec_pre[i]->zero();
            }
          }
        }

        if(params->d_ReferenceTemperature != NULL) {
          setReferenceTemperature(params->d_ReferenceTemperature);
        }

        AMP_INSIST( params->d_db->keyExists("OutputVariable"), "key not found" );
        std::string outVarName = params->d_db->getString("OutputVariable");
        d_outVariable.reset(new AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable, 3>(outVarName, d_MeshAdapter));

        d_isInitialized = false;
      }

    unsigned int MechanicsNonlinearFEOperator :: numberOfDOFMaps() {
      if( d_isActive[Mechanics::TEMPERATURE] || d_isActive[Mechanics::BURNUP] || 
          d_isActive[Mechanics::OXYGEN_CONCENTRATION] || d_isActive[Mechanics::LHGR] ) {
        return 2;
      } else {
        return 1; 
      }
    }

    AMP::LinearAlgebra::Variable::shared_ptr MechanicsNonlinearFEOperator :: getVariableForDOFMap(unsigned int id) {
      AMP_ASSERT( id < (this->numberOfDOFMaps()) );
      if(id == 0) {
        return (d_inpVariables->getVariable(Mechanics::DISPLACEMENT));
      } else {
        if(d_isActive[Mechanics::TEMPERATURE]) {
          return (d_inpVariables->getVariable(Mechanics::TEMPERATURE));
        } else if(d_isActive[Mechanics::BURNUP]) {
          return (d_inpVariables->getVariable(Mechanics::BURNUP));
        } else if(d_isActive[Mechanics::OXYGEN_CONCENTRATION]) {
          return (d_inpVariables->getVariable(Mechanics::OXYGEN_CONCENTRATION));
        } else {
          AMP_ASSERT(d_isActive[Mechanics::LHGR]);
          return (d_inpVariables->getVariable(Mechanics::LHGR));
        }
      }
    }

    AMP::LinearAlgebra::Variable::shared_ptr MechanicsNonlinearFEOperator :: createInputVariable(const std::string & name, int varId) {
      AMP::LinearAlgebra::Variable::shared_ptr inpVar;
      switch(varId) {
        case Mechanics::DISPLACEMENT : {
                                         inpVar.reset(new AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable, 3>(name) );
                                         break;
                                       }
        case Mechanics::TEMPERATURE : {
                                        inpVar.reset(new AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable, 1>(name) );
                                        break;
                                      }
        case Mechanics::BURNUP : {
                                   inpVar.reset(new AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable, 1>(name) );
                                   break;
                                 }
        case Mechanics::OXYGEN_CONCENTRATION : {
                                                 inpVar.reset(new AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable, 1>(name) );
                                                 break;
                                               }
        case Mechanics::LHGR : {
                                 inpVar.reset(new AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable, 1>(name) );
                                 break;
                               }
        default: 
                               assert(false);
      }
      return inpVar;
    }

    void MechanicsNonlinearFEOperator :: preAssembly(const boost::shared_ptr< AMP::LinearAlgebra::Vector >  &u, 
        boost::shared_ptr< AMP::LinearAlgebra::Vector >  &r) {

      AMP_INSIST( (u != NULL), "NULL Input Vector" );

      if(!d_isInitialized) {
        init();
      }

      AMP::LinearAlgebra::Variable::shared_ptr dispVar = d_inpVariables->getVariable(Mechanics::DISPLACEMENT);
      AMP::LinearAlgebra::Vector::shared_ptr dispVector = u->subsetVectorForVariable(dispVar);
      setVector(Mechanics::DISPLACEMENT, dispVector);

      if(d_isActive[Mechanics::TEMPERATURE]) {
        if(!(d_isFrozen[Mechanics::TEMPERATURE])) {
          AMP::LinearAlgebra::Variable::shared_ptr tempVar = d_inpVariables->getVariable(Mechanics::TEMPERATURE); 
          AMP::LinearAlgebra::Vector::shared_ptr tempVector = u->subsetVectorForVariable(tempVar);
          setVector(Mechanics::TEMPERATURE, tempVector);
        }
      }

      if(d_isActive[Mechanics::BURNUP]) {
        if(!(d_isFrozen[Mechanics::BURNUP])) {
          AMP::LinearAlgebra::Variable::shared_ptr burnVar = d_inpVariables->getVariable(Mechanics::BURNUP);
          AMP::LinearAlgebra::Vector::shared_ptr burnVector = u->subsetVectorForVariable(burnVar);
          setVector(Mechanics::BURNUP, burnVector);
        }
      }

      if(d_isActive[Mechanics::OXYGEN_CONCENTRATION]) {
        if(!(d_isFrozen[Mechanics::OXYGEN_CONCENTRATION])) {
          AMP::LinearAlgebra::Variable::shared_ptr oxyVar = d_inpVariables->getVariable(Mechanics::OXYGEN_CONCENTRATION); 
          AMP::LinearAlgebra::Vector::shared_ptr oxyVector = u->subsetVectorForVariable(oxyVar);
          setVector(Mechanics::OXYGEN_CONCENTRATION, oxyVector);
        }
      }

      if(d_isActive[Mechanics::LHGR]) {
        if(!(d_isFrozen[Mechanics::LHGR])) {
          AMP::LinearAlgebra::Variable::shared_ptr lhgrVar = d_inpVariables->getVariable(Mechanics::LHGR);
          AMP::LinearAlgebra::Vector::shared_ptr lhgrVector = u->subsetVectorForVariable(lhgrVar);
          setVector(Mechanics::LHGR, lhgrVector);
        }
      }

      d_outVec = r->subsetVectorForVariable(d_outVariable);
      d_outVec->zero();

      d_materialModel->preNonlinearAssembly();

      if(d_useUpdatedLagrangian == true) {
        d_mechNULElem->zeroOutGaussPointCount();
      }
    }

    void MechanicsNonlinearFEOperator :: postAssembly()
    {
      d_materialModel->postNonlinearAssembly();

      d_outVec->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_ADD );
    }

    void MechanicsNonlinearFEOperator :: preElementOperation( const AMP::Mesh::MeshManager::Adapter::Element & elem, 
        const std::vector<AMP::Mesh::DOFMap::shared_ptr> & dof_maps )
    {
      unsigned int num_local_type0Dofs = 0;
      for(unsigned int i = 0; i < 3; i++) {
        (dof_maps[0])->getDOFs (elem, d_type0DofIndices[i], i);
        num_local_type0Dofs += d_type0DofIndices[i].size();
      }//end for i

      unsigned int num_local_type1Dofs = 0;
      if( d_isActive[Mechanics::TEMPERATURE] || d_isActive[Mechanics::BURNUP] ||
          d_isActive[Mechanics::OXYGEN_CONCENTRATION] || d_isActive[Mechanics::LHGR] ) {
        (dof_maps[1])->getDOFs (elem, d_type1DofIndices);
        num_local_type1Dofs = d_type1DofIndices.size();
      }

      std::vector<std::vector<double> > elementInputVectors(Mechanics::TOTAL_NUMBER_OF_VARIABLES);
      std::vector<std::vector<double> > elementInputVectors_pre(Mechanics::TOTAL_NUMBER_OF_VARIABLES);

      elementInputVectors[Mechanics::DISPLACEMENT].resize(num_local_type0Dofs);
      if(d_useUpdatedLagrangian) {
        elementInputVectors_pre[Mechanics::DISPLACEMENT].resize(num_local_type0Dofs);
      }
      if(d_isActive[Mechanics::TEMPERATURE]) {
        elementInputVectors[Mechanics::TEMPERATURE].resize(num_local_type1Dofs);
        if(d_useUpdatedLagrangian) {
          elementInputVectors_pre[Mechanics::TEMPERATURE].resize(num_local_type1Dofs);
        }
      }
      if(d_isActive[Mechanics::BURNUP]) {
        elementInputVectors[Mechanics::BURNUP].resize(num_local_type1Dofs);
        if(d_useUpdatedLagrangian) {
          elementInputVectors_pre[Mechanics::BURNUP].resize(num_local_type1Dofs);
        }
      }
      if(d_isActive[Mechanics::OXYGEN_CONCENTRATION]) {
        elementInputVectors[Mechanics::OXYGEN_CONCENTRATION].resize(num_local_type1Dofs);
        if(d_useUpdatedLagrangian) {
          elementInputVectors_pre[Mechanics::OXYGEN_CONCENTRATION].resize(num_local_type1Dofs);
        }
      }
      if(d_isActive[Mechanics::LHGR]) {
        elementInputVectors[Mechanics::LHGR].resize(num_local_type1Dofs);
        if(d_useUpdatedLagrangian) {
          elementInputVectors_pre[Mechanics::LHGR].resize(num_local_type1Dofs);
        }
      }

      d_numNodesForCurrentElement = elem.numNodes(); 

      std::vector<double> elementRefXYZ;
      elementRefXYZ.resize(3 * d_numNodesForCurrentElement);

      for(unsigned int r = 0; r < d_numNodesForCurrentElement; r++) {
        for(unsigned int d = 0; d < 3; d++) {
          elementInputVectors[Mechanics::DISPLACEMENT][(3*r) + d] = (d_inVec[Mechanics::DISPLACEMENT])->
            getValueByGlobalID( d_type0DofIndices[d][r] );
          if(d_useUpdatedLagrangian) {
            elementInputVectors_pre[Mechanics::DISPLACEMENT][(3*r) + d] = (d_inVec_pre[Mechanics::DISPLACEMENT])->getValueByGlobalID( d_type0DofIndices[d][r] );
            elementRefXYZ[(3 * r) + d] = d_refXYZ->getValueByGlobalID(d_type0DofIndices[d][r]);
            //AMP::pout<<"elementRefXYZ["<<(3 * r) + d<<"] = "<<elementRefXYZ[(3 * r) + d]<<std::endl;
          }
        }
        if(d_isActive[Mechanics::TEMPERATURE]) {
          elementInputVectors[Mechanics::TEMPERATURE][r] = (d_inVec[Mechanics::TEMPERATURE])->
            getValueByGlobalID( d_type1DofIndices[r] );
          if(d_useUpdatedLagrangian) {
            elementInputVectors_pre[Mechanics::TEMPERATURE][r] = (d_inVec_pre[Mechanics::TEMPERATURE])->getValueByGlobalID( d_type1DofIndices[r] );
          }
        }
        if(d_isActive[Mechanics::BURNUP]) {
          elementInputVectors[Mechanics::BURNUP][r] = (d_inVec[Mechanics::BURNUP])->
            getValueByGlobalID( d_type1DofIndices[r] );
          if(d_useUpdatedLagrangian) {
            elementInputVectors_pre[Mechanics::BURNUP][r] = (d_inVec_pre[Mechanics::BURNUP])->getValueByGlobalID( d_type1DofIndices[r] );
          }
        }
        if(d_isActive[Mechanics::OXYGEN_CONCENTRATION]) {
          elementInputVectors[Mechanics::OXYGEN_CONCENTRATION][r] = (d_inVec[Mechanics::OXYGEN_CONCENTRATION])->
            getValueByGlobalID( d_type1DofIndices[r] );
          if(d_useUpdatedLagrangian) {
            elementInputVectors_pre[Mechanics::OXYGEN_CONCENTRATION][r] = (d_inVec_pre[Mechanics::OXYGEN_CONCENTRATION])->getValueByGlobalID( d_type1DofIndices[r] );
          }
        }
        if(d_isActive[Mechanics::LHGR]) {
          elementInputVectors[Mechanics::LHGR][r] = (d_inVec[Mechanics::LHGR])->
            getValueByGlobalID( d_type1DofIndices[r] );
          if(d_useUpdatedLagrangian) {
            elementInputVectors_pre[Mechanics::LHGR][r] = (d_inVec_pre[Mechanics::LHGR])->getValueByGlobalID( d_type1DofIndices[r] );
          }
        }
      }

      d_elementOutputVector.resize(num_local_type0Dofs);
      for(unsigned int i = 0; i < num_local_type0Dofs; i++) {
        d_elementOutputVector[i] = 0.0;
      }

      const ::Elem* elemPtr = &(elem.getElem());

      if(!d_useUpdatedLagrangian) {
        d_mechNonlinElem->initializeForCurrentElement( elemPtr, d_materialModel );
        d_mechNonlinElem->setElementVectors( elementInputVectors, d_elementOutputVector );
      } else {
        d_mechNULElem->initializeForCurrentElement( elemPtr, d_materialModel );
        d_mechNULElem->setElementVectors( elementInputVectors, elementInputVectors_pre, d_elementOutputVector );
        d_mechNULElem->assignReferenceXYZ(elementRefXYZ);
      }
    }

    void MechanicsNonlinearFEOperator :: postElementOperation()
    {
      for(unsigned int r = 0; r < d_numNodesForCurrentElement; r++) {
        for(unsigned int d = 0; d < 3; d++) {
          d_outVec->addValueByGlobalID( d_type0DofIndices[d][r], d_elementOutputVector[(3*r) + d] );
        }
      }
    }

    void MechanicsNonlinearFEOperator :: init() 
    {
      d_isInitialized = true;

      AMP::Mesh::DOFMap::shared_ptr dof_map;
      if(d_isActive[Mechanics::TEMPERATURE]) {
        dof_map = d_MeshAdapter->getDOFMap(d_inpVariables->getVariable(Mechanics::TEMPERATURE)); 
      }

      AMP::Mesh::DOFMap::shared_ptr dof_map_disp;
      if(d_useUpdatedLagrangian) {
        dof_map_disp = d_MeshAdapter->getDOFMap(d_inpVariables->getVariable(Mechanics::DISPLACEMENT)); 
      }

      AMP::Mesh::MeshManager::Adapter::ElementIterator  el = d_MeshAdapter->beginElement();
      AMP::Mesh::MeshManager::Adapter::ElementIterator  end_el = d_MeshAdapter->endElement();

      d_materialModel->preNonlinearInit(d_resetReusesRadialReturn, d_jacobianReusesRadialReturn);

      if(d_useUpdatedLagrangian) {
        d_mechNULElem->preNonlinearElementInit();
      }

      for( ; el != end_el; ++el) {
        std::vector<unsigned int> dofIndices;
        if(d_isActive[Mechanics::TEMPERATURE]) {
          dof_map->getDOFs (*el, dofIndices);
        }
        unsigned int numDofs = dofIndices.size();

        std::vector<unsigned int> dofIndices_disp[3];
        unsigned int num_dofIndices_disp = 0;
        if(d_useUpdatedLagrangian) {
          for(unsigned int i = 0; i < 3; i++) {
            (dof_map_disp)->getDOFs (*el, dofIndices_disp[i], i);
            num_dofIndices_disp += dofIndices_disp[i].size();
          }//end for i
        }

        std::vector<double> localVector;
        if(d_isActive[Mechanics::TEMPERATURE]) {
          localVector.resize(numDofs);
          for(unsigned int r = 0; r < numDofs; r++) {
            localVector[r] = d_referenceTemperature->getValueByGlobalID(dofIndices[r]);
          }
        }

        std::vector<double> elementRefXYZ;
        if(d_useUpdatedLagrangian) {
          elementRefXYZ.resize(num_dofIndices_disp);
        }

        const ::Elem* elemPtr = &(el->getElem());

        if(d_useUpdatedLagrangian) {
          d_mechNULElem->initializeForCurrentElement( elemPtr, d_materialModel );
          d_mechNULElem->initMaterialModel(localVector);
          d_mechNULElem->initializeReferenceXYZ(elementRefXYZ);
        } else {
          d_mechNonlinElem->initializeForCurrentElement( elemPtr, d_materialModel );
          d_mechNonlinElem->initMaterialModel(localVector);
        }

        if(d_useUpdatedLagrangian) {
          for(unsigned int i = 0; i < 3; i++) {
            for(unsigned int j = 0; j < dofIndices_disp[i].size(); j++) {
              //AMP::pout<<"elementRefXYZ["<<(3 * j) + i<<"] = "<<elementRefXYZ[(3 * j) + i]<<std::endl;
              d_refXYZ->setValueByGlobalID(dofIndices_disp[i][j], elementRefXYZ[(3 * j) + i]);
            }
          }
        }

      }//end for el

      if(d_useUpdatedLagrangian) {
        d_refXYZ->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
      }

      d_materialModel->postNonlinearInit();
    }

    void MechanicsNonlinearFEOperator :: reset(const boost::shared_ptr<OperatorParameters>& params)
    {
      if(!d_isInitialized) {
        init();
      }

      boost::shared_ptr<MechanicsNonlinearFEOperatorParameters> myParams =
        boost::dynamic_pointer_cast<MechanicsNonlinearFEOperatorParameters>(params); 

      AMP_INSIST( ((myParams.get()) != NULL), "Null parameter!" );

      if(d_resetReusesRadialReturn) {
        d_materialModel->globalReset();
      } else {
        unsigned int numDOFMaps = numberOfDOFMaps();
        std::vector<AMP::Mesh::DOFMap::shared_ptr> dof_maps(numDOFMaps);

        for(unsigned int i = 0; i < numDOFMaps; i++) {
          dof_maps[i] = d_MeshAdapter->getDOFMap( getVariableForDOFMap(i) );
        }

        AMP::Mesh::MeshManager::Adapter::ElementIterator  el = d_MeshAdapter->beginElement();
        AMP::Mesh::MeshManager::Adapter::ElementIterator  end_el = d_MeshAdapter->endElement();

        setVector(Mechanics::DISPLACEMENT, myParams->d_EquilibriumDisplacement);

        if(d_isActive[Mechanics::TEMPERATURE]) {
          if( myParams->d_EquilibriumTemperature != NULL ) {
            setVector(Mechanics::TEMPERATURE, myParams->d_EquilibriumTemperature);
          }
        }

        if(d_isActive[Mechanics::BURNUP]) {
          if( myParams->d_EquilibriumBurnup != NULL ) {
            setVector(Mechanics::BURNUP, myParams->d_EquilibriumBurnup);
          }
        }

        if(d_isActive[Mechanics::OXYGEN_CONCENTRATION]) {
          if( myParams->d_EquilibriumOxygenConcentration != NULL ) {
            setVector(Mechanics::OXYGEN_CONCENTRATION, myParams->d_EquilibriumOxygenConcentration);
          }
        }

        if(d_isActive[Mechanics::LHGR]) {
          if( myParams->d_EquilibriumLHGR != NULL ) {
            setVector(Mechanics::LHGR, myParams->d_EquilibriumLHGR);
          }
        }

        d_outVec.reset();

        d_materialModel->preNonlinearReset();

        for( ; el != end_el; ++el) {
          if(d_useUpdatedLagrangian) {
            updateMaterialForUpdatedLagrangianElement<MechanicsNonlinearUpdatedLagrangianElement::RESET>(*el, dof_maps);
          } else {
            updateMaterialForElement<MechanicsNonlinearElement::RESET>(*el, dof_maps);
          }
        }//end for el

        d_materialModel->postNonlinearReset();
      }

      if(d_useUpdatedLagrangian == true) {
        d_mechNULElem->resetElementInfo();
      }

      if(d_useUpdatedLagrangian) {
        //Copy values in pre
        for(unsigned int i = 0; i < Mechanics::TOTAL_NUMBER_OF_VARIABLES; i++) {
          if(d_isActive[i]) {
            d_inVec_pre[i]->copyVector(d_inVec[i]);
            d_inVec_pre[i]->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
          }
        }
      }

      //
      if(d_isActive[Mechanics::TEMPERATURE]) {
        if(d_isFrozen[Mechanics::TEMPERATURE]) {
          if( myParams->d_FrozenTemperature != NULL ) {
            setVector(Mechanics::TEMPERATURE, myParams->d_FrozenTemperature);
          }
        }
      }

      if(d_isActive[Mechanics::BURNUP]) {
        if(d_isFrozen[Mechanics::BURNUP]) {
          if( myParams->d_FrozenBurnup != NULL ) {
            setVector(Mechanics::BURNUP, myParams->d_FrozenBurnup);
          }
        }
      }

      if(d_isActive[Mechanics::OXYGEN_CONCENTRATION]) {
        if(d_isFrozen[Mechanics::OXYGEN_CONCENTRATION]) {
          if( myParams->d_FrozenOxygenConcentration != NULL ) {
            setVector(Mechanics::OXYGEN_CONCENTRATION, myParams->d_FrozenOxygenConcentration);
          }
        }
      }

      if(d_isActive[Mechanics::LHGR]) {
        if(d_isFrozen[Mechanics::LHGR]) {
          if( myParams->d_FrozenLHGR != NULL ) {
            setVector(Mechanics::LHGR, myParams->d_FrozenLHGR);
          }
        }
      }

    }

    boost::shared_ptr<OperatorParameters> MechanicsNonlinearFEOperator ::
      getJacobianParameters(const boost::shared_ptr<AMP::LinearAlgebra::Vector>& u) {
        if(!d_isInitialized) {
          init();
        }

        // set up a database for the linear operator params
        boost::shared_ptr<AMP::InputDatabase> tmp_db (new AMP::InputDatabase("Dummy"));
        tmp_db->putBool("reset_reuses_matrix", true);
        tmp_db->putBool("isAttachedToNonlinearOperator", true);
        tmp_db->putBool("isNonlinearOperatorInitialized", true);

        // create the linear operator params
        boost::shared_ptr<MechanicsLinearFEOperatorParameters> outParams(new
            MechanicsLinearFEOperatorParameters(tmp_db));

        // If updated-lagrangian is being used, then displacement and other variable has to be passed to the linear element level.
        if(d_useUpdatedLagrangian) {
          AMP::LinearAlgebra::Vector::shared_ptr displacementVector = u->subsetVectorForVariable(d_inpVariables->getVariable(Mechanics::DISPLACEMENT));
          outParams->d_dispVec = displacementVector;
          (outParams->d_dispVec)->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
        }

        if(d_jacobianReusesRadialReturn == false) {
          unsigned int numDOFMaps = numberOfDOFMaps();
          std::vector<AMP::Mesh::DOFMap::shared_ptr> dof_maps(numDOFMaps);

          for(unsigned int i = 0; i < numDOFMaps; i++) {
            dof_maps[i] = d_MeshAdapter->getDOFMap( getVariableForDOFMap(i) );
          }

          AMP::LinearAlgebra::Vector::shared_ptr dispVector = u->subsetVectorForVariable (
              d_inpVariables->getVariable(Mechanics::DISPLACEMENT) );
          setVector(Mechanics::DISPLACEMENT, dispVector);

          if(d_isActive[Mechanics::TEMPERATURE]) {
            if(!(d_isFrozen[Mechanics::TEMPERATURE])) {
              AMP::LinearAlgebra::Vector::shared_ptr tempVector = u->subsetVectorForVariable (
                  d_inpVariables->getVariable(Mechanics::TEMPERATURE) );
              setVector(Mechanics::TEMPERATURE, tempVector);
            }
          }

          if(d_isActive[Mechanics::BURNUP]) {
            if(!(d_isFrozen[Mechanics::BURNUP])) {
              AMP::LinearAlgebra::Vector::shared_ptr burnVector = u->subsetVectorForVariable (
                  d_inpVariables->getVariable(Mechanics::BURNUP) );
              setVector(Mechanics::BURNUP, burnVector);
            }
          }

          if(d_isActive[Mechanics::OXYGEN_CONCENTRATION]) {
            if(!(d_isFrozen[Mechanics::OXYGEN_CONCENTRATION])) {
              AMP::LinearAlgebra::Vector::shared_ptr oxyVector = u->subsetVectorForVariable (
                  d_inpVariables->getVariable(Mechanics::OXYGEN_CONCENTRATION) );
              setVector(Mechanics::OXYGEN_CONCENTRATION, oxyVector);
            }
          }

          if(d_isActive[Mechanics::LHGR]) {
            if(!(d_isFrozen[Mechanics::LHGR])) {
              AMP::LinearAlgebra::Variable::shared_ptr lhgrVar = d_inpVariables->getVariable(Mechanics::LHGR);
              AMP::LinearAlgebra::Vector::shared_ptr lhgrVector = u->subsetVectorForVariable(lhgrVar);
              setVector(Mechanics::LHGR, lhgrVector);
            }
          }

          d_outVec.reset();

          d_materialModel->preNonlinearJacobian();

          AMP::Mesh::MeshManager::Adapter::ElementIterator  el = d_MeshAdapter->beginElement();
          AMP::Mesh::MeshManager::Adapter::ElementIterator  end_el = d_MeshAdapter->endElement();

          for( ; el != end_el; ++el) {
            if(d_useUpdatedLagrangian) {
              updateMaterialForUpdatedLagrangianElement<MechanicsNonlinearUpdatedLagrangianElement::JACOBIAN>(*el, dof_maps);
            } else {
              updateMaterialForElement<MechanicsNonlinearElement::JACOBIAN>(*el, dof_maps);
            }
          }//end for el

          d_materialModel->postNonlinearJacobian();
        }

        return outParams;
      }

    void MechanicsNonlinearFEOperator :: printStressAndStrain(AMP::LinearAlgebra::Vector::shared_ptr u,
        const std::string & fname) {
      if(!d_isInitialized) {
        init();
      }

      unsigned int numDOFMaps = this->numberOfDOFMaps();
      std::vector<AMP::Mesh::DOFMap::shared_ptr> dof_maps(numDOFMaps);

      for(unsigned int i = 0; i < numDOFMaps; i++) {
        dof_maps[i] = d_MeshAdapter->getDOFMap( this->getVariableForDOFMap(i) );
      }

      AMP::Mesh::MeshManager::Adapter::ElementIterator  el = d_MeshAdapter->beginElement();
      AMP::Mesh::MeshManager::Adapter::ElementIterator  end_el = d_MeshAdapter->endElement();

      FILE* fp = fopen(fname.c_str(), "w");

      //fprintf(fp, "x, y, z, Stresses(11, 22, 33, 23, 13, 12), Strains(11, 22, 33, 23, 13, 12) \n\n");

      AMP::LinearAlgebra::Variable::shared_ptr dispVar = d_inpVariables->getVariable(Mechanics::DISPLACEMENT);
      AMP::LinearAlgebra::Vector::shared_ptr dispVector = u->subsetVectorForVariable(dispVar);
      setVector(Mechanics::DISPLACEMENT, dispVector);

      if(d_isActive[Mechanics::TEMPERATURE]) {
        if(!(d_isFrozen[Mechanics::TEMPERATURE])) {
          AMP::LinearAlgebra::Variable::shared_ptr tempVar = d_inpVariables->getVariable(Mechanics::TEMPERATURE); 
          AMP::LinearAlgebra::Vector::shared_ptr tempVector = u->subsetVectorForVariable(tempVar);
          setVector(Mechanics::TEMPERATURE, tempVector);
        }
      }

      if(d_isActive[Mechanics::BURNUP]) {
        if(!(d_isFrozen[Mechanics::BURNUP])) {
          AMP::LinearAlgebra::Variable::shared_ptr burnVar = d_inpVariables->getVariable(Mechanics::BURNUP);
          AMP::LinearAlgebra::Vector::shared_ptr burnVector = u->subsetVectorForVariable(burnVar);
          setVector(Mechanics::BURNUP, burnVector);
        }
      }

      if(d_isActive[Mechanics::OXYGEN_CONCENTRATION]) {
        if(!(d_isFrozen[Mechanics::OXYGEN_CONCENTRATION])) {
          AMP::LinearAlgebra::Variable::shared_ptr oxyVar = d_inpVariables->getVariable(Mechanics::OXYGEN_CONCENTRATION); 
          AMP::LinearAlgebra::Vector::shared_ptr oxyVector = u->subsetVectorForVariable(oxyVar);
          setVector(Mechanics::OXYGEN_CONCENTRATION, oxyVector);
        }
      }

      if(d_isActive[Mechanics::LHGR]) {
        if(!(d_isFrozen[Mechanics::LHGR])) {
          AMP::LinearAlgebra::Variable::shared_ptr lhgrVar = d_inpVariables->getVariable(Mechanics::LHGR);
          AMP::LinearAlgebra::Vector::shared_ptr lhgrVector = u->subsetVectorForVariable(lhgrVar);
          setVector(Mechanics::LHGR, lhgrVector);
        }
      }

      d_materialModel->preNonlinearAssembly();

      for( ; el != end_el; ++el) {
        unsigned int num_local_type0Dofs = 0;
        for(unsigned int i = 0; i < 3; i++) {
          (dof_maps[0])->getDOFs (*el, d_type0DofIndices[i], i);
          num_local_type0Dofs += d_type0DofIndices[i].size();
        }//end for i

        unsigned int num_local_type1Dofs = 0;
        if( d_isActive[Mechanics::TEMPERATURE] || d_isActive[Mechanics::BURNUP] ||
            d_isActive[Mechanics::OXYGEN_CONCENTRATION] || d_isActive[Mechanics::LHGR]) {
          (dof_maps[1])->getDOFs (*el, d_type1DofIndices);
          num_local_type1Dofs = d_type1DofIndices.size();
        }

        std::vector<std::vector<double> > elementInputVectors(Mechanics::TOTAL_NUMBER_OF_VARIABLES);
        std::vector<std::vector<double> > elementInputVectors_pre(Mechanics::TOTAL_NUMBER_OF_VARIABLES);

        elementInputVectors[Mechanics::DISPLACEMENT].resize(num_local_type0Dofs);
        if(d_useUpdatedLagrangian) {
          elementInputVectors_pre[Mechanics::DISPLACEMENT].resize(num_local_type0Dofs);
        }
        if(d_isActive[Mechanics::TEMPERATURE]) {
          elementInputVectors[Mechanics::TEMPERATURE].resize(num_local_type1Dofs);
          if(d_useUpdatedLagrangian) {
            elementInputVectors_pre[Mechanics::TEMPERATURE].resize(num_local_type1Dofs);
          }
        }
        if(d_isActive[Mechanics::BURNUP]) {
          elementInputVectors[Mechanics::BURNUP].resize(num_local_type1Dofs);
          if(d_useUpdatedLagrangian) {
            elementInputVectors_pre[Mechanics::BURNUP].resize(num_local_type1Dofs);
          }
        }
        if(d_isActive[Mechanics::OXYGEN_CONCENTRATION]) {
          elementInputVectors[Mechanics::OXYGEN_CONCENTRATION].resize(num_local_type1Dofs);
          if(d_useUpdatedLagrangian) {
            elementInputVectors_pre[Mechanics::OXYGEN_CONCENTRATION].resize(num_local_type1Dofs);
          }
        }
        if(d_isActive[Mechanics::LHGR]) {
          elementInputVectors[Mechanics::LHGR].resize(num_local_type1Dofs);
          if(d_useUpdatedLagrangian) {
            elementInputVectors_pre[Mechanics::LHGR].resize(num_local_type1Dofs);
          }
        }

        d_numNodesForCurrentElement = el->numNodes(); 

        for(unsigned int r = 0; r < d_numNodesForCurrentElement; r++) {
          for(unsigned int d = 0; d < 3; d++) {
            elementInputVectors[Mechanics::DISPLACEMENT][(3*r) + d] = (d_inVec[Mechanics::DISPLACEMENT])->
              getValueByGlobalID( d_type0DofIndices[d][r] );
            if(d_useUpdatedLagrangian) {
              elementInputVectors_pre[Mechanics::DISPLACEMENT][(3*r) + d] = (d_inVec_pre[Mechanics::DISPLACEMENT])->getValueByGlobalID( d_type0DofIndices[d][r] );
            }
          }
          if(d_isActive[Mechanics::TEMPERATURE]) {
            elementInputVectors[Mechanics::TEMPERATURE][r] = (d_inVec[Mechanics::TEMPERATURE])->
              getValueByGlobalID( d_type1DofIndices[r] );
            if(d_useUpdatedLagrangian) {
              elementInputVectors_pre[Mechanics::TEMPERATURE][r] = (d_inVec_pre[Mechanics::TEMPERATURE])->getValueByGlobalID( d_type1DofIndices[r] );
            }
          }
          if(d_isActive[Mechanics::BURNUP]) {
            elementInputVectors[Mechanics::BURNUP][r] = (d_inVec[Mechanics::BURNUP])->
              getValueByGlobalID( d_type1DofIndices[r] );
            if(d_useUpdatedLagrangian) {
              elementInputVectors_pre[Mechanics::BURNUP][r] = (d_inVec_pre[Mechanics::BURNUP])->getValueByGlobalID( d_type1DofIndices[r] );
            }
          }
          if(d_isActive[Mechanics::OXYGEN_CONCENTRATION]) {
            elementInputVectors[Mechanics::OXYGEN_CONCENTRATION][r] = (d_inVec[Mechanics::OXYGEN_CONCENTRATION])->
              getValueByGlobalID( d_type1DofIndices[r] );
            if(d_useUpdatedLagrangian) {
              elementInputVectors_pre[Mechanics::OXYGEN_CONCENTRATION][r] = (d_inVec_pre[Mechanics::OXYGEN_CONCENTRATION])->getValueByGlobalID( d_type1DofIndices[r] );
            }
          }
          if(d_isActive[Mechanics::LHGR]) {
            elementInputVectors[Mechanics::LHGR][r] = (d_inVec[Mechanics::LHGR])->
              getValueByGlobalID( d_type1DofIndices[r] );
            if(d_useUpdatedLagrangian) {
              elementInputVectors_pre[Mechanics::LHGR][r] = (d_inVec_pre[Mechanics::LHGR])->getValueByGlobalID( d_type1DofIndices[r] );
            }
          }
        }

        const ::Elem* elemPtr = &(el->getElem());

        d_mechNonlinElem->initializeForCurrentElement( elemPtr, d_materialModel );

        d_mechNonlinElem->printStressAndStrain(fp, elementInputVectors);
      }//end for el

      d_materialModel->postNonlinearAssembly();

      fprintf(fp, "\n\n");

      fclose(fp);
    }

    void MechanicsNonlinearFEOperator :: computeStressesAndStrains(AMP::LinearAlgebra::Vector::shared_ptr u,
        AMP::LinearAlgebra::Vector::shared_ptr & stress, AMP::LinearAlgebra::Vector::shared_ptr & strain) {
      if(!d_isInitialized) {
        init();
      }

      unsigned int numGaussPts = d_mechNonlinElem->getNumberOfGaussPoints();
      AMP::LinearAlgebra::Variable::shared_ptr stressVar ( new AMP::Mesh::RunTimeIntegrationPointVariable ( "stress" , (6*numGaussPts) ) );
      AMP::LinearAlgebra::Variable::shared_ptr strainVar ( new AMP::Mesh::RunTimeIntegrationPointVariable ( "strain" , (6*numGaussPts) ) );

      stress =  d_MeshAdapter->createVector ( stressVar );
      strain =  d_MeshAdapter->createVector ( strainVar );

      AMP::Mesh::DOFMap::shared_ptr gaussPt_dof_map = d_MeshAdapter->getDOFMap(stressVar);

      unsigned int numDOFMaps = this->numberOfDOFMaps();
      std::vector<AMP::Mesh::DOFMap::shared_ptr> dof_maps(numDOFMaps);

      for(unsigned int i = 0; i < numDOFMaps; i++) {
        dof_maps[i] = d_MeshAdapter->getDOFMap( this->getVariableForDOFMap(i) );
      }

      AMP::Mesh::MeshManager::Adapter::ElementIterator  el = d_MeshAdapter->beginElement();
      AMP::Mesh::MeshManager::Adapter::ElementIterator  end_el = d_MeshAdapter->endElement();

      AMP::LinearAlgebra::Variable::shared_ptr dispVar = d_inpVariables->getVariable(Mechanics::DISPLACEMENT);
      AMP::LinearAlgebra::Vector::shared_ptr dispVector = u->subsetVectorForVariable(dispVar);
      setVector(Mechanics::DISPLACEMENT, dispVector);

      if(d_isActive[Mechanics::TEMPERATURE]) {
        if(!(d_isFrozen[Mechanics::TEMPERATURE])) {
          AMP::LinearAlgebra::Variable::shared_ptr tempVar = d_inpVariables->getVariable(Mechanics::TEMPERATURE); 
          AMP::LinearAlgebra::Vector::shared_ptr tempVector = u->subsetVectorForVariable(tempVar);
          setVector(Mechanics::TEMPERATURE, tempVector);
        }
      }

      if(d_isActive[Mechanics::BURNUP]) {
        if(!(d_isFrozen[Mechanics::BURNUP])) {
          AMP::LinearAlgebra::Variable::shared_ptr burnVar = d_inpVariables->getVariable(Mechanics::BURNUP);
          AMP::LinearAlgebra::Vector::shared_ptr burnVector = u->subsetVectorForVariable(burnVar);
          setVector(Mechanics::BURNUP, burnVector);
        }
      }

      if(d_isActive[Mechanics::OXYGEN_CONCENTRATION]) {
        if(!(d_isFrozen[Mechanics::OXYGEN_CONCENTRATION])) {
          AMP::LinearAlgebra::Variable::shared_ptr oxyVar = d_inpVariables->getVariable(Mechanics::OXYGEN_CONCENTRATION); 
          AMP::LinearAlgebra::Vector::shared_ptr oxyVector = u->subsetVectorForVariable(oxyVar);
          setVector(Mechanics::OXYGEN_CONCENTRATION, oxyVector);
        }
      }

      if(d_isActive[Mechanics::LHGR]) {
        if(!(d_isFrozen[Mechanics::LHGR])) {
          AMP::LinearAlgebra::Variable::shared_ptr lhgrVar = d_inpVariables->getVariable(Mechanics::LHGR);
          AMP::LinearAlgebra::Vector::shared_ptr lhgrVector = u->subsetVectorForVariable(lhgrVar);
          setVector(Mechanics::LHGR, lhgrVector);
        }
      }

      d_materialModel->preNonlinearAssembly();

      for( ; el != end_el; ++el) {
        std::vector<unsigned int> gaussPtIndices;
        gaussPt_dof_map->getDOFs (*el, gaussPtIndices);

        unsigned int num_local_type0Dofs = 0;
        for(unsigned int i = 0; i < 3; i++) {
          (dof_maps[0])->getDOFs (*el, d_type0DofIndices[i], i);
          num_local_type0Dofs += d_type0DofIndices[i].size();
        }//end for i

        unsigned int num_local_type1Dofs = 0;
        if( d_isActive[Mechanics::TEMPERATURE] || d_isActive[Mechanics::BURNUP] ||
            d_isActive[Mechanics::OXYGEN_CONCENTRATION] || d_isActive[Mechanics::LHGR]) {
          (dof_maps[1])->getDOFs (*el, d_type1DofIndices);
          num_local_type1Dofs = d_type1DofIndices.size();
        }

        std::vector<std::vector<double> > elementInputVectors(Mechanics::TOTAL_NUMBER_OF_VARIABLES);
        std::vector<std::vector<double> > elementInputVectors_pre(Mechanics::TOTAL_NUMBER_OF_VARIABLES);

        elementInputVectors[Mechanics::DISPLACEMENT].resize(num_local_type0Dofs);
        if(d_useUpdatedLagrangian) {
          elementInputVectors_pre[Mechanics::DISPLACEMENT].resize(num_local_type0Dofs);
        }
        if(d_isActive[Mechanics::TEMPERATURE]) {
          elementInputVectors[Mechanics::TEMPERATURE].resize(num_local_type1Dofs);
          if(d_useUpdatedLagrangian) {
            elementInputVectors_pre[Mechanics::TEMPERATURE].resize(num_local_type1Dofs);
          }
        }
        if(d_isActive[Mechanics::BURNUP]) {
          elementInputVectors[Mechanics::BURNUP].resize(num_local_type1Dofs);
          if(d_useUpdatedLagrangian) {
            elementInputVectors_pre[Mechanics::BURNUP].resize(num_local_type1Dofs);
          }
        }
        if(d_isActive[Mechanics::OXYGEN_CONCENTRATION]) {
          elementInputVectors[Mechanics::OXYGEN_CONCENTRATION].resize(num_local_type1Dofs);
          if(d_useUpdatedLagrangian) {
            elementInputVectors_pre[Mechanics::OXYGEN_CONCENTRATION].resize(num_local_type1Dofs);
          }
        }
        if(d_isActive[Mechanics::LHGR]) {
          elementInputVectors[Mechanics::LHGR].resize(num_local_type1Dofs);
          if(d_useUpdatedLagrangian) {
            elementInputVectors_pre[Mechanics::LHGR].resize(num_local_type1Dofs);
          }
        }

        d_numNodesForCurrentElement = el->numNodes(); 

        for(unsigned int r = 0; r < d_numNodesForCurrentElement; r++) {
          for(unsigned int d = 0; d < 3; d++) {
            elementInputVectors[Mechanics::DISPLACEMENT][(3*r) + d] = (d_inVec[Mechanics::DISPLACEMENT])->
              getValueByGlobalID( d_type0DofIndices[d][r] );
            if(d_useUpdatedLagrangian) {
              elementInputVectors_pre[Mechanics::DISPLACEMENT][(3*r) + d] = (d_inVec_pre[Mechanics::DISPLACEMENT])->getValueByGlobalID( d_type0DofIndices[d][r] );
            }
          }
          if(d_isActive[Mechanics::TEMPERATURE]) {
            elementInputVectors[Mechanics::TEMPERATURE][r] = (d_inVec[Mechanics::TEMPERATURE])->
              getValueByGlobalID( d_type1DofIndices[r] );
            if(d_useUpdatedLagrangian) {
              elementInputVectors_pre[Mechanics::TEMPERATURE][r] = (d_inVec_pre[Mechanics::TEMPERATURE])->getValueByGlobalID( d_type1DofIndices[r] );
            }
          }
          if(d_isActive[Mechanics::BURNUP]) {
            elementInputVectors[Mechanics::BURNUP][r] = (d_inVec[Mechanics::BURNUP])->
              getValueByGlobalID( d_type1DofIndices[r] );
            if(d_useUpdatedLagrangian) {
              elementInputVectors_pre[Mechanics::BURNUP][r] = (d_inVec_pre[Mechanics::BURNUP])->getValueByGlobalID( d_type1DofIndices[r] );
            }
          }
          if(d_isActive[Mechanics::OXYGEN_CONCENTRATION]) {
            elementInputVectors[Mechanics::OXYGEN_CONCENTRATION][r] = (d_inVec[Mechanics::OXYGEN_CONCENTRATION])->
              getValueByGlobalID( d_type1DofIndices[r] );
            if(d_useUpdatedLagrangian) {
              elementInputVectors_pre[Mechanics::OXYGEN_CONCENTRATION][r] = (d_inVec_pre[Mechanics::OXYGEN_CONCENTRATION])->getValueByGlobalID( d_type1DofIndices[r] );
            }
          }
          if(d_isActive[Mechanics::LHGR]) {
            elementInputVectors[Mechanics::LHGR][r] = (d_inVec[Mechanics::LHGR])->
              getValueByGlobalID( d_type1DofIndices[r] );
            if(d_useUpdatedLagrangian) {
              elementInputVectors_pre[Mechanics::LHGR][r] = (d_inVec_pre[Mechanics::LHGR])->getValueByGlobalID( d_type1DofIndices[r] );
            }
          }
        }

        std::vector<double> elementStressVector(6*numGaussPts);
        std::vector<double> elementStrainVector(6*numGaussPts);

        const ::Elem* elemPtr = &(el->getElem());

        d_mechNonlinElem->initializeForCurrentElement( elemPtr, d_materialModel );

        d_mechNonlinElem->computeStressAndStrain(elementInputVectors, elementStressVector, elementStrainVector);

        for(unsigned int i = 0; i < (6*numGaussPts); i++) {
          stress->setValueByGlobalID(gaussPtIndices[i], elementStressVector[i]);
          strain->setValueByGlobalID(gaussPtIndices[i], elementStrainVector[i]);
        }//end for i
      }//end for el

      d_materialModel->postNonlinearAssembly();
    }

    void MechanicsNonlinearFEOperator :: updateMaterialForElementCommonFunction(const
        AMP::Mesh::MeshManager::Adapter::Element & elem, const std::vector<AMP::Mesh::DOFMap::shared_ptr> & dof_maps,
        std::vector<std::vector<double> > & elementInputVectors, std::vector<std::vector<double> > & elementInputVectors_pre )
    {
      unsigned int num_local_type0Dofs = 0;
      for(unsigned int i = 0; i < 3; i++) {
        (dof_maps[0])->getDOFs (elem, d_type0DofIndices[i], i);
        num_local_type0Dofs += d_type0DofIndices[i].size();
      }//end for i

      unsigned int num_local_type1Dofs = 0;
      if( d_isActive[Mechanics::TEMPERATURE] || d_isActive[Mechanics::BURNUP] ||
          d_isActive[Mechanics::OXYGEN_CONCENTRATION] || d_isActive[Mechanics::LHGR]) {
        (dof_maps[1])->getDOFs (elem, d_type1DofIndices);
        num_local_type1Dofs = d_type1DofIndices.size();
      }

      elementInputVectors[Mechanics::DISPLACEMENT].resize(num_local_type0Dofs);
      if(d_useUpdatedLagrangian) {
        elementInputVectors_pre[Mechanics::DISPLACEMENT].resize(num_local_type0Dofs);
      }
      if(d_isActive[Mechanics::TEMPERATURE]) {
        elementInputVectors[Mechanics::TEMPERATURE].resize(num_local_type1Dofs);
        if(d_useUpdatedLagrangian) {
          elementInputVectors_pre[Mechanics::TEMPERATURE].resize(num_local_type1Dofs);
        }
      }
      if(d_isActive[Mechanics::BURNUP]) {
        elementInputVectors[Mechanics::BURNUP].resize(num_local_type1Dofs);
        if(d_useUpdatedLagrangian) {
          elementInputVectors_pre[Mechanics::BURNUP].resize(num_local_type1Dofs);
        }
      }
      if(d_isActive[Mechanics::OXYGEN_CONCENTRATION]) {
        elementInputVectors[Mechanics::OXYGEN_CONCENTRATION].resize(num_local_type1Dofs);
        if(d_useUpdatedLagrangian) {
          elementInputVectors_pre[Mechanics::OXYGEN_CONCENTRATION].resize(num_local_type1Dofs);
        }
      }
      if(d_isActive[Mechanics::LHGR]) {
        elementInputVectors[Mechanics::LHGR].resize(num_local_type1Dofs);
        if(d_useUpdatedLagrangian) {
          elementInputVectors_pre[Mechanics::LHGR].resize(num_local_type1Dofs);
        }
      }

      d_numNodesForCurrentElement = elem.numNodes(); 

      std::vector<double> elementRefXYZ;
      elementRefXYZ.resize(3 * d_numNodesForCurrentElement);

      for(unsigned int r = 0; r < d_numNodesForCurrentElement; r++) {
        for(unsigned int d = 0; d < 3; d++) {
          elementInputVectors[Mechanics::DISPLACEMENT][(3*r) + d] = (d_inVec[Mechanics::DISPLACEMENT])->
            getValueByGlobalID( d_type0DofIndices[d][r] );
          if(d_useUpdatedLagrangian) {
            elementInputVectors_pre[Mechanics::DISPLACEMENT][(3*r) + d] = (d_inVec_pre[Mechanics::DISPLACEMENT])->getValueByGlobalID( d_type0DofIndices[d][r] );
            elementRefXYZ[(3 * r) + d] = d_refXYZ->getValueByGlobalID(d_type0DofIndices[d][r]);
          }
        }
        if(d_isActive[Mechanics::TEMPERATURE]) {
          elementInputVectors[Mechanics::TEMPERATURE][r] = (d_inVec[Mechanics::TEMPERATURE])->
            getValueByGlobalID( d_type1DofIndices[r] );
          if(d_useUpdatedLagrangian) {
            elementInputVectors_pre[Mechanics::TEMPERATURE][r] = (d_inVec_pre[Mechanics::TEMPERATURE])->getValueByGlobalID( d_type1DofIndices[r] );
          }
        }
        if(d_isActive[Mechanics::BURNUP]) {
          elementInputVectors[Mechanics::BURNUP][r] = (d_inVec[Mechanics::BURNUP])->
            getValueByGlobalID( d_type1DofIndices[r] );
          if(d_useUpdatedLagrangian) {
            elementInputVectors_pre[Mechanics::BURNUP][r] = (d_inVec_pre[Mechanics::BURNUP])->getValueByGlobalID( d_type1DofIndices[r] );
          }
        }
        if(d_isActive[Mechanics::OXYGEN_CONCENTRATION]) {
          elementInputVectors[Mechanics::OXYGEN_CONCENTRATION][r] = (d_inVec[Mechanics::OXYGEN_CONCENTRATION])->
            getValueByGlobalID( d_type1DofIndices[r] );
          if(d_useUpdatedLagrangian) {
            elementInputVectors_pre[Mechanics::OXYGEN_CONCENTRATION][r] = (d_inVec_pre[Mechanics::OXYGEN_CONCENTRATION])->getValueByGlobalID( d_type1DofIndices[r] );
          }
        }
        if(d_isActive[Mechanics::LHGR]) {
          elementInputVectors[Mechanics::LHGR][r] = (d_inVec[Mechanics::LHGR])->
            getValueByGlobalID( d_type1DofIndices[r] );
          if(d_useUpdatedLagrangian) {
            elementInputVectors_pre[Mechanics::LHGR][r] = (d_inVec_pre[Mechanics::LHGR])->getValueByGlobalID( d_type1DofIndices[r] );
          }
        }
      }

      const ::Elem* elemPtr = &(elem.getElem());

      if(d_useUpdatedLagrangian) {
        d_mechNULElem->initializeForCurrentElement( elemPtr, d_materialModel );
        d_mechNULElem->assignReferenceXYZ(elementRefXYZ);
      } else {
        d_mechNonlinElem->initializeForCurrentElement( elemPtr, d_materialModel );
      }
    }

  }
}//end namespace


