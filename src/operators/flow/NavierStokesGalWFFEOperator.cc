
#include "NavierStokesGalWFFEOperator.h"
#include "NavierStokesGalWFLinearFEOperatorParameters.h"
#include "utils/Utilities.h"
#include "utils/InputDatabase.h"

namespace AMP {
namespace Operator {

  NavierStokesGalWFFEOperator :: NavierStokesGalWFFEOperator (
      const boost::shared_ptr<NavierStokesGalWFFEOperatorParameters> & params)
    : NonlinearFEOperator (params) {
      
      AMP_INSIST( ((params.get()) != NULL), "NULL parameter!" );
      AMP_INSIST( (((params->d_db).get()) != NULL), "NULL database!" );

      d_transportModel = params->d_transportModel;

      d_isActive.resize(NavierStokes::TOTAL_NUMBER_OF_VARIABLES);
      d_isFrozen.resize(NavierStokes::TOTAL_NUMBER_OF_VARIABLES);
      d_inVec.resize(NavierStokes::TOTAL_NUMBER_OF_VARIABLES);

      d_coupledFormulation = (params->d_db)->getBoolWithDefault("VELOCITY_PRESSURE_COUPLING", true);

      AMP_INSIST( params->d_db->keyExists("ActiveInputVariables"), "key not found" );
      boost::shared_ptr<AMP::Database> activeInpVar_db = params->d_db->getDatabase("ActiveInputVariables");

      AMP_INSIST(activeInpVar_db->keyExists("VELOCITY"), "VELOCITY must be active");
      AMP_INSIST(activeInpVar_db->keyExists("PRESSURE"), "PRESSURE must be active");

      d_isActive[NavierStokes::VELOCITY] = true;
      d_isActive[NavierStokes::PRESSURE] = true;
      d_isActive[NavierStokes::TEMPERATURE] = activeInpVar_db->keyExists("TEMPERATURE");

      d_isFrozen[NavierStokes::VELOCITY] = false;
      d_isFrozen[NavierStokes::PRESSURE] = false;
      d_isFrozen[NavierStokes::TEMPERATURE] = (params->d_db)->getBoolWithDefault("FREEZE_TEMPERATURE", true);

      d_inpVariables.reset(new AMP::LinearAlgebra::MultiVariable("myInpVar"));
      d_outVariables.reset(new AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable, 4>("myInpVar", d_MeshAdapter));
      for(unsigned int i = 0; i < NavierStokes::TOTAL_NUMBER_OF_VARIABLES; i++) {
        AMP::LinearAlgebra::Variable::shared_ptr dummyVar;
        d_inpVariables->add(dummyVar);
      }//end for i

      std::string tempVarName = activeInpVar_db->getString("VELOCITY");
      AMP::LinearAlgebra::Variable::shared_ptr tempVar(new AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable, 3>(tempVarName, d_MeshAdapter) ); 
      d_inpVariables->setVariable(NavierStokes::VELOCITY, tempVar);

      tempVarName = activeInpVar_db->getString("PRESSURE");
      AMP::LinearAlgebra::Variable::shared_ptr tempVar2(new AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable, 1>(tempVarName, d_MeshAdapter) ); 
      d_inpVariables->setVariable(NavierStokes::PRESSURE, tempVar2);

      if(d_isActive[NavierStokes::TEMPERATURE]) {
        std::string varName = activeInpVar_db->getString("TEMPERATURE");
        AMP::LinearAlgebra::Variable::shared_ptr dummyVar(new AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable, 1>(varName, d_MeshAdapter) );
        d_inpVariables->setVariable(NavierStokes::TEMPERATURE, dummyVar);
        if(d_isFrozen[NavierStokes::TEMPERATURE]) {
          if( params->d_FrozenTemperature != NULL ) {
            setVector(NavierStokes::TEMPERATURE, params->d_FrozenTemperature);
          }
        }
      }

      d_isInitialized = false;
    }

  unsigned int NavierStokesGalWFFEOperator :: numberOfDOFMaps() {
      return 1;
  }

  AMP::LinearAlgebra::Variable::shared_ptr NavierStokesGalWFFEOperator :: getVariableForDOFMap(unsigned int id) {
    AMP_ASSERT( id < (this->numberOfDOFMaps()) );
    if(id == 0) {
      return (d_inpVariables->getVariable(NavierStokes::VELOCITY));
    } else if(id == 1) {
      return (d_inpVariables->getVariable(NavierStokes::PRESSURE));
    } else {
//      if(d_isActive[NavierStokes::TEMPERATURE]) {
      AMP_INSIST(d_isActive[NavierStokes::TEMPERATURE],"if it fails this there's a problem"); 
      return (d_inpVariables->getVariable(NavierStokes::TEMPERATURE));
//      }    
    }
  }

  AMP::LinearAlgebra::Variable::shared_ptr NavierStokesGalWFFEOperator :: createInputVariable(const std::string & name, int varId) {
    AMP::LinearAlgebra::Variable::shared_ptr inpVar;
    switch(varId) {
      case NavierStokes::VELOCITY : {
                                       inpVar.reset(new AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable, 3>(name) );
                                       break;
                                    }
      case NavierStokes::PRESSURE : {
                                       inpVar.reset(new AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable, 1>(name) );
                                       break;
                                    }
      case NavierStokes::TEMPERATURE : {
                                      inpVar.reset(new AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable, 1>(name) );
                                      break;
                                    }
      default: 
                             assert(false);
    }
    return inpVar;
  }

  void NavierStokesGalWFFEOperator :: preAssembly(const boost::shared_ptr<AMP::LinearAlgebra::Vector>  &u, 
      boost::shared_ptr<AMP::LinearAlgebra::Vector>  &r) {
    if(!d_isInitialized) {
      init();
    }

    AMP::LinearAlgebra::Variable::shared_ptr tempVar = d_inpVariables->getVariable(NavierStokes::VELOCITY);
    AMP::LinearAlgebra::Vector::shared_ptr tempVector = u->subsetVectorForVariable(tempVar);
    setVector(NavierStokes::VELOCITY, tempVector);

    tempVar = d_inpVariables->getVariable(NavierStokes::PRESSURE);
    tempVector = u->subsetVectorForVariable(tempVar);
    setVector(NavierStokes::PRESSURE, tempVector);

    if(d_isActive[NavierStokes::TEMPERATURE]) {
      if(!(d_isFrozen[NavierStokes::TEMPERATURE])) {
        AMP::LinearAlgebra::Variable::shared_ptr tempVar = d_inpVariables->getVariable(NavierStokes::TEMPERATURE); 
        AMP::LinearAlgebra::Vector::shared_ptr tempVector = u->subsetVectorForVariable(tempVar);
        setVector(NavierStokes::TEMPERATURE, tempVector);
      }
    }

    d_outVec = r->subsetVectorForVariable(d_outVariables);
    d_outVec->zero();

  }

  void NavierStokesGalWFFEOperator :: postAssembly()
  {
    d_outVec->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_ADD );
  }

  void NavierStokesGalWFFEOperator :: preElementOperation( const AMP::Mesh::MeshManager::Adapter::Element & elem, 
      const std::vector<AMP::Mesh::DOFMap::shared_ptr> & dof_maps )
  {
    unsigned int num_local_type0Dofs = 0;
    for(unsigned int i = 0; i < 3; i++) {
      (dof_maps[0])->getDOFs (elem, d_type0DofIndices[i], i);
      num_local_type0Dofs += d_type0DofIndices[i].size();
    }//end for i

    unsigned int num_local_type1Dofs = 0;
    (dof_maps[1])->getDOFs (elem, d_type1DofIndices);
    num_local_type1Dofs = d_type1DofIndices.size();

    std::vector<std::vector<double> > elementInputVectors(NavierStokes::TOTAL_NUMBER_OF_VARIABLES);

    elementInputVectors[NavierStokes::VELOCITY].resize(num_local_type0Dofs);
    elementInputVectors[NavierStokes::PRESSURE].resize(num_local_type1Dofs);

    if(d_isActive[NavierStokes::TEMPERATURE]) {
      elementInputVectors[NavierStokes::TEMPERATURE].resize(num_local_type1Dofs);
    }

    d_numNodesForCurrentElement = elem.numNodes(); 

    for(unsigned int r = 0; r < d_numNodesForCurrentElement; r++) {
      for(unsigned int d = 0; d < 3; d++) {
        elementInputVectors[NavierStokes::VELOCITY][(3*r) + d] = (d_inVec[NavierStokes::VELOCITY])->
          getValueByGlobalID( d_type0DofIndices[d][r] );
      }
      
      elementInputVectors[NavierStokes::PRESSURE][r] = (d_inVec[NavierStokes::PRESSURE])->
        getValueByGlobalID( d_type1DofIndices[r] );
      
      if(d_isActive[NavierStokes::TEMPERATURE]) {
        elementInputVectors[NavierStokes::TEMPERATURE][r] = (d_inVec[NavierStokes::TEMPERATURE])->
          getValueByGlobalID( d_type1DofIndices[r] );
      }
    }

    d_elementOutputVector.resize(num_local_type0Dofs + num_local_type1Dofs);
    for(unsigned int i = 0; i < num_local_type0Dofs + num_local_type1Dofs; i++) {
      d_elementOutputVector[i] = 0.0;
    }

    const ::Elem* elemPtr = &(elem.getElem());

      d_flowGalWFElem->initializeForCurrentElement( elemPtr, d_transportModel );
      d_flowGalWFElem->setElementVectors( elementInputVectors, d_elementOutputVector );
  }

  void NavierStokesGalWFFEOperator :: postElementOperation()
  {
    for(unsigned int r = 0; r < d_numNodesForCurrentElement; r++) {
      for(unsigned int d = 0; d < 3; d++) {
        d_outVec->addValueByGlobalID( d_type0DofIndices[d][r], d_elementOutputVector[(3*r) + d] );
      }
    }
  }

  void NavierStokesGalWFFEOperator :: init() 
  {
    d_isInitialized = true;
  }

  void NavierStokesGalWFFEOperator :: reset(const boost::shared_ptr<OperatorParameters>& params)
  {
    if(!d_isInitialized) {
      init();
    }

    boost::shared_ptr<NavierStokesGalWFFEOperatorParameters> myParams =
      boost::dynamic_pointer_cast<NavierStokesGalWFFEOperatorParameters>(params); 

    AMP_INSIST( ((myParams.get()) != NULL), "Null parameter!" );

      unsigned int numDOFMaps = numberOfDOFMaps();
      std::vector<AMP::Mesh::DOFMap::shared_ptr> dof_maps(numDOFMaps);

      for(unsigned int i = 0; i < numDOFMaps; i++) {
        dof_maps[i] = d_MeshAdapter->getDOFMap( getVariableForDOFMap(i) );
      }

      d_outVec.reset();

    if(d_isActive[NavierStokes::TEMPERATURE]) {
      if(d_isFrozen[NavierStokes::TEMPERATURE]) {
        if( myParams->d_FrozenTemperature != NULL ) {
          setVector(NavierStokes::TEMPERATURE, myParams->d_FrozenTemperature);
        }
      }
    }

  }

  boost::shared_ptr<OperatorParameters> NavierStokesGalWFFEOperator ::
    getJacobianParameters(const boost::shared_ptr<AMP::LinearAlgebra::Vector>& u) {

      if(!d_isInitialized) {
        init();
      }

      // set up a database for the linear operator params
      boost::shared_ptr<AMP::InputDatabase> tmp_db (new AMP::InputDatabase("Dummy"));
      tmp_db->putBool("isAttachedToNonlinearOperator", true);
      tmp_db->putBool("isNonlinearOperatorInitialized", true);

      if(d_coupledFormulation)
      {
        tmp_db->putString("name", "NavierStokesGalWFLinearFEOperator");
        tmp_db->putString("InputVariable", NavierStokes::VELOCITY);
        tmp_db->putString("OutputVariable", d_outVariables->getName());
        tmp_db->putBool("FixedTemperature", d_isActive[NavierStokes::TEMPERATURE]?false:true);
      }

      // create the linear operator params
      boost::shared_ptr<NavierStokesGalWFLinearFEOperatorParameters> outParams(new
          NavierStokesGalWFLinearFEOperatorParameters(tmp_db));

      (outParams->d_frozenVec).resize(NavierStokes::TOTAL_NUMBER_OF_VARIABLES);
      // add variables to parameters
      if (d_isActive[NavierStokes::TEMPERATURE]) {
        if (d_isFrozen[NavierStokes::TEMPERATURE]) {
            outParams->d_frozenVec[NavierStokes::TEMPERATURE] = d_inVec[NavierStokes::TEMPERATURE];
        } else {
            AMP::LinearAlgebra::Variable::shared_ptr tvar = d_inpVariables->getVariable(NavierStokes::TEMPERATURE);
            AMP::LinearAlgebra::Vector::shared_ptr temperature = u->subsetVectorForVariable(tvar);
            outParams->d_frozenVec[NavierStokes::TEMPERATURE] = temperature;
            outParams->d_frozenVec[NavierStokes::TEMPERATURE]->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
        }
      }

      return outParams;
    
    }

}
}//end namespace
