
#include "ConsMomentumGalWFFEOperator.h"
#include "ConsMomentumGalWFLinearFEOperatorParameters.h"
#include "utils/Utilities.h"
#include "utils/InputDatabase.h"

namespace AMP {
namespace Operator {

  ConsMomentumGalWFFEOperator :: ConsMomentumGalWFFEOperator (
      const boost::shared_ptr<ConsMomentumGalWFFEOperatorParameters> & params)
	  : NonlinearFEOperator (params) {

		  AMP_INSIST( ((params.get()) != NULL), "NULL parameter!" );
		  AMP_INSIST( (((params->d_db).get()) != NULL), "NULL database!" );

		  d_transportModel = params->d_transportModel;

		  d_isActive.resize(NavierStokes::TOTAL_NUMBER_OF_VARIABLES);
		  d_isFrozen.resize(NavierStokes::TOTAL_NUMBER_OF_VARIABLES);
		  d_inVec.resize(NavierStokes::TOTAL_NUMBER_OF_VARIABLES);

		  AMP_INSIST( params->d_db->keyExists("ActiveInputVariables"), "key not found" );
		  boost::shared_ptr<AMP::Database> activeInpVar_db = params->d_db->getDatabase("ActiveInputVariables");

		  AMP_INSIST(activeInpVar_db->keyExists("VELOCITY"), "VELOCITY must be active");
		  AMP_INSIST(activeInpVar_db->keyExists("PRESSURE"), "PRESSURE must be active");

		  d_isActive[NavierStokes::VELOCITY] = true;
		  d_isActive[NavierStokes::PRESSURE] = activeInpVar_db->keyExists("PRESSURE");
		  d_isActive[NavierStokes::TEMPERATURE] = activeInpVar_db->keyExists("TEMPERATURE");

		  d_isFrozen[NavierStokes::VELOCITY] = false;
		  d_isFrozen[NavierStokes::PRESSURE] =  (params->d_db)->getBoolWithDefault("FREEZE_PRESSURE", true);
		  d_isFrozen[NavierStokes::TEMPERATURE] = (params->d_db)->getBoolWithDefault("FREEZE_TEMPERATURE", true);

		  d_inpVariables.reset(new AMP::LinearAlgebra::MultiVariable("myInpVar"));

		  for(unsigned int i = 0; i < NavierStokes::TOTAL_NUMBER_OF_VARIABLES; i++) {
			  AMP::LinearAlgebra::Variable::shared_ptr dummyVar;
			  d_inpVariables->add(dummyVar);
		  }//end for i

		  std::string tempVarName = activeInpVar_db->getString("VELOCITY");
		  AMP::LinearAlgebra::Variable::shared_ptr tempVar(new AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable, 3>(tempVarName, d_MeshAdapter) ); 
		  d_inpVariables->setVariable(NavierStokes::VELOCITY, tempVar);

		  if(d_isActive[NavierStokes::PRESSURE]) {
			  std::string varName = activeInpVar_db->getString("PRESSURE");
			  AMP::LinearAlgebra::Variable::shared_ptr dummyVar(new AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable, 1>(varName, d_MeshAdapter) );
			  d_inpVariables->setVariable(NavierStokes::PRESSURE, dummyVar);
			  if(d_isFrozen[NavierStokes::PRESSURE]) {
				  if( params->d_FrozenPressure != NULL ) {
					  setVector(NavierStokes::PRESSURE, params->d_FrozenPressure);
				  }
			  }
		  }

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

		  AMP_INSIST( params->d_db->keyExists("OutputVariable"), "key not found" );
		  std::string outVar = params->d_db->getString("OutputVariable");
		  d_outVariables.reset(new AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable, 3>(outVar, d_MeshAdapter));
		  d_isInitialized = false;
	  }

  unsigned int ConsMomentumGalWFFEOperator :: numberOfDOFMaps() {
	  if( d_isActive[NavierStokes::TEMPERATURE] || d_isActive[NavierStokes::PRESSURE]) {
		  return 2;
	  } else {
		  return 1; 
	  }
  }

  AMP::LinearAlgebra::Variable::shared_ptr ConsMomentumGalWFFEOperator :: getVariableForDOFMap(unsigned int id) {
	  AMP_ASSERT( id < (this->numberOfDOFMaps()) );
	  if(id == 0) {
		  return (d_inpVariables->getVariable(NavierStokes::VELOCITY));
	  } else {
		  if(d_isActive[NavierStokes::PRESSURE]) {
			  return (d_inpVariables->getVariable(NavierStokes::PRESSURE));
		  //} else if(d_isActive[NavierStokes::TEMPERATURE]) {
		  } else {
			  return (d_inpVariables->getVariable(NavierStokes::TEMPERATURE));
		  }    
	  }
  }

  AMP::LinearAlgebra::Variable::shared_ptr ConsMomentumGalWFFEOperator :: createInputVariable(const std::string & name, int varId) {
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

  void ConsMomentumGalWFFEOperator :: preAssembly(const boost::shared_ptr<AMP::LinearAlgebra::Vector>  &u, 
		  boost::shared_ptr<AMP::LinearAlgebra::Vector>  &r) {
	  if(!d_isInitialized) {
		  init();
	  }

	  AMP::LinearAlgebra::Variable::shared_ptr tempVar = d_inpVariables->getVariable(NavierStokes::VELOCITY);
	  AMP::LinearAlgebra::Vector::shared_ptr tempVector = u->subsetVectorForVariable(tempVar);
	  setVector(NavierStokes::VELOCITY, tempVector);

	  if(d_isActive[NavierStokes::PRESSURE]) {
		  if(!(d_isFrozen[NavierStokes::PRESSURE])) {
			  AMP::LinearAlgebra::Variable::shared_ptr tempVar = d_inpVariables->getVariable(NavierStokes::PRESSURE); 
			  AMP::LinearAlgebra::Vector::shared_ptr tempVector = u->subsetVectorForVariable(tempVar);
			  setVector(NavierStokes::PRESSURE, tempVector);
		  }
	  }

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

  void ConsMomentumGalWFFEOperator :: postAssembly()
  {
	  d_outVec->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_ADD );
  }

  void ConsMomentumGalWFFEOperator :: preElementOperation( const AMP::Mesh::MeshManager::Adapter::Element & elem, 
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

	  if(d_isActive[NavierStokes::TEMPERATURE]) {
		  elementInputVectors[NavierStokes::PRESSURE].resize(num_local_type1Dofs);
          }
	  if(d_isActive[NavierStokes::TEMPERATURE]) {
		  elementInputVectors[NavierStokes::TEMPERATURE].resize(num_local_type1Dofs);
	  }

	  d_numNodesForCurrentElement = elem.numNodes(); 

	  for(unsigned int r = 0; r < d_numNodesForCurrentElement; r++) {
		  for(unsigned int d = 0; d < 3; d++) {
			  elementInputVectors[NavierStokes::VELOCITY][(3*r) + d] = (d_inVec[NavierStokes::VELOCITY])->
				  getValueByGlobalID( d_type0DofIndices[d][r] );
		  }

		  if(d_isActive[NavierStokes::PRESSURE]) {
			  elementInputVectors[NavierStokes::PRESSURE][r] = (d_inVec[NavierStokes::PRESSURE])->
				  getValueByGlobalID( d_type1DofIndices[r] );
		  }

		  if(d_isActive[NavierStokes::TEMPERATURE]) {
			  elementInputVectors[NavierStokes::TEMPERATURE][r] = (d_inVec[NavierStokes::TEMPERATURE])->
				  getValueByGlobalID( d_type1DofIndices[r] );
		  }
	  }

	  d_elementOutputVector.resize(num_local_type0Dofs );
	  for(unsigned int i = 0; i < num_local_type0Dofs ; i++) {
		  d_elementOutputVector[i] = 0.0;
	  }

	  const ::Elem* elemPtr = &(elem.getElem());

	  d_flowGalWFElem->initializeForCurrentElement( elemPtr, d_transportModel );
	  d_flowGalWFElem->setElementVectors( elementInputVectors, d_elementOutputVector );
  }

  void ConsMomentumGalWFFEOperator :: postElementOperation()
  {
	  for(unsigned int r = 0; r < d_numNodesForCurrentElement; r++) {
		  for(unsigned int d = 0; d < 3; d++) {
			  d_outVec->addValueByGlobalID( d_type0DofIndices[d][r], d_elementOutputVector[(3*r) + d] );
		  }
	  }
  }

  void ConsMomentumGalWFFEOperator :: init() 
  {
    d_isInitialized = true;
  }

  void ConsMomentumGalWFFEOperator :: reset(const boost::shared_ptr<OperatorParameters>& params)
  {
	  if(!d_isInitialized) {
		  init();
	  }

	  boost::shared_ptr<ConsMomentumGalWFFEOperatorParameters> myParams =
		  boost::dynamic_pointer_cast<ConsMomentumGalWFFEOperatorParameters>(params); 

	  AMP_INSIST( ((myParams.get()) != NULL), "Null parameter!" );

	  unsigned int numDOFMaps = numberOfDOFMaps();
	  std::vector<AMP::Mesh::DOFMap::shared_ptr> dof_maps(numDOFMaps);

	  for(unsigned int i = 0; i < numDOFMaps; i++) {
		  dof_maps[i] = d_MeshAdapter->getDOFMap( getVariableForDOFMap(i) );
	  }

	  d_outVec.reset();

	  if(d_isActive[NavierStokes::PRESSURE]) {
		  if(d_isFrozen[NavierStokes::PRESSURE]) {
			  if( myParams->d_FrozenPressure != NULL ) {
				  setVector(NavierStokes::PRESSURE, myParams->d_FrozenPressure);
			  }
		  }
	  }


	  if(d_isActive[NavierStokes::TEMPERATURE]) {
		  if(d_isFrozen[NavierStokes::TEMPERATURE]) {
			  if( myParams->d_FrozenTemperature != NULL ) {
				  setVector(NavierStokes::TEMPERATURE, myParams->d_FrozenTemperature);
			  }
		  }
	  }

  }

  boost::shared_ptr<OperatorParameters> ConsMomentumGalWFFEOperator ::
	  getJacobianParameters(const boost::shared_ptr<AMP::LinearAlgebra::Vector>& u) {

		  if(!d_isInitialized) {
			  init();
		  }

		  // set up a database for the linear operator params
		  boost::shared_ptr<AMP::InputDatabase> tmp_db (new AMP::InputDatabase("Dummy"));
		  tmp_db->putBool("isAttachedToNonlinearOperator", true);
		  tmp_db->putBool("isNonlinearOperatorInitialized", true);

		  // create the linear operator params
		  boost::shared_ptr<ConsMomentumGalWFLinearFEOperatorParameters> outParams(new
				  ConsMomentumGalWFLinearFEOperatorParameters(tmp_db));

		  (outParams->d_frozenVec).resize(NavierStokes::TOTAL_NUMBER_OF_VARIABLES);
		  // add variables to parameters
		  if (d_isActive[NavierStokes::PRESSURE]) {
			  if (d_isFrozen[NavierStokes::PRESSURE]) {
				  outParams->d_frozenVec[NavierStokes::PRESSURE] = d_inVec[NavierStokes::PRESSURE];
				  outParams->d_frozenVec[NavierStokes::PRESSURE]->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
			  } else {
				  AMP::LinearAlgebra::Variable::shared_ptr tvar = d_inpVariables->getVariable(NavierStokes::PRESSURE);
				  AMP::LinearAlgebra::Vector::shared_ptr pressure = u->subsetVectorForVariable(tvar);
				  outParams->d_frozenVec[NavierStokes::PRESSURE] = pressure ;
				  outParams->d_frozenVec[NavierStokes::PRESSURE]->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
			  }
		  }
		  if (d_isActive[NavierStokes::TEMPERATURE]) {
			  if (d_isFrozen[NavierStokes::TEMPERATURE]) {
				  outParams->d_frozenVec[NavierStokes::TEMPERATURE] = d_inVec[NavierStokes::TEMPERATURE];
				  outParams->d_frozenVec[NavierStokes::TEMPERATURE]->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
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
