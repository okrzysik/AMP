
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
        d_outVariables.reset(new AMP::LinearAlgebra::MultiVariable("myOutVar"));
        for(unsigned int i = 0; i < NavierStokes::TOTAL_NUMBER_OF_VARIABLES; i++) {
          AMP::LinearAlgebra::Variable::shared_ptr dummyVar;
          d_inpVariables->add(dummyVar);
        }//end for i

        std::string tempVarName = activeInpVar_db->getString("VELOCITY");
        AMP::LinearAlgebra::Variable::shared_ptr tempVar(new AMP::LinearAlgebra::Variable(tempVarName) ); 
        d_inpVariables->setVariable(NavierStokes::VELOCITY, tempVar);

        tempVarName = activeInpVar_db->getString("PRESSURE");
        AMP::LinearAlgebra::Variable::shared_ptr tempVar2(new AMP::LinearAlgebra::Variable(tempVarName) ); 
        d_inpVariables->setVariable(NavierStokes::PRESSURE, tempVar2);

        if(d_isActive[NavierStokes::TEMPERATURE]) {
          std::string varName = activeInpVar_db->getString("TEMPERATURE");
          AMP::LinearAlgebra::Variable::shared_ptr dummyVar(new AMP::LinearAlgebra::Variable(varName) );
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
                                        inpVar.reset(new AMP::LinearAlgebra::Variable(name) );
                                        break;
                                      }
        case NavierStokes::PRESSURE : {
                                        inpVar.reset(new AMP::LinearAlgebra::Variable(name) );
                                        break;
                                      }
        case NavierStokes::TEMPERATURE : {
                                           inpVar.reset(new AMP::LinearAlgebra::Variable(name) );
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

      AMP::LinearAlgebra::Variable::shared_ptr tempVar; 
      AMP::LinearAlgebra::Vector::shared_ptr tempVector; 

      tempVar = d_inpVariables->getVariable(NavierStokes::VELOCITY);
      tempVector = u->subsetVectorForVariable(tempVar);
      setVector(NavierStokes::VELOCITY, tempVector);

      tempVar = d_inpVariables->getVariable(NavierStokes::PRESSURE);
      tempVector = u->subsetVectorForVariable(tempVar);
      setVector(NavierStokes::PRESSURE, tempVector);

      if(d_isActive[NavierStokes::TEMPERATURE]) {
        if(!(d_isFrozen[NavierStokes::TEMPERATURE])) {
          tempVar = d_inpVariables->getVariable(NavierStokes::TEMPERATURE); 
          tempVector = u->subsetVectorForVariable(tempVar);
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

    // preelement operation for this operator assumes that the pressure and
    // temperature use same discretization and order of approximation 
    void NavierStokesGalWFFEOperator :: preElementOperation( const AMP::Mesh::MeshElement & elem) 
    {
      d_currNodes = elem.getElements(AMP::Mesh::Vertex);
      unsigned int numNodesInCurrElem = d_currNodes.size();

      gettype0DofIndicesForCurrentElement(NavierStokes::VELOCITY, d_type0DofIndices);
      gettype1DofIndicesForCurrentElement(NavierStokes::PRESSURE, d_type1DofIndices);

      std::vector<std::vector<double> > elementInputVectors(NavierStokes::TOTAL_NUMBER_OF_VARIABLES);

      elementInputVectors[NavierStokes::VELOCITY].resize(3*numNodesInCurrElem);
      elementInputVectors[NavierStokes::PRESSURE].resize(numNodesInCurrElem);

      if(d_isActive[NavierStokes::TEMPERATURE]) {
        elementInputVectors[NavierStokes::TEMPERATURE].resize(numNodesInCurrElem);
      }

      for(unsigned int r = 0; r < numNodesInCurrElem ; r++) {
        for(unsigned int d = 0; d < 3; d++) {
          elementInputVectors[NavierStokes::VELOCITY][(3*r) + d] = (d_inVec[NavierStokes::VELOCITY])->
            getValueByGlobalID( d_type0DofIndices[r][d] );
        }

        elementInputVectors[NavierStokes::PRESSURE][r] = (d_inVec[NavierStokes::PRESSURE])->
          getValueByGlobalID( d_type1DofIndices[r] );

        if(d_isActive[NavierStokes::TEMPERATURE]) {
          elementInputVectors[NavierStokes::TEMPERATURE][r] = (d_inVec[NavierStokes::TEMPERATURE])->
            getValueByGlobalID( d_type1DofIndices[r] );
        }
      }

      d_elementOutputVector.resize(4*numNodesInCurrElem );
      for(unsigned int i = 0; i < 4*numNodesInCurrElem ; i++) {
        d_elementOutputVector[i] = 0.0;
      }

      d_flowGalWFElem->initializeForCurrentElement( d_currElemPtrs[d_currElemIdx], d_transportModel );
      d_flowGalWFElem->setElementVectors( elementInputVectors, d_elementOutputVector );
    }

    void NavierStokesGalWFFEOperator :: postElementOperation()
    {
      for(unsigned int r = 0; r < d_type0DofIndices.size() ; r++) {
        AMP_ASSERT(d_type0DofIndices[r].size() == 3);
        for(unsigned int d = 0; d < 3; d++) {
          d_outVec->addValueByGlobalID( d_type0DofIndices[r][d], d_elementOutputVector[(4*r) + d] );
        }
      }
      for(unsigned int r = 0; r < d_type1DofIndices.size() ; r++) {
        d_outVec->addValueByGlobalID( d_type1DofIndices[r], d_elementOutputVector[(4*r) + 3] );
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

    void NavierStokesGalWFFEOperator :: gettype0DofIndicesForCurrentElement(int varId, std::vector<std::vector<size_t> > & dofIds) {
      dofIds.resize(d_currNodes.size());
      for(unsigned int j = 0; j < d_currNodes.size(); j++) {
        d_dofMap[varId]->getDOFs(d_currNodes[j].globalID(), dofIds[j]);
      } // end of j
    }

    void NavierStokesGalWFFEOperator :: gettype1DofIndicesForCurrentElement(int varId, std::vector<size_t> & dofIds) {
      for(unsigned int j = 0; j < d_currNodes.size(); j++) {
        d_dofMap[varId]->getDOFs(d_currNodes[j].globalID(), dofIds);
      } // end of j
    }
  }
}//end namespace
