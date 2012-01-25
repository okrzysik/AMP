#include "DiffusionNonlinearFEOperator.h"
#include "DiffusionLinearFEOperatorParameters.h"
#include "ElementOperationParameters.h"
#include "DiffusionLinearElement.h"
#include "DiffusionConstants.h"
#include "utils/Utilities.h"
#include "utils/InputDatabase.h"
#include "materials/Material.h"

#include <cstring>
#include <iostream>

namespace AMP {
namespace Operator {


// Functions moved from header

    AMP::LinearAlgebra::Variable::shared_ptr DiffusionNonlinearFEOperator::createInputVariable(const std::string & name, int varId) {
        if (varId == -1) {
            return d_inpVariables->cloneVariable(name);
        } else {
            return (d_inpVariables->getVariable(varId))->cloneVariable(name);
        }
    }

    AMP::LinearAlgebra::Variable::shared_ptr DiffusionNonlinearFEOperator::createOutputVariable(const std::string & name, int varId) {
        (void) varId;
        return d_outVariable->cloneVariable(name);
    }

    AMP::LinearAlgebra::Variable::shared_ptr DiffusionNonlinearFEOperator::getInputVariable() {
        return d_inpVariables;
    }

    AMP::LinearAlgebra::Variable::shared_ptr DiffusionNonlinearFEOperator::getOutputVariable() {
            return d_outVariable;
    }

    unsigned int DiffusionNonlinearFEOperator::numberOfDOFMaps() {
            return 1;
    }

    AMP::LinearAlgebra::Variable::shared_ptr DiffusionNonlinearFEOperator::getVariableForDOFMap(unsigned int id) {
      (void) id;
            return d_inpVariables->getVariable(d_PrincipalVariable);
    }

    unsigned int DiffusionNonlinearFEOperator::getPrincipalVariableId(){return d_PrincipalVariable;}

    std::vector<unsigned int> DiffusionNonlinearFEOperator::getNonPrincipalVariableIds(){
        std::vector<unsigned int> ids;
        for (size_t i=0; i<Diffusion::NUMBER_VARIABLES; i++) {
            if (i != d_PrincipalVariable and d_isActive[i]) ids.push_back(i);
        }
        return ids;
    }

    boost::shared_ptr<DiffusionTransportModel> DiffusionNonlinearFEOperator::getTransportModel(){
        return d_transportModel;
    }

    std::vector<AMP::LinearAlgebra::Vector::shared_ptr> DiffusionNonlinearFEOperator::getFrozen(){return d_Frozen;}

    void DiffusionNonlinearFEOperator::setVector(unsigned int id, AMP::LinearAlgebra::Vector::shared_ptr &frozenVec)
    {
      d_Frozen[id] = frozenVec->subsetVectorForVariable(d_inpVariables->getVariable(id));
      (d_Frozen[id])->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
    }




DiffusionNonlinearFEOperator::DiffusionNonlinearFEOperator(
        const boost::shared_ptr<DiffusionNonlinearFEOperatorParameters> &params) :
    NonlinearFEOperator(params), d_Frozen(Diffusion::NUMBER_VARIABLES)
{
    AMP_INSIST( ((params.get()) != NULL), "NULL parameter!" );

    d_diffNonlinElem = boost::dynamic_pointer_cast<DiffusionNonlinearElement>(
            d_elemOp);

    AMP_INSIST( ((d_diffNonlinElem.get()) != NULL), "d_elemOp is not of type DiffusionNonlinearElement" );

    d_transportModel = params->d_transportModel;

    d_isActive.resize(Diffusion::NUMBER_VARIABLES);
    d_isFrozen.resize(Diffusion::NUMBER_VARIABLES);
    d_inVec.resize(Diffusion::NUMBER_VARIABLES);

    boost::shared_ptr<AMP::Database> activeVariables_db =
            params->d_db->getDatabase("ActiveInputVariables");

    for (size_t var=0; var<Diffusion::NUMBER_VARIABLES; var++) {
        std::string namespec = activeVariables_db->getStringWithDefault(
                Diffusion::names[var], "not_specified");
        bool isthere = namespec != "not_specified";
        d_isActive[var] = isthere;
        std::string frozenName = "Freeze"+Diffusion::names[var];
        d_isFrozen[var] = params->d_db->getBoolWithDefault(frozenName, false);
        if (d_isFrozen[var]) AMP_INSIST(d_isActive[var], "can not freeze a variable unless it is active");
    }

    d_numberActive = 0; d_numberFrozen = 0;
    for (unsigned int i = 0; i < Diffusion::NUMBER_VARIABLES; i++)
    {
        if (d_isActive[i]) d_numberActive++;
        if (d_isFrozen[i]) d_numberFrozen++;
    }

    AMP_INSIST(params->d_db->keyExists("PrincipalVariable"), "must specify PrincipalVariable");
    int value = params->d_db->getInteger("PrincipalVariable");
    AMP_INSIST(value>=0, "can not specify negative PrincipalVariable");
    AMP_INSIST(static_cast<unsigned int>(value)<Diffusion::NUMBER_VARIABLES, "PrincipalVariable is too large");
    d_PrincipalVariable = static_cast<unsigned int>(value);

    AMP_INSIST(d_isActive[d_PrincipalVariable], "must have Principal_Variable active");
    d_diffNonlinElem->setPrincipalVariable(d_PrincipalVariable);

    resetFrozen(params);

    d_inpVariables.reset(new AMP::LinearAlgebra::MultiVariable("InputVariables"));
    for(unsigned int i = 0; i < Diffusion::NUMBER_VARIABLES; i++) {
        AMP::LinearAlgebra::Variable::shared_ptr dummyVar;
        d_inpVariables->add(dummyVar);
    }//end for i

    for (unsigned int var = 0; var < Diffusion::NUMBER_VARIABLES; var++) {
        if (d_isActive[var])
        {
            std::string name = activeVariables_db->getString(Diffusion::names[var]);
            AMP::LinearAlgebra::Variable::shared_ptr dummyVar(new AMP::LinearAlgebra::Variable(name));
            d_inpVariables->setVariable(var,dummyVar);
            if (d_isFrozen[var])
            {
                d_inVec[var] = d_Frozen[var];
                if(d_inVec[var]!=NULL) (d_inVec[var])->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
            }
        }
        else
        {
            AMP::LinearAlgebra::Variable::shared_ptr dummyVar;
            d_inpVariables->add(dummyVar);
        }
    }

    d_outVariable.reset(new AMP::LinearAlgebra::Variable(params->d_db->getString("OutputVariable")));

    init(params);
}

void DiffusionNonlinearFEOperator::preAssembly(
        const boost::shared_ptr<AMP::LinearAlgebra::Vector> &u, boost::shared_ptr<AMP::LinearAlgebra::Vector> &r)
{
  AMP_INSIST( (u != NULL), "NULL Input Vector!" );

  if( d_iDebugPrintInfoLevel > 7 )
    {
      AMP::pout << "DiffusionNonlinearFEOperator::preAssembly, entering" << std::endl;
    }
  
  for (unsigned int var = 0; var < Diffusion::NUMBER_VARIABLES; var++)
  {
    if (d_isActive[var])
    {
      if(d_isFrozen[var]) {
        d_inVec[var] = d_Frozen[var];
      } else {
        AMP::LinearAlgebra::Variable::shared_ptr tvar = d_inpVariables->getVariable(var);
        d_inVec[var] = u->subsetVectorForVariable(tvar);
        AMP_ASSERT(d_inVec[var]);
      }

      AMP_ASSERT(d_inVec[var]!=NULL);
      (d_inVec[var])->makeConsistent(AMP::LinearAlgebra::Vector::CONSISTENT_SET);
      if (d_iDebugPrintInfoLevel > 5)
        std::cout << "Max Value inside preAssembly: "
          << d_inVec[var]->max() << std::endl;
    }
  }

  d_outVec = r->subsetVectorForVariable(d_outVariable);
  d_outVec->zero();

  d_transportModel->preNonlinearAssembly();

  if( d_iDebugPrintInfoLevel > 7 )
  {
    AMP::pout << "DiffusionNonlinearFEOperator::preAssembly, leaving" << std::endl;
  }
}

void DiffusionNonlinearFEOperator::postAssembly()
{
  if( d_iDebugPrintInfoLevel > 7 )
  {
    AMP::pout << "DiffusionNonlinearFEOperator::postAssembly, entering" << std::endl;
  }

  d_transportModel->postNonlinearAssembly();
  d_outVec->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_ADD );

  if( d_iDebugPrintInfoLevel > 7 )
  {
    AMP::pout << "DiffusionNonlinearFEOperator::postAssembly, leaving" << std::endl;
  }
}

void DiffusionNonlinearFEOperator::preElementOperation(
    const AMP::Mesh::MeshElement & elem )
{
AMP_ERROR("Not converted yet");
/*
  if( d_iDebugPrintInfoLevel > 7 )
    {
      AMP::pout << "DiffusionNonlinearFEOperator::preElementOperation, entering" << std::endl;     
    }
  unsigned int num_local_Dofs = 0;
  AMP_ASSERT(dof_maps.size() == 1);
  (dof_maps[0])->getDOFs(elem, d_DofIndices);
  num_local_Dofs = d_DofIndices.size();
  
  std::vector<std::vector<double> > elementInputVectors(Diffusion::NUMBER_VARIABLES);

  d_numNodesForCurrentElement = elem.numNodes();
  
  for (unsigned int var = 0; var < Diffusion::NUMBER_VARIABLES; var++)
    {
      if (d_isActive[var])
    {
      elementInputVectors[var].resize(num_local_Dofs);
      for (unsigned int i = 0; i < d_numNodesForCurrentElement; i++)
        {
          elementInputVectors[var][i]
        = (d_inVec[var])->getValueByGlobalID(d_DofIndices[i]);
        }
    }
    }
  
  d_elementOutputVector.resize(num_local_Dofs);
  for (unsigned int i = 0; i < num_local_Dofs; i++)
    {
      d_elementOutputVector[i] = 0.0;
    }
  
  d_diffNonlinElem->setElementVectors(elementInputVectors,
                      d_elementOutputVector);
  
  const ::Elem* elemPtr = &(elem.getElem());
  d_diffNonlinElem->initializeForCurrentElement(elemPtr, d_transportModel);

  if( d_iDebugPrintInfoLevel > 7 )
    {
      AMP::pout << "DiffusionNonlinearFEOperator::preElementOperation, leaving" << std::endl;     
    }
*/
}

void DiffusionNonlinearFEOperator::postElementOperation()
{
  if( d_iDebugPrintInfoLevel > 7 )
    {
      AMP::pout << "DiffusionNonlinearFEOperator::postElementOperation, entering" << std::endl;
    }
  
  for (unsigned int i = 0; i < d_numNodesForCurrentElement; i++)
    {
      d_outVec->addValueByGlobalID(d_DofIndices[i],
                   d_elementOutputVector[i]);
    }

  if( d_iDebugPrintInfoLevel > 7 )
    {
      AMP::pout << "DiffusionNonlinearFEOperator::postElementOperation, leaving" << std::endl;
    }
}

void DiffusionNonlinearFEOperator::init(const boost::shared_ptr<
        DiffusionNonlinearFEOperatorParameters>& params)
{
  if( d_iDebugPrintInfoLevel > 7 )
    {
      AMP::pout << "DiffusionNonlinearFEOperator::init, entering" << std::endl;
    }
  
  (void) params;
  AMP::Mesh::MeshIterator el = d_Mesh->getIterator(AMP::Mesh::Volume,0);
  AMP::Mesh::MeshIterator end_el = el.end();
  
AMP_ERROR("Not converted yet");
/*
  for (; el != end_el; ++el)
    {
      const ::Elem* elemPtr = &(el->getElem());
      d_diffNonlinElem->initializeForCurrentElement(elemPtr, d_transportModel);
      d_diffNonlinElem->initTransportModel();
    }//end for el
*/  

  if( d_iDebugPrintInfoLevel > 7 )
    {
      AMP::pout << "DiffusionNonlinearFEOperator::init, leaving" << std::endl;
    }
}

void DiffusionNonlinearFEOperator::reset(
        const boost::shared_ptr<OperatorParameters>& params)
{
    boost::shared_ptr<DiffusionNonlinearFEOperatorParameters> dnlparams_sp =
            boost::dynamic_pointer_cast<DiffusionNonlinearFEOperatorParameters,
                                        OperatorParameters>(params);

    if (d_PrincipalVariable == Diffusion::TEMPERATURE)
        d_inVec[d_PrincipalVariable] = dnlparams_sp->d_FrozenTemperature;
    if (d_PrincipalVariable == Diffusion::CONCENTRATION)
        d_inVec[d_PrincipalVariable] = dnlparams_sp->d_FrozenConcentration;
    if (d_PrincipalVariable == Diffusion::BURNUP)
        d_inVec[d_PrincipalVariable] = dnlparams_sp->d_FrozenBurnup;
    (d_inVec[d_PrincipalVariable])->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );

    resetFrozen(dnlparams_sp);
    for (unsigned int var = 0; var < Diffusion::NUMBER_VARIABLES; var++) {
        if (d_isActive[var])
        {
            if (d_isFrozen[var])
            {
                d_inVec[var] = d_Frozen[var];
                (d_inVec[var])->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
            }
        }
    }
}

boost::shared_ptr<OperatorParameters> DiffusionNonlinearFEOperator::getJacobianParameters(
        const boost::shared_ptr<AMP::LinearAlgebra::Vector>& u)
{
    boost::shared_ptr<AMP::InputDatabase> tmp_db(new AMP::InputDatabase("Dummy"));

    // set up a database for the linear operator params
    tmp_db->putString("name", "DiffusionLinearFEOperator");
    tmp_db->putString("InputVariable", Diffusion::names[d_PrincipalVariable]);
    tmp_db->putString("OutputVariable", d_outVariable->getName());
    tmp_db->putBool("FixedTemperature", d_isActive[Diffusion::TEMPERATURE]?false:true);
    tmp_db->putBool("FixedConcentration", d_isActive[Diffusion::CONCENTRATION]?false:true);
    tmp_db->putBool("FixedBurnup", d_isActive[Diffusion::BURNUP]?false:true);

    // create the linear operator params
    boost::shared_ptr<DiffusionLinearFEOperatorParameters> outParams(
            new DiffusionLinearFEOperatorParameters(tmp_db));

    // create the linear element object
    boost::shared_ptr<AMP::InputDatabase> elem_db(new AMP::InputDatabase("Dummy"));
    tmp_db->putBool("TransportAtGaussPoints", d_diffNonlinElem->getTransportAtGauss());
    boost::shared_ptr<ElementOperationParameters> eparams(new ElementOperationParameters(elem_db));
    boost::shared_ptr<DiffusionLinearElement> linearElement(new DiffusionLinearElement(eparams));

    // add miscellaneous to output parameters
    outParams->d_transportModel = d_transportModel;
    outParams->d_elemOp = linearElement;
    outParams->d_Mesh = d_Mesh;

    // add variables to parameters
    if (d_isActive[Diffusion::TEMPERATURE]) {
        if (d_isFrozen[Diffusion::TEMPERATURE]) {
            outParams->d_temperature = d_Frozen[Diffusion::TEMPERATURE];
        } else {
            AMP::LinearAlgebra::Variable::shared_ptr tvar = d_inpVariables->getVariable(Diffusion::TEMPERATURE);
            AMP::LinearAlgebra::Vector::shared_ptr temperature = u->subsetVectorForVariable(tvar);
            outParams->d_temperature = temperature;
            outParams->d_temperature->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
        }
    }

    if (d_isActive[Diffusion::CONCENTRATION]) {
        if (d_isFrozen[Diffusion::CONCENTRATION]) {
            outParams->d_concentration = d_Frozen[Diffusion::CONCENTRATION];
        } else {
            AMP::LinearAlgebra::Variable::shared_ptr cvar = d_inpVariables->getVariable(Diffusion::CONCENTRATION);
            AMP::LinearAlgebra::Vector::shared_ptr concentration = u->subsetVectorForVariable(cvar);
            outParams->d_concentration = concentration;
            outParams->d_concentration->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
        }
    }

    if (d_isActive[Diffusion::BURNUP]) {
        if (d_isFrozen[Diffusion::BURNUP]) {
            outParams->d_burnup = d_Frozen[Diffusion::BURNUP];
        } else {
            AMP::LinearAlgebra::Variable::shared_ptr bvar = d_inpVariables->getVariable(Diffusion::BURNUP);
            AMP::LinearAlgebra::Vector::shared_ptr burnup = u->subsetVectorForVariable(bvar);
            outParams->d_burnup = burnup;
            outParams->d_burnup->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
        }
    }

    return outParams;
}

void DiffusionNonlinearFEOperator::resetFrozen(
        const boost::shared_ptr<DiffusionNonlinearFEOperatorParameters> &params)
{
    using namespace Diffusion;
    for (size_t var=0; var<Diffusion::NUMBER_VARIABLES; var++) {
        if (d_isActive[var] and d_isFrozen[var]) {
            if (var == TEMPERATURE)   d_Frozen[var]=params->d_FrozenTemperature;
            if (var == CONCENTRATION) d_Frozen[var]=params->d_FrozenConcentration;
            if (var == BURNUP)           d_Frozen[var]=params->d_FrozenBurnup;
        }
    }
}

bool DiffusionNonlinearFEOperator::isValidInput(AMP::LinearAlgebra::Vector::shared_ptr &u)
{
    AMP::Materials::PropertyPtr property = d_transportModel->getProperty();
    std::vector<std::string> names = property->get_arguments();
    size_t nnames = names.size();
    std::string argname;
    bool found = false;
    if (d_PrincipalVariable == Diffusion::TEMPERATURE) {
        for (size_t i=0; i<nnames; i++) {
            if (names[i] == "temperature") {
                argname = "temperature";
                found = true;
                break;
            }
        }
    }
    if (d_PrincipalVariable == Diffusion::CONCENTRATION) {
        for (size_t i=0; i<nnames; i++) {
            if (names[i] == "concentration") {
                argname = "concentration";
                found = true;
                break;
            }
        }
    }

    bool result = true;
    if (found) {
		AMP::LinearAlgebra::Vector::shared_ptr uinp = u->subsetVectorForVariable(d_inpVariables->getVariable(d_PrincipalVariable));
		size_t nvals = 0, nit=0;
		for (AMP::LinearAlgebra::Vector::iterator val = uinp->begin(); val != uinp->end(); val++) {
			nvals++;
		}
		std::vector<double> vals(nvals);
		for (AMP::LinearAlgebra::Vector::iterator val = uinp->begin(); val != uinp->end(); val++) {
			vals[nit] = *val;
			nit++;
		}
		result = property->in_range(argname, vals);
    }

    return result;
}


}
}//end namespace
