#include "DiffusionLinearFEOperator.h"
#include "utils/Utilities.h"
#include "ampmesh/MeshVariable.h"

namespace AMP {
namespace Operator {

DiffusionLinearFEOperator::DiffusionLinearFEOperator(const boost::shared_ptr<
        DiffusionLinearFEOperatorParameters> & params) :
    LinearFEOperator(params) {
    AMP_INSIST( ((params.get()) != NULL), "NULL parameter" );

    d_diffLinElem = boost::dynamic_pointer_cast<DiffusionLinearElement>(d_elemOp);

    AMP_INSIST( ((d_diffLinElem.get()) != NULL), "d_elemOp is not of type DiffusionLinearElement" );

    d_useConstantTemperature = params->d_db->getBoolWithDefault("FixedTemperature", false);
    d_useConstantConcentration = params->d_db->getBoolWithDefault("FixedConcentration", false);
    d_useConstantBurnup = params->d_db->getBoolWithDefault("FixedBurnup", false);

    std::string inpVar = params->d_db->getString("InputVariable");
    d_inpVariable.reset(new AMP::Mesh::NodalScalarVariable(inpVar,d_MeshAdapter));

    std::string outVar = params->d_db->getString("OutputVariable");
    d_outVariable.reset(new AMP::Mesh::NodalScalarVariable(outVar,d_MeshAdapter));

    reset(params);
}

void DiffusionLinearFEOperator::preAssembly(const boost::shared_ptr<
        OperatorParameters>& oparams) {
    boost::shared_ptr<DiffusionLinearFEOperatorParameters> params =
        boost::dynamic_pointer_cast<DiffusionLinearFEOperatorParameters>(oparams);

  if( d_iDebugPrintInfoLevel > 7 )
    {
      AMP::pout << "DiffusionLinearFEOperator::preAssembly, entering" << std::endl;
    }

  d_transportModel = params->d_transportModel;
  
  if(d_temperature.get() == NULL and params->d_temperature.get() != NULL) {
    d_temperature = params->d_temperature->cloneVector();
  }
  if(d_temperature.get() != NULL){
    if(params->d_temperature.get() != NULL) {
      d_temperature->copyVector(params->d_temperature);
      d_temperature->makeConsistent(AMP::LinearAlgebra::Vector::CONSISTENT_SET);
    }
    else {
      d_temperature.reset();
    }
//    std::cout << d_temperature << std::endl;
  }
  
  if(d_concentration.get() == NULL and params->d_concentration.get() != NULL) {
    d_concentration = params->d_concentration->cloneVector();
  }
  if(d_concentration.get() != NULL){
    if(params->d_concentration.get() != NULL) {
      d_concentration->copyVector(params->d_concentration);
      d_concentration->makeConsistent(AMP::LinearAlgebra::Vector::CONSISTENT_SET);
    }
    else {
      d_concentration.reset();
    }
  }
  
  if(d_burnup.get() == NULL and params->d_burnup.get() != NULL) {
    d_burnup = params->d_burnup->cloneVector();
  }
  if(d_burnup.get() != NULL){
    if(params->d_burnup.get() != NULL) {
      d_burnup->copyVector(params->d_burnup);
      d_burnup->makeConsistent(AMP::LinearAlgebra::Vector::CONSISTENT_SET);
    }
    else {
      d_burnup.reset();
    }
  }
  
  d_matrix->zero();

  d_transportModel->preLinearAssembly();
  
  if( d_iDebugPrintInfoLevel > 7 )
    {
      AMP::pout << "DiffusionLinearFEOperator::preAssembly, leaving" << std::endl;
    }
}

void DiffusionLinearFEOperator::postAssembly() {
    
  if( d_iDebugPrintInfoLevel > 7 )
    {
      AMP::pout << "DiffusionLinearFEOperator::postAssembly, entering" << std::endl;
    }
  
    d_transportModel->postLinearAssembly();

    d_matrix->makeConsistent();

  if( d_iDebugPrintInfoLevel > 7 )
    {
      AMP::pout << "DiffusionLinearFEOperator::postAssembly, leaving" << std::endl;
    }
}

void DiffusionLinearFEOperator::preElementOperation(
        const AMP::Mesh::MeshManager::Adapter::Element & elem, const std::vector<
                AMP::Mesh::DOFMap::shared_ptr> & dof_maps) {

    if( d_iDebugPrintInfoLevel > 7 ) {
      AMP::pout << "DiffusionLinearFEOperator::preElementOperation, entering" << std::endl;
    }

    (dof_maps[0])->getDOFs(elem, d_dofIndices);
    unsigned int num_local_dofs = d_dofIndices.size();

    d_elementStiffnessMatrix.resize(num_local_dofs);
    for (unsigned int r = 0; r < num_local_dofs; r++) {
        d_elementStiffnessMatrix[r].resize(num_local_dofs);
        for (unsigned int c = 0; c < num_local_dofs; c++) {
            d_elementStiffnessMatrix[r][c] = 0;
        }
    }

    std::vector<double> localTemperature(num_local_dofs);
    std::vector<double> localConcentration(num_local_dofs);
    std::vector<double> localBurnup(num_local_dofs);

    if (d_useConstantTemperature or d_temperature.get() == NULL) {
        localTemperature.resize(0);
    } else {
            //    AMP::pout << d_temperature << std::endl;
        if (localTemperature.size() == 0) localTemperature.resize(num_local_dofs);
        for (unsigned int r = 0; r < num_local_dofs; r++) {
            localTemperature[r] = d_temperature->getValueByGlobalID(d_dofIndices[r]);
        }
    }

    if (d_useConstantConcentration or d_concentration.get() == NULL) {
        localConcentration.resize(0);
    } else {
        if (localConcentration.size() == 0) localConcentration.resize(num_local_dofs);
        for (unsigned int r = 0; r < num_local_dofs; r++) {
            localConcentration[r] = d_concentration->getValueByGlobalID(d_dofIndices[r]);
        }
    }

    if (d_useConstantBurnup or d_burnup.get() == NULL) {
        localBurnup.resize(0);
    } else {
        if (localBurnup.size() == 0) localBurnup.resize(num_local_dofs);
        for (unsigned int r = 0; r < num_local_dofs; r++) {
            localBurnup[r] = d_burnup->getValueByGlobalID(d_dofIndices[r]);
        }
    }

    const ::Elem* elemPtr = &(elem.getElem());

    d_diffLinElem->initializeForCurrentElement(elemPtr, d_transportModel);

    d_diffLinElem->setElementStiffnessMatrix(d_elementStiffnessMatrix);

    d_diffLinElem->setElementVectors(num_local_dofs, localTemperature, localConcentration,
            localBurnup);

  if( d_iDebugPrintInfoLevel > 7 )
    {
      AMP::pout << "DiffusionLinearFEOperator::preElementOperation, leaving" << std::endl;
    }
}

void DiffusionLinearFEOperator::postElementOperation() {

  if( d_iDebugPrintInfoLevel > 7 )
    {
      AMP::pout << "DiffusionLinearFEOperator::postElementOperation, entering" << std::endl;
    }
    unsigned int num_local_dofs = d_dofIndices.size();

    for (unsigned int r = 0; r < num_local_dofs; r++) {
        for (unsigned int c = 0; c < num_local_dofs; c++) {
            d_matrix->addValueByGlobalID(d_dofIndices[r], d_dofIndices[c],
                    d_elementStiffnessMatrix[r][c]);
        }
    }

  if( d_iDebugPrintInfoLevel > 7 )
    {
      AMP::pout << "DiffusionLinearFEOperator::postElementOperation, leaving" << std::endl;
    }
}

}
}//end namespace



