#include "RobinVectorCorrection.h"
#include "RobinMatrixCorrectionParameters.h"
#include "utils/Utilities.h"
#include "utils/ProfilerApp.h"
#include "utils/InputDatabase.h"

/* Libmesh files */

#include "enum_order.h"
#include "enum_fe_family.h"
#include "enum_quadrature_type.h"
#include "auto_ptr.h"
#include "string_to_enum.h"

#include <string>

namespace AMP {
namespace Operator {

void RobinVectorCorrection::reset(const boost::shared_ptr<OperatorParameters>& params)
{
  NeumannVectorCorrection::reset(params);
  
  AMP_INSIST( ((params.get()) != NULL), "NULL parameters" );
  AMP_INSIST( (((params->d_db).get()) != NULL), "NULL database" );
  
  AMP_INSIST( (params->d_db)->keyExists("fConductivity"), "Missing key: effective convective coefficient" );
  d_hef = (params->d_db)->getDouble("fConductivity");
  
  d_skipParams = (params->d_db)->getBoolWithDefault("skip_params", true);
  
  AMP_INSIST( (params->d_db)->keyExists("alpha"), "Missing key: prefactor alpha" );
  d_alpha = (params->d_db)->getDouble("alpha");
  
  AMP_INSIST( (params->d_db)->keyExists("beta"), "Missing key: prefactor beta" );
  d_beta = (params->d_db)->getDouble("beta");
  
  AMP_INSIST( (params->d_db)->keyExists("gamma"), "Missing key: total prefactor gamma" );
  d_gamma = (params->d_db)->getDouble("gamma");
  
}
  
void
RobinVectorCorrection::apply(const AMP::LinearAlgebra::Vector::shared_ptr &f,
                 const AMP::LinearAlgebra::Vector::shared_ptr &u,
                 AMP::LinearAlgebra::Vector::shared_ptr &r,
                 const double a,
                 const double b)
{
  PROFILE_START("apply");
  AMP_INSIST( ((r.get()) != NULL), "NULL Residual Vector" );
  AMP_INSIST( ((u.get()) != NULL), "NULL Solution Vector" );

  AMP::LinearAlgebra::Vector::shared_ptr rInternal = this->subsetInputVector(r);
  AMP::LinearAlgebra::Vector::shared_ptr uInternal = this->subsetInputVector(u);

  uInternal->makeConsistent ( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
  //rInternal->makeConsistent ( AMP::LinearAlgebra::Vector::CONSISTENT_SET );

  std::vector<std::string> variableNames;
  size_t numVar = 0 ;
  if(d_robinPhysicsModel.get() != NULL)
  {
    variableNames = d_robinPhysicsModel->getVariableName();
    numVar = variableNames.size();
  }

  d_elementInputVec.resize( numVar + 1);
  d_elementInputVec[0] = d_variableFlux;

  if(d_robinPhysicsModel.get() != NULL)
  {

    for (size_t i=0; i<variableNames.size(); i++)
    {
      std::string cview = variableNames[i] + " view";
      if(d_Frozen.get() != NULL)
      {
        if( d_Frozen->select ( AMP::LinearAlgebra::VS_ByVariableName ( variableNames[i] ) , cview ) != NULL )
        {
          d_elementInputVec[i+1] = d_Frozen->select ( AMP::LinearAlgebra::VS_ByVariableName ( variableNames[i] ) , cview );
        }
        else
        {
          d_elementInputVec[i+1] = uInternal->select ( AMP::LinearAlgebra::VS_ByVariableName ( variableNames[i] ) , cview );
        }
      }
      else
      {
        d_elementInputVec[i+1] = uInternal->select ( AMP::LinearAlgebra::VS_ByVariableName ( variableNames[i] ) , cview );
      }
      AMP_INSIST ( d_elementInputVec[i+1] , "Did not find vector" );
      (d_elementInputVec[i+1])->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
    }

    //#define DEBUG_GAP_PRINT
#ifdef DEBUG_GAP_PRINT
    if (d_iDebugPrintInfoLevel==100)
    {
      std::cout << "processing robin boundary operator " << d_InstanceID << "\n";
    }
#endif

  }

  AMP::Discretization::DOFManager::shared_ptr dofManager = rInternal->getDOFManager();
  if(d_robinPhysicsModel.get() != NULL) {
//    for(unsigned int m = 0; m < d_elementInputVec.size(); m++)
//      AMP_ASSERT(*dofManager==*(d_elementInputVec[m]->getDOFManager()));
  }
  AMP_ASSERT(*dofManager==*(uInternal->getDOFManager()));
  AMP_ASSERT(*dofManager==*(rInternal->getDOFManager()));
  if ( d_variableFlux.get()!=NULL && !d_isFluxGaussPtVector )
    AMP_ASSERT(*dofManager==*(d_variableFlux->getDOFManager()));

  std::vector<size_t> gpDofs;
  AMP::Discretization::DOFManager::shared_ptr gpDOFManager; 
  if( !d_isConstantFlux && d_isFluxGaussPtVector){
    gpDOFManager = d_variableFlux->getDOFManager();
  }

  unsigned int numIds = d_boundaryIds.size();
  std::vector<size_t> dofs;
  PROFILE_START("integration loop");
  for (unsigned int nid = 0; nid < numIds; nid++)
  {

    AMP::Mesh::MeshIterator bnd1     = d_Mesh->getBoundaryIDIterator( AMP::Mesh::Face, d_boundaryIds[nid], 0 );
    AMP::Mesh::MeshIterator end_bnd1 = bnd1.end();

    for (; bnd1 != end_bnd1; ++bnd1)
    {
      PROFILE_START("prepare element",2);

      boost::shared_ptr < ::FEType > d_feType ( new ::FEType(d_feTypeOrder, d_feFamily) );
      boost::shared_ptr < ::FEBase > d_fe( (::FEBase::build(2, (*d_feType))).release() );

      if(d_qruleOrderName == "DEFAULT") {
        d_qruleOrder = d_feType->default_quadrature_order();
      } else {
        d_qruleOrder = Utility::string_to_enum<libMeshEnums::Order>(d_qruleOrderName);
      }
      boost::shared_ptr < ::QBase > d_qrule( (::QBase::build(d_qruleType, 2, d_qruleOrder)).release() );

      // Get the nodes for the current element
      d_currNodes = bnd1->getElements(AMP::Mesh::Vertex);
      unsigned int numNodesInCurrElem = d_currNodes.size();

      // Get the dofs for the vectors
      std::vector<AMP::Mesh::MeshElementID> ids(d_currNodes.size());
      for (size_t i=0; i<d_currNodes.size(); i++)
        ids[i] = d_currNodes[i].globalID();
      dofManager->getDOFs( ids, dofs );
      AMP_ASSERT(dofs.size()==numNodesInCurrElem);

      if( !d_isConstantFlux && d_isFluxGaussPtVector){
        gpDOFManager->getDOFs (bnd1->globalID(), gpDofs);
      }
      // Get the libmesh element
      d_currElemPtr = libmeshElements.getElement( bnd1->globalID() );

      d_fe->attach_quadrature_rule( d_qrule.get() );

      d_phi = &(d_fe->get_phi());
      d_JxW = &(d_fe->get_JxW());

      d_fe->reinit ( d_currElemPtr );

      const std::vector<Real> & JxW = (*d_JxW);
      const std::vector<std::vector<Real> > & phi = (*d_phi);
      unsigned int numGaussPts = d_qrule->n_points(); 
      PROFILE_STOP("prepare element",2);

      std::vector<std::vector<double> > inputArgs(d_elementInputVec.size(),std::vector<double>(numNodesInCurrElem));
      std::vector<std::vector<double> > inputArgsAtGpts(d_elementInputVec.size(),std::vector<double>(numGaussPts));
      std::vector<double> beta(numGaussPts,d_beta);
      std::vector<double> gamma(numGaussPts,d_gamma);
      PROFILE_START("get conductance",2);
      if(d_robinPhysicsModel.get() != NULL) 
      {
        unsigned int startIdx = 0;
        if(d_isFluxGaussPtVector){
          startIdx = 1;
        }
        for(unsigned int m = startIdx; m < d_elementInputVec.size(); m++){
          d_elementInputVec[m]->getValuesByGlobalID( dofs.size(), &dofs[0], &inputArgs[m][0] );
          for (size_t qp = 0; qp < numGaussPts; qp++){ 
            for (size_t n = 0; n < numNodesInCurrElem ; n++) {
              inputArgsAtGpts[m][qp] += phi[n][qp] * inputArgs[m][n];
            }
          }
        }
        
        d_robinPhysicsModel->getConductance(beta, gamma, inputArgsAtGpts);
      }
      PROFILE_STOP("get conductance",2);
      PROFILE_START("perform integration",2);
      std::vector<double> values(dofs.size(),0.0);
      std::vector<double> gpValues(gpDofs.size(),0.0);
      std::vector<double> addValues(dofs.size(),0.0);
      uInternal->getValuesByGlobalID( dofs.size(), &dofs[0], &values[0] );
      for (unsigned int qp = 0; qp < d_qrule->n_points(); qp++)
      {
        Real phi_val = 0.0;
        for (unsigned int l = 0; l < numNodesInCurrElem ; l++)
          phi_val += phi[l][qp] * values[l];
        for (unsigned int j = 0; j < numNodesInCurrElem ; j++)
          addValues[j] += (JxW[qp] * phi[j][qp] * beta[qp] * phi_val);
      }//end for qp
      if (d_IsCoupledBoundary[nid])
      {
        if(!d_isFluxGaussPtVector){
          d_variableFlux->getValuesByGlobalID( dofs.size(), &dofs[0], &values[0] );
        }else{
          d_variableFlux->getValuesByGlobalID( gpDofs.size(), &gpDofs[0], &gpValues[0] );
        }
        for (unsigned int qp = 0; qp < d_qrule->n_points(); qp++)
        {
          Real phi_val = 0.0;
          if(!d_isFluxGaussPtVector){
            for (unsigned int l = 0; l < numNodesInCurrElem ; l++)
              phi_val += phi[l][qp] * values[l];
          }else{
            phi_val =  gpValues[qp];
          }
          for (unsigned int j = 0; j < numNodesInCurrElem ; j++)
            addValues[j] += -1 * (JxW[qp] * phi[j][qp] * gamma[qp] * phi_val);
        }//end for qp
      }//coupled
      rInternal->addValuesByGlobalID( dofs.size(), &dofs[0], &addValues[0] );
      PROFILE_STOP("perform integration",2);

    }//end for bnd
  }//end for nid
  PROFILE_STOP("integration loop");

  rInternal->makeConsistent(AMP::LinearAlgebra::Vector::CONSISTENT_ADD);
  //std::cout << rInternal << std::endl;

  if (f.get() == NULL)
  {
    rInternal->scale(a);
  }
  else
  {
    AMP::LinearAlgebra::Vector::shared_ptr fInternal = this->subsetOutputVector(f);
    if (fInternal.get() == NULL)
    {
      rInternal->scale(a);
    }
    else
    {

      rInternal->axpby(b, a, fInternal);
    }
  }
  PROFILE_STOP("apply");
}

boost::shared_ptr<OperatorParameters> RobinVectorCorrection::getJacobianParameters(const boost::shared_ptr<AMP::LinearAlgebra::Vector>&)
{
  boost::shared_ptr<AMP::InputDatabase> tmp_db(new AMP::InputDatabase("Dummy"));
  tmp_db->putBool("skip_params", d_skipParams);
  tmp_db->putBool("skip_rhs_correction", true);
  tmp_db->putBool("skip_matrix_correction", false);

  if (!d_skipParams)
  {
    tmp_db->putDouble("alpha", d_alpha);
    tmp_db->putDouble("beta", d_beta);
    tmp_db->putDouble("gamma", d_gamma);
    tmp_db->putDouble("fConductivity", d_hef);

    int numIds = d_boundaryIds.size();
    tmp_db->putInteger("number_of_ids", numIds);

    char key[100];
    for (int j = 0; j < numIds; j++)
    {

      sprintf(key, "id_%d", j);
      tmp_db->putInteger(key, d_boundaryIds[j]);

      sprintf(key, "number_of_dofs_%d", j);
      int numDofIds = d_dofIds[j].size();
      tmp_db->putInteger(key, numDofIds);
      
      for (int i = 0; i < numDofIds; i++)
        {
          sprintf(key, "dof_%d_%d", j, i);
          tmp_db->putInteger(key, d_dofIds[j][i]);
          
          sprintf(key, "value_%d_%d", j, i);
          tmp_db->putDouble(key, d_neumannValues[j][i]);
        }
    }
    }

  boost::shared_ptr<RobinMatrixCorrectionParameters> outParams(
      new RobinMatrixCorrectionParameters(tmp_db));

  outParams->d_robinPhysicsModel = d_robinPhysicsModel;
  outParams->d_elementInputVec   = d_elementInputVec;
  
  return outParams;
}

}
}

