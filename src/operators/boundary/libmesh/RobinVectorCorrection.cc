#include "RobinVectorCorrection.h"
#include "RobinMatrixCorrectionParameters.h"
#include "utils/Utilities.h"
#include "utils/ProfilerApp.h"
#include "utils/InputDatabase.h"

/* Libmesh files */

#include "libmesh/enum_order.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/auto_ptr.h"
#include "libmesh/string_to_enum.h"

#include <string>

namespace AMP {
namespace Operator {

void RobinVectorCorrection::reset(const boost::shared_ptr<OperatorParameters>& params)
{
    NeumannVectorCorrection::reset(params);

    AMP_INSIST( ((params.get()) != NULL), "NULL parameters" );
    AMP_INSIST( (((params->d_db).get()) != NULL), "NULL database" );

    d_skipParams = (params->d_db)->getBoolWithDefault("skip_params", false);

    AMP_INSIST( (params->d_db)->keyExists("alpha"), "Missing key: prefactor alpha" );
    d_alpha = (params->d_db)->getDouble("alpha");

    AMP_INSIST( (params->d_db)->keyExists("beta"), "Missing key: prefactor beta" );
    d_beta = (params->d_db)->getDouble("beta");

    AMP_INSIST( (params->d_db)->keyExists("gamma"), "Missing key: total prefactor gamma" );
    d_gamma = (params->d_db)->getDouble("gamma");
}

  
void RobinVectorCorrection::apply(AMP::LinearAlgebra::Vector::const_shared_ptr f,
                 AMP::LinearAlgebra::Vector::const_shared_ptr u,
                 AMP::LinearAlgebra::Vector::shared_ptr r,
                 const double a,
                 const double b)
{
    PROFILE_START("apply");
    AMP_INSIST( ((r.get()) != NULL), "NULL Residual Vector" );
    AMP_INSIST( ((u.get()) != NULL), "NULL Solution Vector" );

    AMP::LinearAlgebra::Vector::shared_ptr rInternal = this->subsetInputVector(r);
    AMP::LinearAlgebra::Vector::const_shared_ptr uInternal = this->subsetInputVector(u);

    AMP_ASSERT(uInternal->getUpdateStatus()==AMP::LinearAlgebra::Vector::UNCHANGED);
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

    if(d_robinPhysicsModel.get() != NULL) {
        for (size_t i=0; i<variableNames.size(); i++) {
            std::string cview = variableNames[i] + " view";
            if(d_Frozen.get() != NULL) {
                if( d_Frozen->select ( AMP::LinearAlgebra::VS_ByVariableName ( variableNames[i] ) , cview ) != NULL ) {
                    d_elementInputVec[i+1] = d_Frozen->constSelect ( AMP::LinearAlgebra::VS_ByVariableName ( variableNames[i] ) , cview );
                } else {
                    d_elementInputVec[i+1] = uInternal->constSelect ( AMP::LinearAlgebra::VS_ByVariableName ( variableNames[i] ) , cview );
                }
            } else {
                d_elementInputVec[i+1] = uInternal->constSelect ( AMP::LinearAlgebra::VS_ByVariableName ( variableNames[i] ) , cview );
            }
            AMP_INSIST ( d_elementInputVec[i+1] , "Did not find vector '"+variableNames[i]+"'" );
            AMP_ASSERT(d_elementInputVec[i+1]->getUpdateStatus()==AMP::LinearAlgebra::Vector::UNCHANGED);
        }

        if (d_iDebugPrintInfoLevel==100)
        {
            std::cout << "processing robin boundary operator " << d_InstanceID << "\n";
        }
    }

    // Get the DOF managers
    AMP::Discretization::DOFManager::shared_ptr dofManager = rInternal->getDOFManager();
    AMP::Discretization::DOFManager::shared_ptr gpDOFManager; 
    if ( d_isFluxGaussPtVector && d_variableFlux!=NULL )
        gpDOFManager = d_variableFlux->getDOFManager();

    // Check that the DOF managers match for the different vectors
    AMP_ASSERT(*dofManager==*(uInternal->getDOFManager()));
    AMP_ASSERT(*dofManager==*(rInternal->getDOFManager()));
    if ( !d_isFluxGaussPtVector ) {
        if ( d_variableFlux.get()!= NULL  )
            AMP_ASSERT(*dofManager==*(d_variableFlux->getDOFManager()));
    }

    unsigned int numIds = d_boundaryIds.size();
    std::vector<size_t> gpDofs;
    std::vector<size_t> dofs;
    std::vector<std::vector<size_t> > dofIndices;
    std::vector<size_t> dofsElementVec;
    PROFILE_START("integration loop");
    for (unsigned int nid = 0; nid < numIds; nid++)
    {
        unsigned int numDofIds = d_dofIds[nid].size();

        for(unsigned int k = 0; k < numDofIds; k++)
        {

          AMP::Mesh::MeshIterator bnd1     = d_Mesh->getBoundaryIDIterator( AMP::Mesh::Face, d_boundaryIds[nid], 0 );
          const AMP::Mesh::MeshIterator end_bnd1 = bnd1.end();

          for (; bnd1 != end_bnd1; ++bnd1)
          {
            PROFILE_START("prepare element",2);

            // Get the nodes for the current element
            d_currNodes = bnd1->getElements(AMP::Mesh::Vertex);
            unsigned int numNodesInCurrElem = d_currNodes.size();

            dofIndices.resize(numNodesInCurrElem);
            // Get the dofs for the vectors
            std::vector<AMP::Mesh::MeshElementID> ids(d_currNodes.size());
            for (size_t i=0; i<d_currNodes.size(); i++)
              ids[i] = d_currNodes[i].globalID();

            for(unsigned int i = 0; i < numNodesInCurrElem ; i++) 
              dofManager->getDOFs(d_currNodes[i].globalID(), dofIndices[i]);

            for (size_t n = 0; n < dofIndices.size() ; n++)
              dofs[n] = dofIndices[n][d_dofIds[nid][k]];

            AMP_ASSERT(dofs.size()==numNodesInCurrElem);

            if(d_isFluxGaussPtVector && d_IsCoupledBoundary[nid]){
              gpDOFManager->getDOFs (bnd1->globalID(), gpDofs);
              AMP_ASSERT(gpDofs.size()>0);
            }

            // Get the current libmesh element
            const libMesh::FEBase* fe = d_libmeshElements.getFEBase( bnd1->globalID() );
            const libMesh::QBase* rule = d_libmeshElements.getQBase( bnd1->globalID() );
            AMP_ASSERT(fe!=NULL);
            AMP_ASSERT(rule!=NULL);
            const unsigned int numGaussPts = rule->n_points(); 

            const std::vector<Real> JxW = fe->get_JxW();
            const std::vector<std::vector<Real> > phi = fe->get_phi();
            PROFILE_STOP("prepare element",2);

            std::vector<std::vector<double> > inputArgs(d_elementInputVec.size(),std::vector<double>(numNodesInCurrElem));
            std::vector<std::vector<double> > inputArgsAtGpts(d_elementInputVec.size(),std::vector<double>(numGaussPts));
            std::vector<double> beta(numGaussPts,d_beta);
            std::vector<double> gamma(numGaussPts,d_gamma);
            PROFILE_START("get conductance",2);
            if(d_robinPhysicsModel.get() != NULL) 
            {
              unsigned int startIdx = 0;
              if(d_isFluxGaussPtVector && d_IsCoupledBoundary[nid]){
                d_variableFlux->getValuesByGlobalID( gpDofs.size(), &gpDofs[0], &inputArgsAtGpts[0][0] );
                startIdx = 1;
              }
              for(unsigned int m = startIdx; m < d_elementInputVec.size(); m++){
                // Note: elementInputVecs may use different DOFManagers from u and r internal
                d_elementInputVec[m]->getDOFManager()->getDOFs( ids, dofsElementVec );
                AMP_ASSERT(dofsElementVec.size()==dofs.size());
                d_elementInputVec[m]->getValuesByGlobalID( dofs.size(), &dofsElementVec[0], &inputArgs[m][0] );
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
            for (unsigned int qp = 0; qp<numGaussPts; qp++)
            {
              Real phi_val = 0.0;
              for (unsigned int l = 0; l < numNodesInCurrElem ; l++)
                phi_val += phi[l][qp] * values[l];
              for (unsigned int j = 0; j < numNodesInCurrElem ; j++)
                addValues[j] += (JxW[qp] * phi[j][qp] * beta[qp] * phi_val);
            }//end for qp
            if ( d_IsCoupledBoundary[nid] )
            {
              if ( !d_isFluxGaussPtVector ){
                d_variableFlux->getValuesByGlobalID( dofs.size(), &dofs[0], &values[0] );
              } else {
                d_variableFlux->getValuesByGlobalID( gpDofs.size(), &gpDofs[0], &gpValues[0] );
              }
              for (unsigned int qp = 0; qp<numGaussPts; qp++)
              {
                Real phi_val = 0.0;
                if ( !d_isFluxGaussPtVector ){
                  for (unsigned int l = 0; l < numNodesInCurrElem ; l++)
                    phi_val += phi[l][qp] * values[l];
                } else {
                  phi_val =  gpValues[qp];
                }
                for (unsigned int j = 0; j < numNodesInCurrElem ; j++)
                  addValues[j] += -1 * (JxW[qp] * phi[j][qp] * gamma[qp] * phi_val);
              }//end for qp
            }//coupled
            rInternal->addValuesByGlobalID( dofs.size(), &dofs[0], &addValues[0] );
            PROFILE_STOP("perform integration",2);

          }//end for bnd
        }//end for dof ids
    }//end for nid
    PROFILE_STOP("integration loop");

    rInternal->makeConsistent(AMP::LinearAlgebra::Vector::CONSISTENT_ADD);
    //std::cout << rInternal << std::endl;

    if (f.get() == NULL) {
      rInternal->scale(a);
    } else {
      AMP::LinearAlgebra::Vector::const_shared_ptr fInternal = this->subsetOutputVector(f);
      if (fInternal.get() == NULL) {
        rInternal->scale(a);
      } else {
        rInternal->axpby(b, a, fInternal);
      }
    }
    PROFILE_STOP("apply");
}


boost::shared_ptr<OperatorParameters> RobinVectorCorrection::getJacobianParameters(const boost::shared_ptr<AMP::LinearAlgebra::Vector>&)
{
  boost::shared_ptr<AMP::InputDatabase> tmp_db(new AMP::InputDatabase("Dummy"));
  tmp_db->putBool("skip_params", true);
  tmp_db->putBool("skip_rhs_correction", true);
  tmp_db->putBool("skip_matrix_correction", false);
  tmp_db->putBool("IsFluxGaussPtVector", d_isFluxGaussPtVector );
  boost::shared_ptr<RobinMatrixCorrectionParameters> outParams(
      new RobinMatrixCorrectionParameters(tmp_db));

  outParams->d_robinPhysicsModel = d_robinPhysicsModel;
  outParams->d_elementInputVec   = d_elementInputVec;
  outParams->d_variableFlux      = d_variableFlux;

  return outParams;
}


}
}

