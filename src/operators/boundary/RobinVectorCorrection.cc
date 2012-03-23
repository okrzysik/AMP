#include "RobinVectorCorrection.h"
#include "RobinMatrixCorrectionParameters.h"
#include "utils/Utilities.h"
#include "ampmesh/MeshUtils.h"
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
  d_beta.resize(1);
  d_beta[0] = (params->d_db)->getDouble("beta");
  
  AMP_INSIST( (params->d_db)->keyExists("gamma"), "Missing key: total prefactor gamma" );
  d_gamma.resize(1);
  d_gamma[0] = (params->d_db)->getDouble("gamma");
  
}
  
void
RobinVectorCorrection::apply(const AMP::LinearAlgebra::Vector::shared_ptr &f,
                 const AMP::LinearAlgebra::Vector::shared_ptr &u,
                 AMP::LinearAlgebra::Vector::shared_ptr &r,
                 const double a,
                 const double b)
{
  AMP::Mesh::DOFMap::shared_ptr dof_map = d_MeshAdapter->getDOFMap(d_variable);

  AMP_INSIST( ((r.get()) != NULL), "NULL Residual Vector" );
  AMP_INSIST( ((u.get()) != NULL), "NULL Solution Vector" );

  AMP::LinearAlgebra::Vector::shared_ptr rInternal = r->subsetVectorForVariable(d_variable);
  AMP::LinearAlgebra::Vector::shared_ptr uInternal = u->subsetVectorForVariable(d_variable);

  AMP::LinearAlgebra::Vector::shared_ptr uOnMesh = u->select ( AMP::Mesh::VS_ByMesh ( d_MeshAdapter ) , u->getVariable()->getName() );
  rInternal = rInternal->select ( AMP::Mesh::VS_ByMesh ( d_MeshAdapter ) , rInternal->getVariable()->getName() );

  uInternal->makeConsistent ( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
  //rInternal->makeConsistent ( AMP::LinearAlgebra::Vector::CONSISTENT_SET );

  std::vector<std::string> variableNames;
  size_t numVar = 0 ;
  if(d_robinPhysicsModel.get() != NULL)
  {
    variableNames = d_robinPhysicsModel->getVariableName();
    numVar = variableNames.size();
  }

  std::vector<AMP::LinearAlgebra::Vector::shared_ptr> elementInputVec;
  elementInputVec.resize( numVar + 1);
  elementInputVec[0] = d_variableFlux;

  if(d_robinPhysicsModel.get() != NULL)
  {

    for (size_t i=0; i<variableNames.size(); i++)
    {
      std::string cview = variableNames[i] + " view";
      if(d_Frozen.get() != NULL)
      {
        if( d_Frozen->select ( AMP::LinearAlgebra::VS_ByVariableName ( variableNames[i] ) , cview ) != NULL )
        {
          elementInputVec[i+1] = d_Frozen->select ( AMP::LinearAlgebra::VS_ByVariableName ( variableNames[i] ) , cview );
        }
        else
        {
          elementInputVec[i+1] = uOnMesh->select ( AMP::LinearAlgebra::VS_ByVariableName ( variableNames[i] ) , cview );
        }
      }
      else
      {
        elementInputVec[i+1] = uOnMesh->select ( AMP::LinearAlgebra::VS_ByVariableName ( variableNames[i] ) , cview );
      }
      AMP_INSIST ( elementInputVec[i+1] , "Did not find vector" );
      (elementInputVec[i+1])->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
    }

    //#define DEBUG_GAP_PRINT
#ifdef DEBUG_GAP_PRINT
    if (d_iDebugPrintInfoLevel==100)
    {
      std::cout << "processing robin boundary operator " << d_InstanceID << "\n";
    }
#endif

  }

  unsigned int numIds = d_boundaryIds.size();


  for (unsigned int nid = 0; nid < numIds; nid++)
  {

    AMP::Mesh::MeshManager::Adapter::BoundarySideIterator bnd1 = d_MeshAdapter->beginSideBoundary(d_boundaryIds[nid]);
    AMP::Mesh::MeshManager::Adapter::BoundarySideIterator end_bnd1 = d_MeshAdapter->endSideBoundary(d_boundaryIds[nid]);

    for (; bnd1 != end_bnd1; ++bnd1)
    {

      AMP::Mesh::MeshManager::Adapter::Element cur_side = *bnd1;

      
      d_feType.reset( new ::FEType(d_feTypeOrder, d_feFamily) );
      d_fe.reset( (::FEBase::build(2, (*d_feType))).release() );
      
      d_phi = &(d_fe->get_phi());
      d_JxW = &(d_fe->get_JxW());
      
      d_qrule.reset( (::QBase::build(d_qruleType, 2, d_qruleOrder)).release() );
      d_fe->attach_quadrature_rule( d_qrule.get() );
      
      d_fe->reinit(&cur_side.getElem());

      const std::vector<Real> & JxW = (*d_JxW);
      const std::vector<std::vector<Real> > & phi = (*d_phi);

      std::vector<unsigned int> bndGlobalIds;
      dof_map->getDOFs(cur_side, bndGlobalIds);
      double temp;
      std::vector<std::vector<double> > inputArgs(elementInputVec.size()) ;

      for (unsigned int qp = 0; qp < d_qrule->n_points(); qp++)
      {
        temp = 0;
        Real phi_val = 0.0;

        for (unsigned int l = 0; l < bndGlobalIds.size(); l++)
        {
          if(d_robinPhysicsModel.get() != NULL)
          {
            for(unsigned int m = 0; m < elementInputVec.size(); m++)
            {
              inputArgs[m].resize(1);
              inputArgs[m][0] = ( elementInputVec[m]->getValueByGlobalID(bndGlobalIds[l]) );
            }
            d_robinPhysicsModel->getConductance(d_beta, d_gamma, inputArgs);
          }
          phi_val += phi[l][qp] * d_beta[0] * uInternal->getValueByGlobalID(bndGlobalIds[l]);
#ifdef DEBUG_GAP_PRINT
          if (d_iDebugPrintInfoLevel == 100)
          {
            std::cout << bndGlobalIds[l] << " " << qp << " " << d_beta[0] << " ";
            for (unsigned int m = 0; m < elementInputVec.size(); m++)
            {
              std::cout << inputArgs[m][0] << " ";
            }
            std::cout << "\n";
          }
#endif
        }
        for (unsigned int j = 0; j < bndGlobalIds.size(); j++)
        {
          temp = (JxW[qp] * phi[j][qp] * phi_val);
          rInternal->addValueByGlobalID(bndGlobalIds[j], temp);
        }//end for j
      }//end for qp

      if (d_IsCoupledBoundary[nid])
      {
        for (unsigned int qp = 0; qp < d_qrule->n_points(); qp++)
        {
          temp = 0;
          Real phi_val = 0.0;

          for (unsigned int l = 0; l < bndGlobalIds.size(); l++)
          {

            if(d_robinPhysicsModel.get() != NULL)
            {
              for(unsigned int m = 0; m < elementInputVec.size(); m++)
              {
                inputArgs[m][0] = ( elementInputVec[m]->getValueByGlobalID(bndGlobalIds[l]) );
              }
              d_robinPhysicsModel->getConductance(d_beta, d_gamma, inputArgs);
            }

            phi_val += phi[l][qp] * d_gamma[0] * d_variableFlux->getValueByGlobalID(bndGlobalIds[l]);
          }

          for (unsigned int j = 0; j < bndGlobalIds.size(); j++)
          {
            temp = (JxW[qp] * phi[j][qp] * phi_val);
            temp = temp * -1;
            rInternal->addValueByGlobalID(bndGlobalIds[j], temp);
          }//end for j
        }//end for qp
      }//coupled

    }//end for bnd
  }//end for nid

  rInternal->makeConsistent(AMP::LinearAlgebra::Vector::CONSISTENT_ADD);
  //std::cout << rInternal << std::endl;

  if (f.get() == NULL)
  {
    rInternal->scale(a);
  }
  else
  {
    AMP::LinearAlgebra::Vector::shared_ptr fInternal = f->subsetVectorForVariable(d_variable);
    if (fInternal.get() == NULL)
    {
      rInternal->scale(a);
    }
    else
    {

      rInternal->axpby(b, a, fInternal);
    }
  }
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
    tmp_db->putDouble("beta", d_beta[0]);
    tmp_db->putDouble("gamma", d_gamma[0]);
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

