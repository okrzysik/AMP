
#include "RobinMatrixCorrection.h"
#include "RobinMatrixCorrectionParameters.h"
#include "utils/Utilities.h"
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

RobinMatrixCorrection :: RobinMatrixCorrection(const boost::shared_ptr<RobinMatrixCorrectionParameters> & params)
  : BoundaryOperator (params)
{
  d_feTypeOrderName = (params->d_db)->getStringWithDefault("FE_ORDER", "FIRST");
  
  d_feTypeOrder = Utility::string_to_enum<libMeshEnums::Order>(d_feTypeOrderName);
  
  d_feFamilyName = (params->d_db)->getStringWithDefault("FE_FAMILY", "LAGRANGE");
  
  d_feFamily = Utility::string_to_enum<libMeshEnums::FEFamily>(d_feFamilyName);
  
  d_qruleTypeName = (params->d_db)->getStringWithDefault("QRULE_TYPE", "QGAUSS");
  
  d_qruleType = Utility::string_to_enum<libMeshEnums::QuadratureType>(d_qruleTypeName);
  
  // This is a boundary integral
  const unsigned int dimension = 2;
  
  d_feType.reset( new ::FEType(d_feTypeOrder, d_feFamily) );
  
  d_fe.reset( (::FEBase::build(dimension, (*d_feType))).release() );
  
  d_qruleOrderName = (params->d_db)->getStringWithDefault("QRULE_ORDER", "DEFAULT");
  
  if(d_qruleOrderName == "DEFAULT")
    {
      d_qruleOrder = d_feType->default_quadrature_order();
    }
  else
    {
      d_qruleOrder = Utility::string_to_enum<libMeshEnums::Order>(d_qruleOrderName);
    }
  
  d_qrule.reset( (::QBase::build(d_qruleType, dimension, d_qruleOrder)).release() );
  
  d_fe->attach_quadrature_rule( d_qrule.get() );
  
  d_JxW = &(d_fe->get_JxW());
  
  d_dphi = &(d_fe->get_dphi());
  
  d_variable = params->d_variable;
  
  d_NeumannParams.reset(new AMP::Operator::NeumannVectorCorrectionParameters( params->d_db ));
  d_NeumannParams->d_variable = params->d_variable;
  d_NeumannParams->d_MeshAdapter = params->d_MeshAdapter;
  d_NeumannParams->d_variableFlux = params->d_variableFlux;
  d_NeumannCorrection.reset(new NeumannVectorCorrection (d_NeumannParams));
  
  reset(params);
}
  
void RobinMatrixCorrection :: reset(const boost::shared_ptr<OperatorParameters>& params)
{
  boost::shared_ptr<RobinMatrixCorrectionParameters> myparams =
    boost::dynamic_pointer_cast<RobinMatrixCorrectionParameters>(params);
  
  AMP_INSIST( ((myparams.get()) != NULL), "NULL parameters" );
  AMP_INSIST( (((myparams->d_db).get()) != NULL), "NULL database" );
  
  AMP_INSIST( (((myparams->d_db).get()) != NULL), "NULL database" );
  bool skipParams = (myparams->d_db)->getBoolWithDefault("skip_params", true);
  
  if(!skipParams)
    {
      AMP_INSIST( (myparams->d_db)->keyExists("alpha"), "Missing key: prefactor alpha" );
      d_alpha    = (myparams->d_db)->getDouble("alpha");
      AMP_INSIST( d_alpha != 0.0, "prefactor alpha must be != 0.0" );
      
      AMP_INSIST( (myparams->d_db)->keyExists("beta"), "Missing key: prefactor beta" );
      d_beta.resize(1);
      d_beta[0]     = (myparams->d_db)->getDouble("beta");
      
      AMP_INSIST( (myparams->d_db)->keyExists("gamma"), "Missing key: total prefactor gamma" );
      d_gamma.resize(1);
      d_gamma[0]    = (myparams->d_db)->getDouble("gamma");
      
      AMP_INSIST( (myparams->d_db)->keyExists("fConductivity"), "Missing key: effective convective coefficient" );
      d_hef   = (myparams->d_db)->getDouble("fConductivity");
      
      AMP_INSIST( (params->d_db)->keyExists("number_of_ids"), "Key ''number_of_ids'' is missing!" );
      int numIds = (params->d_db)->getInteger("number_of_ids");
      
      d_boundaryIds.resize(numIds);
      d_dofIds.resize(numIds);
      d_robinValues.resize(numIds);
      
      char key[100];
      for(int j = 0; j < numIds; j++)
    {
      
      sprintf(key, "id_%d", j);
      AMP_INSIST( (myparams->d_db)->keyExists(key), "Key is missing!" );
      d_boundaryIds[j] = (myparams->d_db)->getInteger(key);
      
      sprintf(key, "number_of_dofs_%d", j);
      AMP_INSIST( (myparams->d_db)->keyExists(key), "Key is missing!" );
      int numDofIds = (myparams->d_db)->getInteger(key);
      
      d_dofIds[j].resize(numDofIds);
      d_robinValues[j].resize(numDofIds);
      for(int i = 0; i < numDofIds; i++)
        {
          sprintf(key, "dof_%d_%d", j, i);
          AMP_INSIST( (myparams->d_db)->keyExists(key), "Key is missing!" );
          d_dofIds[j][i] = (myparams->d_db)->getInteger(key);
          
          sprintf(key, "value_%d_%d", j, i);
          AMP_INSIST( (myparams->d_db)->keyExists(key), "Key is missing!" );
          d_robinValues[j][i] = (myparams->d_db)->getDouble(key);
        }
    }
    }
  
  if(myparams->d_db->isDatabase("RobinPhysicsModel"))
    {
      d_robinPhysicsModel = myparams->d_robinPhysicsModel;
    }
  
  (d_NeumannParams->d_db)->putBool("constant_flux",myparams->d_db->getBoolWithDefault("constant_flux",true));
  d_NeumannParams->d_variableFlux = myparams->d_variableFlux;
  d_NeumannParams->d_robinPhysicsModel = myparams->d_robinPhysicsModel ;
  (d_NeumannParams->d_db)->putDouble("gamma",d_gamma[0]);
  d_NeumannCorrection->reset(d_NeumannParams);
  
  bool skipMatrixCorrection = (myparams->d_db)->getBoolWithDefault("skip_matrix_correction", false);
  if(!skipMatrixCorrection)
    {
      AMP::LinearAlgebra::Matrix::shared_ptr inputMatrix = myparams->d_inputMatrix;
      AMP_INSIST( ((inputMatrix.get()) != NULL), "NULL matrix" );
      
      AMP::Mesh::DOFMap::shared_ptr dof_map = d_MeshAdapter->getDOFMap(d_variable);
      
      std::vector<std::string> variableNames;
      if(d_robinPhysicsModel.get() != NULL)
    {
      variableNames = d_robinPhysicsModel->getVariableName();
    }
      
      unsigned int numIds = d_boundaryIds.size();
      
      for(unsigned int nid = 0; nid < numIds; nid++)
    {
      AMP::Mesh::MeshManager::Adapter::BoundarySideIterator bnd1 = d_MeshAdapter->beginSideBoundary( d_boundaryIds[nid] );
      AMP::Mesh::MeshManager::Adapter::BoundarySideIterator end_bnd1 = d_MeshAdapter->endSideBoundary( d_boundaryIds[nid] );
      for( ; bnd1 != end_bnd1; ++bnd1)
        {
          
          AMP::Mesh::MeshManager::Adapter::Element cur_side = *bnd1;

          
          d_feType.reset( new ::FEType(d_feTypeOrder, d_feFamily) );
          d_fe.reset( (::FEBase::build(2, (*d_feType))).release() );
          
          d_phi = &(d_fe->get_phi());
          d_JxW = &(d_fe->get_JxW());
          
          d_qrule.reset( (::QBase::build(d_qruleType, 2, d_qruleOrder)).release() );
          d_fe->attach_quadrature_rule( d_qrule.get() );
          
          d_fe->reinit ( &cur_side.getElem() );
          
          const std::vector<Real> & JxW = (*d_JxW);
          const std::vector<std::vector<Real> > & phi = (*d_phi);
          
          std::vector<unsigned int> bndGlobalIds;
          dof_map->getDOFs(cur_side, bndGlobalIds);
          double temp;
          
          for(unsigned int qp = 0; qp < d_qrule->n_points(); qp++)
        {
          std::vector<std::vector<double> > inputArgs(1) ;
          for (unsigned int j=0; j<bndGlobalIds.size(); j++)
            {
              for (unsigned int i=0; i<bndGlobalIds.size(); i++)
            {
              temp =  d_beta[0] * ( JxW[qp]*phi[j][qp]*phi[i][qp] ) ;
              inputMatrix->addValueByGlobalID ( bndGlobalIds[j], bndGlobalIds[i], temp );
            }//end for i
            }//end for j
        }//end for qp
        }//end for bnd
    }// end for nid
      inputMatrix->makeConsistent();
    }//skip matrix
}
  
}
}

