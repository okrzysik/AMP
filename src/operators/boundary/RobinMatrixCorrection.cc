
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

#include "face_quad4.h"
#include "node.h"

#include <string>

namespace AMP {
namespace Operator {

RobinMatrixCorrection :: RobinMatrixCorrection(const boost::shared_ptr<RobinMatrixCorrectionParameters> & params)
  : BoundaryOperator (params)
{
  std::string feTypeOrderName = (params->d_db)->getStringWithDefault("FE_ORDER", "FIRST");
  std::string feFamilyName = (params->d_db)->getStringWithDefault("FE_FAMILY", "LAGRANGE");
  std::string qruleTypeName = (params->d_db)->getStringWithDefault("QRULE_TYPE", "QGAUSS");
  d_qruleOrderName = (params->d_db)->getStringWithDefault("QRULE_ORDER", "DEFAULT");
  
  d_feTypeOrder = Utility::string_to_enum<libMeshEnums::Order>(feTypeOrderName);
  d_feFamily = Utility::string_to_enum<libMeshEnums::FEFamily>(feFamilyName);
  d_qruleType = Utility::string_to_enum<libMeshEnums::QuadratureType>(qruleTypeName);
  
  d_variable = params->d_variable;
  
  d_dofManager = (params->d_DofMap);

  d_NeumannParams.reset(new AMP::Operator::NeumannVectorCorrectionParameters( params->d_db ));
  d_NeumannParams->d_variable = params->d_variable;
  d_NeumannParams->d_Mesh = params->d_Mesh;
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

    std::vector<std::string> variableNames;

    if(d_robinPhysicsModel.get() != NULL)
    {
      variableNames = d_robinPhysicsModel->getVariableName();
    }

    unsigned int numIds = d_boundaryIds.size();

    for(unsigned int nid = 0; nid < numIds; nid++)
    {
      AMP::Mesh::MeshIterator bnd1     = d_Mesh->getIDsetIterator( AMP::Mesh::Face, d_boundaryIds[nid], 0 );
      AMP::Mesh::MeshIterator end_bnd1 = bnd1.end();
      for( ; bnd1 != end_bnd1; ++bnd1)
      {

        boost::shared_ptr < ::FEType > d_feType ( new ::FEType(d_feTypeOrder, d_feFamily) );
        boost::shared_ptr < ::FEBase > d_fe( (::FEBase::build(2, (*d_feType))).release() );

        if(d_qruleOrderName == "DEFAULT") {
          d_qruleOrder = d_feType->default_quadrature_order();
        } else {
          d_qruleOrder = Utility::string_to_enum<libMeshEnums::Order>(d_qruleOrderName);
        }
        boost::shared_ptr < ::QBase > d_qrule( (::QBase::build(d_qruleType, 2, d_qruleOrder)).release() );

        d_currNodes = bnd1->getElements(AMP::Mesh::Vertex);

        unsigned int numNodesInCurrElem = d_currNodes.size();

        createCurrentLibMeshElement();

        getDofIndicesForCurrentElement();

        d_fe->attach_quadrature_rule( d_qrule.get() );

        d_phi = &(d_fe->get_phi());
        d_JxW = &(d_fe->get_JxW());

        d_fe->reinit ( d_currElemPtr );

        const std::vector<Real> & JxW = (*d_JxW);
        const std::vector<std::vector<Real> > & phi = (*d_phi);

        double temp;

        for(unsigned int qp = 0; qp < d_qrule->n_points(); qp++)
        {
          std::vector<std::vector<double> > inputArgs(1) ;
          for (unsigned int j=0; j < numNodesInCurrElem ; j++)
          {
            for (unsigned int i=0; i < numNodesInCurrElem ; i++)
            {
              temp =  d_beta[0] * ( JxW[qp]*phi[j][qp]*phi[i][qp] ) ;
              inputMatrix->addValueByGlobalID ( d_dofIndices[j], d_dofIndices[i], temp );
            }//end for i
          }//end for j
        }//end for qp
        
        destroyCurrentLibMeshElement();

      }//end for bnd

    }// end for nid
    
    inputMatrix->makeConsistent();

  }//skip matrix
  
}

void RobinMatrixCorrection :: createCurrentLibMeshElement() {
  d_currElemPtr = new ::Quad4;
  for(size_t j = 0; j < d_currNodes.size(); j++) {
    std::vector<double> pt = d_currNodes[j].coord();
    d_currElemPtr->set_node(j) = new ::Node(pt[0], pt[1], pt[2], j);
  }//end for j
}

void RobinMatrixCorrection :: destroyCurrentLibMeshElement() {
  for(size_t j = 0; j < d_currElemPtr->n_nodes(); j++) {
    delete (d_currElemPtr->get_node(j));
    d_currElemPtr->set_node(j) = NULL;
  }//end for j
  delete d_currElemPtr;
  d_currElemPtr = NULL;
}

void RobinMatrixCorrection :: getDofIndicesForCurrentElement() {
  d_dofIndices.resize(d_currNodes.size());
  std::vector<AMP::Mesh::MeshElementID> globalIDs(d_currNodes.size()); 
  for(unsigned int j = 0; j < d_currNodes.size(); j++) {
    globalIDs[j] = d_currNodes[j].globalID();
  } // end of j
  d_dofManager->getDOFs(globalIDs, d_dofIndices);
}

}
}

