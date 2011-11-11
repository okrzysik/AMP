
#include "NeumannVectorCorrection.h"
#include "NeumannVectorCorrectionParameters.h"
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


    // Constructor
    NeumannVectorCorrection::NeumannVectorCorrection(const boost::shared_ptr<NeumannVectorCorrectionParameters> & params)
        : BoundaryOperator (params)
      {
              d_params = params;

          std::string feTypeOrderName = (params->d_db)->getStringWithDefault("FE_ORDER", "FIRST");

          libMeshEnums::Order feTypeOrder = Utility::string_to_enum<libMeshEnums::Order>(feTypeOrderName);

          std::string feFamilyName = (params->d_db)->getStringWithDefault("FE_FAMILY", "LAGRANGE");

          libMeshEnums::FEFamily feFamily = Utility::string_to_enum<libMeshEnums::FEFamily>(feFamilyName);

          std::string qruleTypeName = (params->d_db)->getStringWithDefault("QRULE_TYPE", "QGAUSS");

          libMeshEnums::QuadratureType qruleType = Utility::string_to_enum<libMeshEnums::QuadratureType>(qruleTypeName);

          const unsigned int dimension = 2;

          d_feType.reset( new ::FEType(feTypeOrder, feFamily) );

          d_fe.reset( (::FEBase::build(dimension, (*d_feType))).release() );

          std::string qruleOrderName = (params->d_db)->getStringWithDefault("QRULE_ORDER", "DEFAULT");

          libMeshEnums::Order qruleOrder;

          if(qruleOrderName == "DEFAULT") {
              qruleOrder = d_feType->default_quadrature_order();
          } else {
              qruleOrder = Utility::string_to_enum<libMeshEnums::Order>(qruleOrderName);
          }

          d_qrule.reset( (::QBase::build(qruleType, dimension, qruleOrder)).release() );

          d_fe->attach_quadrature_rule( d_qrule.get() );

          d_JxW = &(d_fe->get_JxW());

          d_variable = params->d_variable;

          reset(params);
      }


  void NeumannVectorCorrection :: reset(const boost::shared_ptr<OperatorParameters>& params)
  {
    boost::shared_ptr<NeumannVectorCorrectionParameters> myparams = 
      boost::dynamic_pointer_cast<NeumannVectorCorrectionParameters>(params);

    AMP_INSIST( ((myparams.get()) != NULL), "NULL parameters" );
    AMP_INSIST( (((myparams->d_db).get()) != NULL), "NULL database" );

    AMP_INSIST( (myparams->d_db)->keyExists("number_of_ids"), "Key ''number_of_ids'' is missing!" );
    d_numIds = (myparams->d_db)->getInteger("number_of_ids");

    d_isConstantFlux = (myparams->d_db)->getBoolWithDefault("constant_flux", true);

    d_boundaryIds.resize(d_numIds);
    d_dofIds.resize(d_numIds);
    d_neumannValues.resize(d_numIds);
    d_IsCoupledBoundary.resize(d_numIds);
    d_numDofIds.resize(d_numIds);

    char key[100];
    for(int j = 0; j < d_numIds; j++) {
      sprintf(key, "id_%d", j);
      AMP_INSIST( (myparams->d_db)->keyExists(key), "Key is missing!" );
      d_boundaryIds[j] = (myparams->d_db)->getInteger(key);

      sprintf(key, "number_of_dofs_%d", j);
      AMP_INSIST( (myparams->d_db)->keyExists(key), "Key is missing!" );
      d_numDofIds[j] = (myparams->d_db)->getInteger(key);

      sprintf(key, "IsCoupledBoundary_%d", j);
      d_IsCoupledBoundary[j] = (params->d_db)->getBoolWithDefault(key, false);

      d_dofIds[j].resize(d_numDofIds[j]);
      d_neumannValues[j].resize(d_numDofIds[j]);
      for(int i = 0; i < d_numDofIds[j]; i++) {
        sprintf(key, "dof_%d_%d", j, i);
        AMP_INSIST( (myparams->d_db)->keyExists(key), "Key is missing!" );
        d_dofIds[j][i] = (myparams->d_db)->getInteger(key);

        if(d_isConstantFlux){
          sprintf(key, "value_%d_%d", j, i);
          AMP_INSIST( (myparams->d_db)->keyExists(key), "Key is missing!" );
          d_neumannValues[j][i] = (myparams->d_db)->getDouble(key);
        }else{
          d_variableFlux = myparams->d_variableFlux;      
        }
      }//end for i
    }//end for j

#if 0    
    if(myparams->d_db->isDatabase("RobinPhysicsModel") ){
#else
    if(myparams->d_robinPhysicsModel) {
#endif
      
      d_robinPhysicsModel = myparams->d_robinPhysicsModel;
    }
    computeRHScorrection(myparams);
  }

void
NeumannVectorCorrection :: computeRHScorrection(const boost::shared_ptr<NeumannVectorCorrectionParameters> & params)
{
AMP_ERROR("NeumannVectorCorrection.cc Not fixed yet");
/*
  d_gamma.resize(1);
  d_gamma[0] = 1.0;
  
  if(!d_isConstantFlux)
    {
      d_variableFlux->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
    }
  
  if((params->d_db)->keyExists("gamma"))
    {
      d_gamma[0] = (params->d_db)->getDouble("gamma");
    }
  
  AMP::Mesh::DOFMap::shared_ptr dof_map = d_MeshAdapter->getDOFMap(d_variable);
  
  AMP::LinearAlgebra::Vector::shared_ptr rInternal = d_MeshAdapter->createVector(d_variable);
  rInternal->zero();
  
  unsigned int numIds = d_boundaryIds.size();
  
  for(unsigned int j = 0; j < numIds; j++)
    {
      if(!d_IsCoupledBoundary[j])
    {
      unsigned int numDofIds = d_dofIds[j].size();
      
      for(unsigned int k = 0; k < numDofIds; k++)
        {
          
          unsigned int dofId;
          dofId = d_dofIds[j][k];
          
          //This could be changed to BoundaryOwnNodeIterator
          AMP::Mesh::MeshManager::Adapter::BoundarySideIterator bnd = d_MeshAdapter->beginSideBoundary( d_boundaryIds[j] );
          AMP::Mesh::MeshManager::Adapter::BoundarySideIterator end_bnd = d_MeshAdapter->endSideBoundary( d_boundaryIds[j] );
          
          int count =0;
          for( ; bnd != end_bnd; ++bnd)
        {
          count++;
          
          AMP::Mesh::MeshManager::Adapter::Element cur_side = *bnd;
          
          d_phi = &(d_fe->get_phi());
          d_JxW = &(d_fe->get_JxW());
          
          d_fe->reinit ( &cur_side.getElem() );
          
          std::vector<unsigned int> bndGlobalIds;
          dof_map->getDOFs(cur_side, bndGlobalIds, dofId);
          
          std::vector<double> flux(bndGlobalIds.size(), 0.0);
          const std::vector<std::vector<Real> > &phi = *d_phi;
          const std::vector<Real> &djxw = *d_JxW;
          
          for(unsigned int i = 0; i < bndGlobalIds.size(); i++)
            {
              for(unsigned int qp = 0; qp < d_qrule->n_points(); qp++) 
            {
              std::vector<std::vector<double> > temp(1) ;
              // fixing so it will work in both constant and non constant case
              // also the if is based on a lousy conditional
#if 0        
              if(d_isConstantFlux)
                {
                  flux[i] +=  (d_gamma[0])*djxw[qp]*phi[i][qp]*d_neumannValues[j][k];
                }
              else
                {
                  temp[0].push_back(d_variableFlux->getValueByGlobalID(bndGlobalIds[i]) );
                  if(params->d_db->isDatabase("RobinPhysicsModel"))
                {
                  d_robinPhysicsModel->getConductance(d_gamma, d_gamma, temp); 
                }
                  flux[i] +=  (d_gamma[0])*djxw[qp]*phi[i][qp]*temp[0][0];
                }
#else
              // there must be a better way to write this!!
              if(d_isConstantFlux)
                {
                  temp[0].push_back(d_neumannValues[j][k]);            
                }
              else
                {
                  temp[0].push_back(d_variableFlux->getValueByGlobalID(bndGlobalIds[i]) );
                }
              
              if(d_robinPhysicsModel)
                {
                  d_robinPhysicsModel->getConductance(d_gamma, d_gamma, temp); 
                }
              
              flux[i] +=  (d_gamma[0])*djxw[qp]*phi[i][qp]*temp[0][0];
#endif
              
            }
              
            }//end for i
          rInternal->addValuesByGlobalID((int)bndGlobalIds.size() , (int *)&(bndGlobalIds[0]), &(flux[0]));
          
        }//end for bnd
          
        }//end for k
    }//coupled
    }//end for j
  ////////////////////////////////////////////////////////////////
  
  rInternal->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_ADD );
  
  d_rhsCorrectionAdd = rInternal->cloneVector();
  d_rhsCorrectionAdd->copyVector(rInternal);
*/
}
 
 void NeumannVectorCorrection :: apply(const AMP::LinearAlgebra::Vector::shared_ptr &f, const AMP::LinearAlgebra::Vector::shared_ptr &u, AMP::LinearAlgebra::Vector::shared_ptr  &r, const double a , const double b )
  {
    (void) f; (void) u; (void) r; (void) a; (void) b; 
    //Do Nothing
  }

  boost::shared_ptr<OperatorParameters> NeumannVectorCorrection :: 
    getJacobianParameters(const boost::shared_ptr<AMP::LinearAlgebra::Vector>& ) {
      boost::shared_ptr<AMP::InputDatabase> tmp_db (new AMP::InputDatabase("Dummy"));

      tmp_db->putString("FE_ORDER","FIRST");
      tmp_db->putString("FE_FAMILY","LAGRANGE");
      tmp_db->putString("QRULE_TYPE","QGAUSS");
      tmp_db->putInteger("DIMENSION",2);
      tmp_db->putString("QRULE_ORDER","DEFAULT");
      tmp_db->putInteger("number_of_ids",d_numIds);
      tmp_db->putBool("constant_flux", d_isConstantFlux);

      char key[100];
      for(int j = 0; j < d_numIds; j++) {
        sprintf(key, "id_%d", j);
        tmp_db->putInteger(key,d_boundaryIds[j]);
        sprintf(key, "number_of_dofs_%d", j);
        tmp_db->putInteger(key,d_numDofIds[j]);
        sprintf(key, "IsCoupledBoundary_%d", j);
        tmp_db->putBool(key, d_IsCoupledBoundary[j]);

        for(int i = 0; i < d_numDofIds[j]; i++) {
          sprintf(key, "dof_%d_%d", j, i);
          tmp_db->putInteger(key,d_dofIds[j][i]);
          if(d_isConstantFlux){
            sprintf(key, "value_%d_%d", j, i);
            tmp_db->putInteger(key,d_neumannValues[j][i]);
          }else{
            //d_variableFlux ??
          }
        }
      }

      tmp_db->putBool("skip_params", true);

      boost::shared_ptr<NeumannVectorCorrectionParameters> outParams(new NeumannVectorCorrectionParameters(tmp_db));

      return outParams;
    }


  void setFrozenVector ( AMP::LinearAlgebra::Vector::shared_ptr f ) {
AMP_ERROR("NeumannVectorCorrection.cc Not fixed yet");
/*
        if ( d_Frozen )
        {
          d_Frozen->castTo<AMP::LinearAlgebra::MultiVector>().addVector ( f );
        }
        else
        {
          d_Frozen = AMP::LinearAlgebra::MultiVector::view ( f );
        }
        d_Frozen = d_Frozen->select ( AMP::Mesh::VS_ByMesh ( d_MeshAdapter ) , d_Frozen->getVariable()->getName() );
*/
  }



}
}

