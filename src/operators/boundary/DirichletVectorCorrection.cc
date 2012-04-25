
#include "DirichletVectorCorrection.h"
#include "DirichletMatrixCorrectionParameters.h"
#include "utils/Utilities.h"
#include "utils/InputDatabase.h"

namespace AMP {
  namespace Operator {

    void DirichletVectorCorrection :: reset(const boost::shared_ptr<OperatorParameters> & tmpParams)
    {
      boost::shared_ptr<DirichletVectorCorrectionParameters> params = 
        boost::dynamic_pointer_cast<DirichletVectorCorrectionParameters>(tmpParams);

      AMP_INSIST( ((params.get()) != NULL), "NULL parameters" );
      AMP_INSIST( (((params->d_db).get()) != NULL), "NULL database" );

      bool skipParams = (params->d_db)->getBoolWithDefault("skip_params", false);

      if(!skipParams) {
        d_scalingFactor = (params->d_db)->getDoubleWithDefault("SCALING_FACTOR", 1.0);
        d_setResidual = (params->d_db)->getBoolWithDefault("setResidual", false);
        d_isAttachedToVolumeOperator = (params->d_db)->getBoolWithDefault("isAttachedToVolumeOperator", false);
        d_valuesType = (params->d_db)->getIntegerWithDefault("valuesType", 1);
        AMP_INSIST( ((d_valuesType == 1) || (d_valuesType == 2)), "Wrong value.");

        AMP_INSIST( (params->d_db)->keyExists("number_of_ids"), "Key ''number_of_ids'' is missing!" );
        int numIds = (params->d_db)->getInteger("number_of_ids");

        d_boundaryIds.resize(numIds);
        d_dofIds.resize(numIds);

        if(d_valuesType == 1) {
          d_dirichletValues1.resize(numIds);
        }

        char key[100];
        for(int j = 0; j < numIds; j++) {
          sprintf(key, "id_%d", j);
          AMP_INSIST( (params->d_db)->keyExists(key), "Key is missing!" );
          d_boundaryIds[j] = (params->d_db)->getInteger(key);

          sprintf(key, "number_of_dofs_%d", j);
          AMP_INSIST( (params->d_db)->keyExists(key), "Key is missing!" );
          int numDofIds = (params->d_db)->getInteger(key);

          d_dofIds[j].resize(numDofIds);

          if(d_valuesType == 1) {
            d_dirichletValues1[j].resize(numDofIds);
          }
          for(int i = 0; i < numDofIds; i++) {
            sprintf(key, "dof_%d_%d", j, i);
            AMP_INSIST( (params->d_db)->keyExists(key), "Key is missing!" );
            d_dofIds[j][i] = (params->d_db)->getInteger(key);

            if(d_valuesType == 1) {
              sprintf(key, "value_%d_%d", j, i);
              d_dirichletValues1[j][i] = (params->d_db)->getDoubleWithDefault(key, 0.0);
            }
          }//end for i
        }//end for j
      }
    }

    //This is an in-place apply
    void DirichletVectorCorrection :: apply(const AMP::LinearAlgebra::Vector::shared_ptr &, const AMP::LinearAlgebra::Vector::shared_ptr &u,
        AMP::LinearAlgebra::Vector::shared_ptr  &r, const double a, const double ) {
      AMP::LinearAlgebra::Vector::shared_ptr rInternal = mySubsetVector(r, d_variable);

      AMP_ASSERT( ((*(r->getUpdateStatus())) == AMP::LinearAlgebra::Vector::UNCHANGED) || 
          ((*(r->getUpdateStatus())) == AMP::LinearAlgebra::Vector::LOCAL_CHANGED) );
      AMP_ASSERT( ((*(rInternal->getUpdateStatus())) == AMP::LinearAlgebra::Vector::UNCHANGED) || 
          ((*(rInternal->getUpdateStatus())) == AMP::LinearAlgebra::Vector::LOCAL_CHANGED) );

      if(d_iDebugPrintInfoLevel>3)
      {
        AMP::pout << "L2 Norm of rInternal entering DirichletVectorCorrection::apply is : " << rInternal->L2Norm() << std::endl;
      }

      if(d_setResidual) {
        this->applyResidual(u, rInternal);
      } else if(d_isAttachedToVolumeOperator) {
        this->applyZeroValues(rInternal);
      } else {
        this->applyNonZeroValues(rInternal);
      }

      rInternal->scale(a);

      AMP_ASSERT( ((*(r->getUpdateStatus())) == AMP::LinearAlgebra::Vector::UNCHANGED) || 
          ((*(r->getUpdateStatus())) == AMP::LinearAlgebra::Vector::LOCAL_CHANGED) );
      AMP_ASSERT( ((*(rInternal->getUpdateStatus())) == AMP::LinearAlgebra::Vector::UNCHANGED) || 
          ((*(rInternal->getUpdateStatus())) == AMP::LinearAlgebra::Vector::LOCAL_CHANGED) );

      if(d_iDebugPrintInfoLevel>3)
      {
        AMP::pout << "L2 Norm of rInternal leaving DirichletVectorCorrection::apply is : " << rInternal->L2Norm() << std::endl;
      }
    }

    void DirichletVectorCorrection :: applyZeroValues(AMP::LinearAlgebra::Vector::shared_ptr r) {
      AMP::LinearAlgebra::Vector::shared_ptr rInternal = mySubsetVector(r, d_variable);
      AMP::Discretization::DOFManager::shared_ptr dof_map = rInternal->getDOFManager();
      size_t numIds = d_boundaryIds.size();
      AMP_ASSERT( ((*(r->getUpdateStatus())) == AMP::LinearAlgebra::Vector::UNCHANGED) || 
          ((*(r->getUpdateStatus())) == AMP::LinearAlgebra::Vector::LOCAL_CHANGED) );
      AMP_ASSERT( ((*(rInternal->getUpdateStatus())) == AMP::LinearAlgebra::Vector::UNCHANGED) || 
          ((*(rInternal->getUpdateStatus())) == AMP::LinearAlgebra::Vector::LOCAL_CHANGED) );
      for(size_t j = 0; j < numIds; j++) {
        AMP::Mesh::MeshIterator bnd = d_Mesh->getIDsetIterator( AMP::Mesh::Vertex, d_boundaryIds[j], 0 );
        AMP::Mesh::MeshIterator end_bnd = bnd.end();

        for( ; bnd != end_bnd; ++bnd) {
          std::vector<size_t> bndGlobalIds;
          dof_map->getDOFs(bnd->globalID(), bndGlobalIds);
          for(size_t i = 0; i < d_dofIds[j].size(); i++) {
            rInternal->setLocalValueByGlobalID(bndGlobalIds[d_dofIds[j][i]], 0.0);
          }//end for i
        }//end for bnd
      }//end for j
      AMP_ASSERT( ((*(r->getUpdateStatus())) == AMP::LinearAlgebra::Vector::UNCHANGED) || 
          ((*(r->getUpdateStatus())) == AMP::LinearAlgebra::Vector::LOCAL_CHANGED) );
      AMP_ASSERT( ((*(rInternal->getUpdateStatus())) == AMP::LinearAlgebra::Vector::UNCHANGED) || 
          ((*(rInternal->getUpdateStatus())) == AMP::LinearAlgebra::Vector::LOCAL_CHANGED) );
    }

    void DirichletVectorCorrection :: applyNonZeroValues(AMP::LinearAlgebra::Vector::shared_ptr r) {
      AMP::LinearAlgebra::Vector::shared_ptr rInternal = mySubsetVector(r, d_variable);
      AMP::Discretization::DOFManager::shared_ptr dof_map = rInternal->getDOFManager();
      size_t numIds = d_boundaryIds.size();
      AMP_ASSERT( ((*(r->getUpdateStatus())) == AMP::LinearAlgebra::Vector::UNCHANGED) || 
          ((*(r->getUpdateStatus())) == AMP::LinearAlgebra::Vector::LOCAL_CHANGED) );
      AMP_ASSERT( ((*(rInternal->getUpdateStatus())) == AMP::LinearAlgebra::Vector::UNCHANGED) || 
          ((*(rInternal->getUpdateStatus())) == AMP::LinearAlgebra::Vector::LOCAL_CHANGED) );
      for(size_t j = 0; j < numIds; j++) {
        AMP::Mesh::MeshIterator bnd = d_Mesh->getIDsetIterator( AMP::Mesh::Vertex, d_boundaryIds[j], 0 );
        AMP::Mesh::MeshIterator end_bnd = bnd.end();

        for( ; bnd != end_bnd; ++bnd) {
          std::vector<size_t> bndGlobalIds;
          dof_map->getDOFs(bnd->globalID(), bndGlobalIds);

          for(size_t i = 0; i < d_dofIds[j].size(); ++i) {
            double dVal;
            if(d_valuesType == 1) {
              dVal = d_dirichletValues1[j][i];
            } else {
              dVal = d_dirichletValues2->getLocalValueByGlobalID( bndGlobalIds[d_dofIds[j][i]] );
            }
            rInternal->setLocalValueByGlobalID(bndGlobalIds[d_dofIds[j][i]], d_scalingFactor*dVal);
          }//end for i
        }//end for bnd
      }//end for j
      AMP_ASSERT( ((*(r->getUpdateStatus())) == AMP::LinearAlgebra::Vector::UNCHANGED) || 
          ((*(r->getUpdateStatus())) == AMP::LinearAlgebra::Vector::LOCAL_CHANGED) );
      AMP_ASSERT( ((*(rInternal->getUpdateStatus())) == AMP::LinearAlgebra::Vector::UNCHANGED) || 
          ((*(rInternal->getUpdateStatus())) == AMP::LinearAlgebra::Vector::LOCAL_CHANGED) );
    }

    void DirichletVectorCorrection :: applyResidual(AMP::LinearAlgebra::Vector::shared_ptr u, 
        AMP::LinearAlgebra::Vector::shared_ptr r) {
      AMP::LinearAlgebra::Vector::shared_ptr uInternal = mySubsetVector(u, d_variable);
      AMP::Discretization::DOFManager::shared_ptr dof_map = uInternal->getDOFManager();
      size_t numIds = d_boundaryIds.size();
      AMP_ASSERT( ((*(r->getUpdateStatus())) == AMP::LinearAlgebra::Vector::UNCHANGED) || 
          ((*(r->getUpdateStatus())) == AMP::LinearAlgebra::Vector::LOCAL_CHANGED) );
      for(size_t j = 0; j < numIds; j++) {
        AMP::Mesh::MeshIterator bnd = d_Mesh->getIDsetIterator( AMP::Mesh::Vertex, d_boundaryIds[j], 0 );
        AMP::Mesh::MeshIterator end_bnd = bnd.end();

        for( ; bnd != end_bnd; ++bnd) {
          std::vector<size_t> bndGlobalIds;
          dof_map->getDOFs(bnd->globalID(), bndGlobalIds);
          for(size_t i = 0; i < d_dofIds[j].size(); i++) {
            double uVal = uInternal->getLocalValueByGlobalID( bndGlobalIds[d_dofIds[j][i]] );
            double dVal;
            if(d_valuesType == 1) {
              dVal = d_dirichletValues1[j][i];
            } else {
              dVal = d_dirichletValues2->getLocalValueByGlobalID( bndGlobalIds[d_dofIds[j][i]] );
            }
            r->setLocalValueByGlobalID(bndGlobalIds[d_dofIds[j][i]], d_scalingFactor*(uVal - dVal));
          }//end for i
        }//end for bnd
      }//end for j
      AMP_ASSERT( ((*(r->getUpdateStatus())) == AMP::LinearAlgebra::Vector::UNCHANGED) || 
          ((*(r->getUpdateStatus())) == AMP::LinearAlgebra::Vector::LOCAL_CHANGED) );
    }

    boost::shared_ptr<OperatorParameters> DirichletVectorCorrection :: 
      getJacobianParameters(const boost::shared_ptr<AMP::LinearAlgebra::Vector>& ) {
        boost::shared_ptr<AMP::InputDatabase> tmp_db (new AMP::InputDatabase("Dummy"));
        tmp_db->putBool("skip_params", true);

        boost::shared_ptr<DirichletMatrixCorrectionParameters> outParams(new DirichletMatrixCorrectionParameters(tmp_db));

        return outParams;
      }

  }
}



