
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
              AMP_INSIST( (params->d_db)->keyExists(key), "Key is missing!" );
              d_dirichletValues1[j][i] = (params->d_db)->getDouble(key);
            }
          }//end for i
        }//end for j
      }
    }

    //This is an in-place apply
    void DirichletVectorCorrection :: apply(const AMP::LinearAlgebra::Vector::shared_ptr &, const AMP::LinearAlgebra::Vector::shared_ptr &u,
        AMP::LinearAlgebra::Vector::shared_ptr  &r, const double a, const double ) {

      AMP::LinearAlgebra::Vector::shared_ptr rInternal = r->subsetVectorForVariable(d_variable);

      AMP_ASSERT(rInternal != NULL);

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

      if(d_iDebugPrintInfoLevel>3)
      {
        AMP::pout << "L2 Norm of rInternal leaving DirichletVectorCorrection::apply is : " << rInternal->L2Norm() << std::endl;
      }
    }

    void DirichletVectorCorrection :: applyZeroValues(AMP::LinearAlgebra::Vector::shared_ptr r) {
      AMP_ASSERT(d_MeshAdapter != NULL);

      AMP::Mesh::DOFMap::shared_ptr dof_map = d_MeshAdapter->getDOFMap(d_variable);

      AMP_ASSERT(dof_map != NULL);

      AMP::LinearAlgebra::Vector::shared_ptr rInternal = r->subsetVectorForVariable(d_variable);

      AMP_ASSERT(rInternal != NULL);

      unsigned int numIds = d_boundaryIds.size();

      for(unsigned int j = 0; j < numIds; j++) {
        AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator bnd = d_MeshAdapter->beginOwnedBoundary( d_boundaryIds[j] );
        AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator end_bnd = d_MeshAdapter->endOwnedBoundary( d_boundaryIds[j] );

        for( ; bnd != end_bnd; ++bnd) {
          std::vector<unsigned int> bndGlobalIds;
          dof_map->getDOFs(*bnd, bndGlobalIds, d_dofIds[j]);

          for(unsigned int i = 0; i < bndGlobalIds.size(); i++) {
            rInternal->setLocalValueByGlobalID(bndGlobalIds[i], 0.0);
          }//end for i
        }//end for bnd
      }//end for j

    }

    void DirichletVectorCorrection :: applyNonZeroValues(AMP::LinearAlgebra::Vector::shared_ptr r) {
      AMP_ASSERT(d_MeshAdapter != NULL);

      AMP::Mesh::DOFMap::shared_ptr dof_map = d_MeshAdapter->getDOFMap(d_variable);

      AMP_ASSERT(dof_map != NULL);

      AMP::LinearAlgebra::Vector::shared_ptr rInternal = r->subsetVectorForVariable(d_variable);

      AMP_ASSERT(rInternal != NULL);

      if(d_valuesType == 2) {
        AMP_ASSERT(d_dirichletValues2 != NULL);
      }

      unsigned int numIds = d_boundaryIds.size();

      for(unsigned int j = 0; j < numIds; j++) {
        AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator bnd = d_MeshAdapter->beginOwnedBoundary( d_boundaryIds[j] );
        AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator end_bnd = d_MeshAdapter->endOwnedBoundary( d_boundaryIds[j] );

        for( ; bnd != end_bnd; ++bnd) {
          std::vector<unsigned int> bndGlobalIds;
          dof_map->getDOFs(*bnd, bndGlobalIds, d_dofIds[j]);

          for(unsigned int i = 0; i < bndGlobalIds.size(); i++) {
            double dVal;
            if(d_valuesType == 1) {
              dVal = d_dirichletValues1[j][i];
            } else {
              dVal = d_dirichletValues2->getLocalValueByGlobalID( bndGlobalIds[i] );
            }
            rInternal->setLocalValueByGlobalID(bndGlobalIds[i], d_scalingFactor*dVal);
          }//end for i
        }//end for bnd
      }//end for j

    }

    void DirichletVectorCorrection :: applyResidual(AMP::LinearAlgebra::Vector::shared_ptr u, AMP::LinearAlgebra::Vector::shared_ptr r) {
      AMP_ASSERT(d_MeshAdapter != NULL);

      AMP::Mesh::DOFMap::shared_ptr dof_map = d_MeshAdapter->getDOFMap(d_variable);

      AMP_ASSERT(dof_map != NULL);

      AMP::LinearAlgebra::Vector::shared_ptr uInternal = u->subsetVectorForVariable(d_variable);

      AMP_ASSERT(uInternal != NULL);

      unsigned int numIds = d_boundaryIds.size();

      for(unsigned int j = 0; j < numIds; j++) {
        AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator bnd = d_MeshAdapter->beginOwnedBoundary( d_boundaryIds[j] );
        AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator end_bnd = d_MeshAdapter->endOwnedBoundary( d_boundaryIds[j] );

        for( ; bnd != end_bnd; ++bnd) {
          std::vector<unsigned int> bndGlobalIds;
          dof_map->getDOFs(*bnd, bndGlobalIds, d_dofIds[j]);

          for(unsigned int i = 0; i < bndGlobalIds.size(); i++) {
            double uVal = uInternal->getLocalValueByGlobalID( bndGlobalIds[i] );
            double dVal;
            if(d_valuesType == 1) {
              dVal = d_dirichletValues1[j][i];
            } else {
              dVal = d_dirichletValues2->getLocalValueByGlobalID( bndGlobalIds[i] );
            }
            r->setLocalValueByGlobalID(bndGlobalIds[i], d_scalingFactor*(uVal - dVal));
          }//end for i
        }//end for bnd
      }//end for j

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

