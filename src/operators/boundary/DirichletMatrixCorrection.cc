
#include "DirichletMatrixCorrection.h"
#include "utils/Utilities.h"
#include "utils/InputDatabase.h"

namespace AMP {
namespace Operator {

  void DirichletMatrixCorrection :: reset(const boost::shared_ptr<OperatorParameters>& params)
  {
AMP_ERROR("DirichletMatrixCorrection.cc Not fixed yet");
/*
    boost::shared_ptr<DirichletMatrixCorrectionParameters> myParams = 
      boost::dynamic_pointer_cast<DirichletMatrixCorrectionParameters>(params);

    AMP_INSIST( ((myParams.get()) != NULL), "NULL parameters" );

    parseParams(myParams);

    computeRHScorrection(myParams);

    AMP::LinearAlgebra::Matrix::shared_ptr inputMatrix = myParams->d_inputMatrix;
    AMP_INSIST( ((inputMatrix.get()) != NULL), "NULL matrix" );

    AMP::Mesh::DOFMap::shared_ptr dof_map = d_MeshAdapter->getDOFMap ( d_variable );

    unsigned int numIds = d_boundaryIds.size();

    for(unsigned int k = 0; k < numIds; k++) {
      //This has to be a BoundarySharedNodeIterator and can not be a BoundaryOwnNodeIterator since
      //1) we don't have access to "SharedElements".
      //2) even we did have access to "SharedElements" and did a
      //BoundaryOwnNodeIteration we can't avoid the makeConsistent since the dof
      //of a sharedElement may be neither shared nor owned.
      AMP::Mesh::MeshManager::Adapter::BoundaryNodeIterator bnd = d_MeshAdapter->beginBoundary( d_boundaryIds[k] );
      AMP::Mesh::MeshManager::Adapter::BoundaryNodeIterator end_bnd = d_MeshAdapter->endBoundary( d_boundaryIds[k] );

      for( ; bnd != end_bnd; ++bnd) {
        AMP::Mesh::MeshManager::Adapter::NodeElementIterator el =  d_MeshAdapter->beginElementForNode ( *bnd );
        AMP::Mesh::MeshManager::Adapter::NodeElementIterator end_el = d_MeshAdapter->endElementForNode ( *bnd );

        std::vector<unsigned int> bndGlobalIds;
        dof_map->getDOFs(*bnd, bndGlobalIds, d_dofIds[k]);

        for( ; el != end_el; ++el) {
          std::vector<unsigned int> dofIndices;
          dof_map->getDOFs(*el, dofIndices);

          for(unsigned int j = 0; j < bndGlobalIds.size(); j++) {
            for(unsigned int i = 0; i < dofIndices.size(); i++) {
              if(bndGlobalIds[j] == dofIndices[i]) {
                if(d_zeroDirichletBlock) {
                  inputMatrix->setValueByGlobalID ( bndGlobalIds[j], bndGlobalIds[j], 0.0 );
                } else {
                  inputMatrix->setValueByGlobalID ( bndGlobalIds[j], bndGlobalIds[j], 1.0 );
                }
              } else {
                inputMatrix->setValueByGlobalID ( bndGlobalIds[j], dofIndices[i] , 0.0 );
                if(d_symmetricCorrection) {
                  inputMatrix->setValueByGlobalID ( dofIndices[i], bndGlobalIds[j] , 0.0 );
                }
              }
            }//end for i
          }//end for j
        }//end for el
      }//end for bnd
    }//end for k

    //This does consistent for both "Sum-into" and "set".
    inputMatrix->makeConsistent();
*/
  }

  void DirichletMatrixCorrection :: parseParams(const boost::shared_ptr<DirichletMatrixCorrectionParameters> & params)
  {
    AMP_INSIST( (((params->d_db).get()) != NULL), "NULL database" );
    bool skipParams = (params->d_db)->getBoolWithDefault("skip_params", false);

    if(!skipParams) {
      d_symmetricCorrection = (params->d_db)->getBoolWithDefault("symmetric_correction", true);
      d_zeroDirichletBlock = (params->d_db)->getBoolWithDefault("zero_dirichlet_block", false);

      d_skipRHSsetCorrection = (params->d_db)->getBoolWithDefault("skip_rhs_correction", true);
      d_skipRHSaddCorrection = (params->d_db)->getBoolWithDefault("skip_rhs_add_correction", d_skipRHSsetCorrection);

      if(d_symmetricCorrection == false) {
        d_skipRHSaddCorrection = true;
      }

      AMP_INSIST( (params->d_db)->keyExists("number_of_ids"), "Key ''number_of_ids'' is missing!" );
      int numIds = (params->d_db)->getInteger("number_of_ids");

      d_boundaryIds.resize(numIds);
      d_dofIds.resize(numIds);

      char key[100];
      for(int j = 0; j < numIds; j++) {
        sprintf(key, "id_%d", j);
        AMP_INSIST( (params->d_db)->keyExists(key), "Key is missing!" );
        d_boundaryIds[j] = (params->d_db)->getInteger(key);

        sprintf(key, "number_of_dofs_%d", j);
        AMP_INSIST( (params->d_db)->keyExists(key), "Key is missing!" );
        int numDofIds = (params->d_db)->getInteger(key);

        d_dofIds[j].resize(numDofIds);
        for(int i = 0; i < numDofIds; i++) {
          sprintf(key, "dof_%d_%d", j, i);
          AMP_INSIST( (params->d_db)->keyExists(key), "Key is missing!" );
          d_dofIds[j][i] = (params->d_db)->getInteger(key);
        }//end for i
      }//end for j

      if(!d_skipRHSsetCorrection) {
        d_dirichletValues.resize(numIds);
        for(int j = 0; j < numIds; j++) {
          int numDofIds = d_dofIds[j].size();
          d_dirichletValues[j].resize(numDofIds);

          for(int i = 0; i < numDofIds; i++) {
            sprintf(key, "value_%d_%d", j, i);
            AMP_INSIST( (params->d_db)->keyExists(key), "Key is missing!" );
            d_dirichletValues[j][i] = (params->d_db)->getDouble(key);
          }//end for i
        }//end for j
      }
    }
  }

  void DirichletMatrixCorrection :: computeRHScorrection(const boost::shared_ptr<DirichletMatrixCorrectionParameters> & params)
  {
    AMP_INSIST( (((params->d_db).get()) != NULL), "NULL database" );

    if(!d_skipRHSsetCorrection) {
      int numIds = d_dofIds.size();
      char key[100];
      boost::shared_ptr<AMP::InputDatabase> tmp_db(new AMP::InputDatabase("Dummy"));
      tmp_db->putBool("skip_params", false);
      tmp_db->putBool("isAttachedToVolumeOperator", false);      
      tmp_db->putInteger("number_of_ids", numIds);
      tmp_db->putInteger("print_info_level", d_iDebugPrintInfoLevel);      
      for(int j = 0; j < numIds; j++) {
        int numDofIds = d_dofIds[j].size();

        sprintf(key, "id_%d", j);
        tmp_db->putInteger(key, d_boundaryIds[j]);

        sprintf(key, "number_of_dofs_%d", j);
        tmp_db->putInteger(key, numDofIds);

        for(int i = 0; i < numDofIds; i++) {
          sprintf(key, "dof_%d_%d", j, i);
          tmp_db->putInteger(key, d_dofIds[j][i]);

          sprintf(key, "value_%d_%d", j, i);
          tmp_db->putDouble(key, d_dirichletValues[j][i]);
        }//end for i
      }//end for j

      boost::shared_ptr<DirichletVectorCorrectionParameters> setDispOpParams(new DirichletVectorCorrectionParameters(tmp_db));
      setDispOpParams->d_variable = d_variable;
      setDispOpParams->d_Mesh = d_Mesh;

      if(d_rhsCorrectionSet.get() == NULL) {
        d_rhsCorrectionSet.reset(new DirichletVectorCorrection(setDispOpParams));
      } else {
        d_rhsCorrectionSet->reset(setDispOpParams);
      }

      if(!d_skipRHSaddCorrection) {
AMP_ERROR("DirichletMatrixCorrection.cc Not fixed yet");
/*
        if(d_dispVals.get() == NULL) {
          d_dispVals = d_MeshAdapter->createVector(d_variable);
        }

        d_dispVals->zero();

        AMP::LinearAlgebra::Vector::shared_ptr emptyVec;
        d_rhsCorrectionSet->apply(emptyVec, emptyVec, d_dispVals, 1.0, 0.0);

        if(d_rhsCorrectionAdd.get() == NULL) {
          d_rhsCorrectionAdd = d_dispVals->cloneVector();
        }

        AMP::LinearAlgebra::Matrix::shared_ptr inputMatrix = params->d_inputMatrix;
        AMP_INSIST( ((inputMatrix.get()) != NULL), "NULL matrix" );

        inputMatrix->mult(d_dispVals, d_rhsCorrectionAdd);

        d_rhsCorrectionAdd->scale(-1.0);
*/
      }
    }
  }

}
}

