#include "MassMatrixCorrection.h"
#include "utils/Utilities.h"
#include "utils/InputDatabase.h"

namespace AMP {
namespace Operator {
    
    void MassMatrixCorrection :: resetBoundaryIds(const boost::shared_ptr<MassMatrixCorrectionParameters> & params)
    {
        AMP_INSIST( (((params->d_db).get()) != NULL), "NULL database" );
        bool skipParams = (params->d_db)->getBoolWithDefault("skip_params", true);
        d_bSetIdentityOnDiagonal =  (params->d_db)->getBoolWithDefault("setIdentityOnDiagonal", false);
        
        if(!skipParams) {
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
        }
    }
    
        void MassMatrixCorrection :: reset(const boost::shared_ptr<OperatorParameters>& params)
    {


        boost::shared_ptr<MassMatrixCorrectionParameters> myParams = 
        boost::dynamic_pointer_cast<MassMatrixCorrectionParameters>(params);
        
        AMP_INSIST( ((myParams.get()) != NULL), "NULL parameters" );
        
        resetBoundaryIds(myParams);
        
        double diagVal = d_bSetIdentityOnDiagonal? 1.0:0.0;
        
        AMP::LinearAlgebra::Matrix::shared_ptr inputMatrix = myParams->d_inputMatrix;
        AMP_INSIST( ((inputMatrix.get()) != NULL), "NULL matrix" );
        
        AMP_ERROR("This needs to be fixed");
/*        AMP::Mesh::DOFMap::shared_ptr dof_map = d_MeshAdapter->getDOFMap ( d_variable );
        
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
                                inputMatrix->setValueByGlobalID ( bndGlobalIds[j], bndGlobalIds[j], diagVal );
                            } else {
                                inputMatrix->setValueByGlobalID ( bndGlobalIds[j], dofIndices[i] , 0.0 );
                                inputMatrix->setValueByGlobalID ( dofIndices[i], bndGlobalIds[j] , 0.0 );
                            }
                        }//end for i
                    }//end for j
                }//end for el
            }//end for bnd
        }//end for k
*/

        //This does consistent for both "Sum-into" and "set".
        inputMatrix->makeConsistent();
    }
    
}
}

