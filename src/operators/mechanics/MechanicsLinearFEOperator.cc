
#include "MechanicsLinearFEOperator.h"
#include "utils/Utilities.h"

namespace AMP {
  namespace Operator {

    MechanicsLinearFEOperator :: MechanicsLinearFEOperator (
        const boost::shared_ptr<MechanicsLinearFEOperatorParameters> & params)
      : LinearFEOperator (params) {
        AMP_INSIST( ((params.get()) != NULL), "NULL parameter" );

        d_useUpdatedLagrangian = (params->d_db)->getBoolWithDefault("USE_UPDATED_LAGRANGIAN", false);

        if(d_useUpdatedLagrangian) {
          d_mechLinULElem = boost::dynamic_pointer_cast<MechanicsLinearUpdatedLagrangianElement>(d_elemOp);
        } else {
          d_mechLinElem = boost::dynamic_pointer_cast<MechanicsLinearElement>(d_elemOp);
        }

        if(d_useUpdatedLagrangian) {
          AMP_INSIST( ((d_mechLinULElem.get()) != NULL), "d_elemOp is not of type MechanicsLinearUpdatedLagrangianElement" );
        } else {
          AMP_INSIST( ((d_mechLinElem.get()) != NULL), "d_elemOp is not of type MechanicsLinearElement" );
        }

        d_materialModel = params->d_materialModel;

        AMP_INSIST( params->d_db->keyExists("InputVariable"), "key not found" );
        std::string inpVarName = params->d_db->getString("InputVariable");
        d_inpVariable.reset(new AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable, 3>(inpVarName, d_MeshAdapter) );

        AMP_INSIST( params->d_db->keyExists("OutputVariable"), "key not found" );
        std::string outVarName = params->d_db->getString("OutputVariable");
        d_outVariable.reset(new AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable, 3>(outVarName, d_MeshAdapter) );

        if(d_useUpdatedLagrangian) {
          d_refXYZ = d_MeshAdapter->createVector(d_inpVariable);
          d_refXYZ->zero();
          
          AMP::Mesh::DOFMap::shared_ptr dof_map = d_MeshAdapter->getDOFMap(d_inpVariable);

          AMP::Mesh::MeshManager::Adapter::ElementIterator  el = d_MeshAdapter->beginElement();
          AMP::Mesh::MeshManager::Adapter::ElementIterator  end_el = d_MeshAdapter->endElement();

          for( ; el != end_el; ++el) {
            for(unsigned int i = 0; i < 3; i++) {
              dof_map->getDOFs (*el, d_dofIndices[i], i);
            }//end for i

            const ::Elem* elemPtr = &(el->getElem());

            unsigned int numNodesInCurrElem = el->numNodes();

            std::vector<double> elementRefXYZ;
            elementRefXYZ.resize(3 * numNodesInCurrElem);

            d_mechLinULElem->initializeForCurrentElement( elemPtr, d_materialModel );
            d_mechLinULElem->initializeReferenceXYZ(elementRefXYZ);

            for(unsigned int i = 0; i < 3; i++) {
              for(unsigned int j = 0; j < numNodesInCurrElem; j++) {
                d_refXYZ->setValueByGlobalID(d_dofIndices[i][j], elementRefXYZ[(3 * j) + i]);
              } // end of j
            } // end of i

          } // end of el

          d_refXYZ->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );

        } // end of UpdatedLagrangian condition.

        bool isAttachedToNonlinearOperator = params->d_db->getBoolWithDefault("isAttachedToNonlinearOperator", false);

        if(isAttachedToNonlinearOperator) {
          bool isNonlinearOperatorInitialized = params->d_db->getBoolWithDefault("isNonlinearOperatorInitialized", false);
          if(isNonlinearOperatorInitialized) {
            reset(params);
          } else {
            d_matrix = d_MeshAdapter->createMatrix ( d_inpVariable, d_outVariable );
          }
        } else {
          reset(params);
        }

      }

    void MechanicsLinearFEOperator :: preAssembly(const boost::shared_ptr<OperatorParameters>& oparams) 
    {
      if(d_useUpdatedLagrangian) {
        boost::shared_ptr<MechanicsLinearFEOperatorParameters> params = boost::dynamic_pointer_cast<MechanicsLinearFEOperatorParameters>(oparams);

        AMP_INSIST( (params != NULL), "NULL params" );

        if((d_dispVec == NULL) and (params->d_dispVec != NULL)) {
          d_dispVec = (params->d_dispVec)->cloneVector();
        }
        if(d_dispVec != NULL) {
          if(params->d_dispVec != NULL) {
            d_dispVec->copyVector(params->d_dispVec);
            d_dispVec->makeConsistent(AMP::LinearAlgebra::Vector::CONSISTENT_SET);
          } else {
            d_dispVec.reset();
          }
        }

      }

      d_matrix->zero();

      d_materialModel->preLinearAssembly();
    }

    void MechanicsLinearFEOperator :: postAssembly()
    {
      d_materialModel->postLinearAssembly();

      d_matrix->makeConsistent ();
    }

    void MechanicsLinearFEOperator :: preElementOperation( const AMP::Mesh::MeshManager::Adapter::Element & elem, 
        const std::vector<AMP::Mesh::DOFMap::shared_ptr> & dof_maps ) {
      unsigned int num_local_dofs = 0;
      for(unsigned int i = 0; i < 3; i++) {
        (dof_maps[0])->getDOFs (elem, d_dofIndices[i], i);
        num_local_dofs += d_dofIndices[i].size();
      }//end for i

      std::vector<std::vector<double> > elementInputVectors(Mechanics::TOTAL_NUMBER_OF_VARIABLES);
      std::vector<double> elementRefXYZ;

      if(d_useUpdatedLagrangian) {
        d_numNodesForCurrentElement = elem.numNodes();

        elementInputVectors[Mechanics::DISPLACEMENT].resize(num_local_dofs);
        elementRefXYZ.resize(num_local_dofs);

        for(unsigned int r = 0; r < d_numNodesForCurrentElement; r++) {
          for(unsigned int d = 0; d < 3; d++) {
            if(d_dispVec != NULL) {
              elementInputVectors[Mechanics::DISPLACEMENT][(3*r) + d] = d_dispVec->getValueByGlobalID( d_dofIndices[d][r] );
            } else {
              elementInputVectors[Mechanics::DISPLACEMENT][(3*r) + d] = 0.0;
            }
            elementRefXYZ[(3 * r) + d] = d_refXYZ->getValueByGlobalID(d_dofIndices[d][r]);
          }
        }
      }

      d_elementStiffnessMatrix.resize(num_local_dofs);
      for(unsigned int r = 0; r < num_local_dofs; r++) {
        d_elementStiffnessMatrix[r].resize(num_local_dofs);
        for(unsigned int c = 0; c < num_local_dofs; c++) {
          d_elementStiffnessMatrix[r][c] = 0.0;
        }
      }

      const ::Elem* elemPtr = &(elem.getElem());

      if(!d_useUpdatedLagrangian) {
        d_mechLinElem->initializeForCurrentElement( elemPtr, d_materialModel );
        d_mechLinElem->setElementStiffnessMatrix( d_elementStiffnessMatrix );
      } else {
        d_mechLinULElem->initializeForCurrentElement( elemPtr, d_materialModel );
        d_mechLinULElem->setElementVectors( elementInputVectors );
        d_mechLinULElem->setElementStiffnessMatrix( d_elementStiffnessMatrix );
        d_mechLinULElem->assignReferenceXYZ(elementRefXYZ);
      }
    }

    void MechanicsLinearFEOperator :: postElementOperation() {
      unsigned int num_local_dofs_per_component = d_dofIndices[0].size();

      for(unsigned int r = 0; r < num_local_dofs_per_component; r++) {
        for(unsigned int dr = 0; dr < 3; dr++) {
          for(unsigned int c = 0; c < num_local_dofs_per_component; c++) {
            for(unsigned int dc = 0; dc < 3; dc++) {
              d_matrix->addValueByGlobalID( d_dofIndices[dr][r], d_dofIndices[dc][c], 
                  d_elementStiffnessMatrix[(3*r) + dr][(3*c) + dc] );
            }
          }
        }
      }

    }

    void MechanicsLinearFEOperator :: printStressAndStrain(AMP::LinearAlgebra::Vector::shared_ptr disp, 
        const std::string & fname) {

      AMP::Mesh::DOFMap::shared_ptr dof_map = d_MeshAdapter->getDOFMap(d_inpVariable);

      AMP::Mesh::MeshManager::Adapter::ElementIterator  el = d_MeshAdapter->beginElement();
      AMP::Mesh::MeshManager::Adapter::ElementIterator  end_el = d_MeshAdapter->endElement();

      FILE* fp = fopen(fname.c_str(), "w");

      fprintf(fp, "x, y, z, Stresses(11, 22, 33, 23, 13, 12), Strains(11, 22, 33, 23, 13, 12) \n\n");

      disp->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );

      d_materialModel->preLinearAssembly();

      for( ; el != end_el; ++el) {
        for(unsigned int i = 0; i < 3; i++) {
          dof_map->getDOFs (*el, d_dofIndices[i], i);
        }//end for i

        const ::Elem* elemPtr = &(el->getElem());

        unsigned int numNodesInCurrElem = el->numNodes(); 

        std::vector<double> elementInputVector(3*numNodesInCurrElem);

        for(unsigned int r = 0; r < numNodesInCurrElem; r++) {
          for(unsigned int d = 0; d < 3; d++) {
            elementInputVector[(3*r) + d] = disp->getValueByGlobalID( d_dofIndices[d][r] );
          }
        }

        d_mechLinElem->initializeForCurrentElement( elemPtr, d_materialModel );

        d_mechLinElem->printStressAndStrain(fp, elementInputVector);

      }//end for el

      d_materialModel->postLinearAssembly();

      fprintf(fp, "\n\n");

      fclose(fp);
    }

    void MechanicsLinearFEOperator :: computeStressesAndStrains(AMP::LinearAlgebra::Vector::shared_ptr disp,
        AMP::LinearAlgebra::Vector::shared_ptr & stress, AMP::LinearAlgebra::Vector::shared_ptr & strain) {

      unsigned int numGaussPts = d_mechLinElem->getNumberOfGaussPoints();
      AMP::LinearAlgebra::Variable::shared_ptr stressVar ( new AMP::Mesh::RunTimeIntegrationPointVariable ( "stress" , (6*numGaussPts) ) );
      AMP::LinearAlgebra::Variable::shared_ptr strainVar ( new AMP::Mesh::RunTimeIntegrationPointVariable ( "strain" , (6*numGaussPts) ) );

      stress =  d_MeshAdapter->createVector ( stressVar );
      strain =  d_MeshAdapter->createVector ( strainVar );

      AMP::Mesh::DOFMap::shared_ptr gaussPt_dof_map = d_MeshAdapter->getDOFMap(stressVar);

      AMP::Mesh::DOFMap::shared_ptr dof_map = d_MeshAdapter->getDOFMap(d_inpVariable);

      AMP::Mesh::MeshManager::Adapter::ElementIterator  el = d_MeshAdapter->beginElement();
      AMP::Mesh::MeshManager::Adapter::ElementIterator  end_el = d_MeshAdapter->endElement();

      disp->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );

      d_materialModel->preLinearAssembly();

      for( ; el != end_el; ++el) {
        for(unsigned int i = 0; i < 3; i++) {
          dof_map->getDOFs (*el, d_dofIndices[i], i);
        }//end for i

        std::vector<unsigned int> gaussPtIndices;
        gaussPt_dof_map->getDOFs (*el, gaussPtIndices);

        const ::Elem* elemPtr = &(el->getElem());

        unsigned int numNodesInCurrElem = el->numNodes(); 

        std::vector<double> elementInputVector(3*numNodesInCurrElem);

        std::vector<double> elementStressVector(6*numGaussPts);
        std::vector<double> elementStrainVector(6*numGaussPts);

        for(unsigned int r = 0; r < numNodesInCurrElem; r++) {
          for(unsigned int d = 0; d < 3; d++) {
            elementInputVector[(3*r) + d] = disp->getValueByGlobalID( d_dofIndices[d][r] );
          }
        }

        d_mechLinElem->initializeForCurrentElement( elemPtr, d_materialModel );

        d_mechLinElem->computeStressAndStrain(elementInputVector, elementStressVector, elementStrainVector);

        for(unsigned int i = 0; i < (6*numGaussPts); i++) {
          stress->setValueByGlobalID(gaussPtIndices[i], elementStressVector[i]);
          strain->setValueByGlobalID(gaussPtIndices[i], elementStrainVector[i]);
        }//end for i
      }//end for el

      d_materialModel->postLinearAssembly();
    }

  }
}//end namespace


