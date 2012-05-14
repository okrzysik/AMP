
#include "operators/boundary/TractionBoundaryOperator.h"
#include "face_quad4.h"

namespace AMP {
  namespace Operator {

    TractionBoundaryOperator :: TractionBoundaryOperator(const boost::shared_ptr<TractionBoundaryOperatorParameters> & params)
      : BoundaryOperator(params) {
        AMP_INSIST( params->d_db->keyExists("InputVariable"), "key not found");
        AMP_INSIST( params->d_db->keyExists("OutputVariable"), "key not found");
        std::string inpVarName = params->d_db->getString("InputVariable");
        std::string outVarName = params->d_db->getString("OutputVariable");
        d_inputVar.reset(new AMP::LinearAlgebra::Variable(inpVarName));
        d_outputVar.reset(new AMP::LinearAlgebra::Variable(outVarName));
        d_residualMode = params->d_db->getBool("ResidualMode");
        d_boundaryId = params->d_db->getInteger("BoundaryID");
        setTraction(params->d_traction);
      }

    void TractionBoundaryOperator :: addRHScorrection(AMP::LinearAlgebra::Vector::shared_ptr rhs) {
      if(!d_residualMode) {
        AMP::LinearAlgebra::Vector::shared_ptr myRhs = mySubsetVector(rhs, d_outputVar);
        if(d_correction == NULL) {
          d_correction = myRhs->cloneVector();
        }
        computeCorrection();
        myRhs->add(myRhs, d_correction);
      }
    }

    void TractionBoundaryOperator :: apply(const AMP::LinearAlgebra::Vector::shared_ptr &, const AMP::LinearAlgebra::Vector::shared_ptr &,
        AMP::LinearAlgebra::Vector::shared_ptr &r, const double, const double) {
      if(d_residualMode) {
        AMP::LinearAlgebra::Vector::shared_ptr rInternal = mySubsetVector(r, d_outputVar);
        if(d_correction == NULL) {
          d_correction = rInternal->cloneVector();
        }
        computeCorrection();
        rInternal->subtract(rInternal, d_correction);
      }
    }


    void TractionBoundaryOperator :: computeCorrection() {
      AMP::Discretization::DOFManager::shared_ptr inpDofMap = d_traction->getDOFManager();
      AMP::Discretization::DOFManager::shared_ptr outDofMap = d_correction->getDOFManager();
      d_correction->zero();

      AMP::Mesh::MeshIterator bnd     = d_Mesh->getBoundaryIDIterator(AMP::Mesh::Face, d_boundaryId, 0);
      AMP::Mesh::MeshIterator end_bnd = bnd.end();

      for( ; bnd != end_bnd; ++bnd) {
        d_currNodes = bnd->getElements(AMP::Mesh::Vertex);
        unsigned int numNodesInCurrElem = d_currNodes.size();
        createCurrentLibMeshElement();

        std::vector<size_t> inpDofIndices;
        inpDofMap->getDOFs(bnd->globalID(), inpDofIndices);

        std::vector<std::vector<size_t> > outDofIndices(numNodesInCurrElem);
        for(unsigned int i = 0; i < numNodesInCurrElem ; ++i) {
          outDofMap->getDOFs(d_currNodes[i].globalID(), outDofIndices[i]);
        }//end i

        destroyCurrentLibMeshElement();
      }//end bnd

      d_correction->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_ADD );
    }

    AMP::LinearAlgebra::Vector::shared_ptr TractionBoundaryOperator :: mySubsetVector(AMP::LinearAlgebra::Vector::shared_ptr vec, 
        AMP::LinearAlgebra::Variable::shared_ptr var) {
      if(vec != NULL) {
        if(d_Mesh.get() != NULL) {
          AMP::LinearAlgebra::VS_Mesh meshSelector(var->getName(), d_Mesh);
          AMP::LinearAlgebra::Vector::shared_ptr meshSubsetVec = vec->select(meshSelector, var->getName());
          return meshSubsetVec->subsetVectorForVariable(var);
        } else {
          return vec->subsetVectorForVariable(var);
        }
      } else {
        return vec;
      }
    }

    void TractionBoundaryOperator :: createCurrentLibMeshElement() {
      d_currElemPtr = new ::Quad4;
      for(size_t j = 0; j < d_currNodes.size(); j++) {
        std::vector<double> pt = d_currNodes[j].coord();
        d_currElemPtr->set_node(j) = new ::Node(pt[0], pt[1], pt[2], j);
      }//end for j
    }

    void TractionBoundaryOperator :: destroyCurrentLibMeshElement() {
      for(size_t j = 0; j < d_currElemPtr->n_nodes(); j++) {
        delete (d_currElemPtr->get_node(j));
        d_currElemPtr->set_node(j) = NULL;
      }//end for j
      delete d_currElemPtr;
      d_currElemPtr = NULL;
    }

  }
}

