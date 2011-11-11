/*
#include "MultiPelletMechanicsSolver.h"
#include "utils/AMP_MPI.h"
#include "operators/map/NodeToNodeMap.h"
#include "LinearOperator.h"
#include <cassert>

namespace AMP {
namespace Solver {

  void MultiPelletMechanicsSolver :: solveSerial(boost::shared_ptr<AMP::LinearAlgebra::Vector> f, 
      boost::shared_ptr<AMP::LinearAlgebra::Vector>  u) {
    int numSolvers = d_columnSolver->getNumberOfSolvers();

    //First pellet
    for(int s = 0; s < numSolvers; s++) {
      if(d_meshIds[s] == 1) {
        boost::shared_ptr<AMP::Solver::SolverStrategy> currSolver = d_columnSolver->getSolver(s);
        boost::shared_ptr<AMP::Operator::Operator> currOp = d_columnOperator->getOperator(s);
        AMP::LinearAlgebra::Variable::shared_ptr inputVar = currOp->getInputVariable();
        AMP::LinearAlgebra::Variable::shared_ptr outputVar = currOp->getOutputVariable();
        AMP::LinearAlgebra::Vector::shared_ptr subUvec = u->subsetVectorForVariable(inputVar);
        AMP::LinearAlgebra::Vector::shared_ptr subFvec = f->subsetVectorForVariable(outputVar);
        currSolver->solve(subFvec, subUvec);

        if(d_iDebugPrintInfoLevel > 1) {
          AMP::pout<<"Pellet #"<<1<<" After solve U L2-norm = "<<std::setprecision(15)<<subUvec->L2Norm()<<std::endl;
          AMP::pout<<"Pellet #"<<1<<" After solve U Max-norm = "<<std::setprecision(15)<<subUvec->maxNorm()<<std::endl;
        }

        break;
      }//end if
    }//end for s

    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    int numMaps = d_n2nmaps->getNumberOfOperators();
    for(size_t id = 2; id <= d_totalNumberOfPellets; id++) {
      d_frozenVector->zero();

      for(int s = 0; s < numSolvers; s++) {
        if(d_meshIds[s] == (id - 1)) {
          boost::shared_ptr<AMP::Operator::Operator> currOp = d_columnOperator->getOperator(s);
          AMP::LinearAlgebra::Variable::shared_ptr currOpVar = currOp->getInputVariable();
          for(int m = 0; m < numMaps; m++) {
            boost::shared_ptr<AMP::Operator::NodeToNodeMap> currMap = boost::dynamic_pointer_cast<
              AMP::Operator::NodeToNodeMap>(d_n2nmaps->getOperator(m));
            AMP::LinearAlgebra::Variable::shared_ptr currMapVar = currMap->getInputVariable();
            if((*currMapVar) == (*currOpVar)) {
              AMP_ASSERT(currMap->getMeshAdapter() == d_meshAdapters[s]);
              // if(currMap->getMeshAdapter() == d_meshAdapters[s])
              {
                if(currMap->isMaster() == true) {
                  currMap->applyStart(nullVec, u, nullVec, 1.0, 0.0);
                  break;
                }
              }
            }
          }//end for m
          break;
        }//end if
      }//end for s

      for(int s = 0; s < numSolvers; s++) {
        if(d_meshIds[s] == id) {
          boost::shared_ptr<AMP::Operator::Operator> currOp = d_columnOperator->getOperator(s);
          AMP::LinearAlgebra::Variable::shared_ptr currOpVar = currOp->getInputVariable();
          for(int m = 0; m < numMaps; m++) {
            boost::shared_ptr<AMP::Operator::NodeToNodeMap> currMap = boost::dynamic_pointer_cast<
              AMP::Operator::NodeToNodeMap>(d_n2nmaps->getOperator(m));
            AMP::LinearAlgebra::Variable::shared_ptr currMapVar = currMap->getInputVariable();
            if((*currMapVar) == (*currOpVar)) {
              AMP_ASSERT(currMap->getMeshAdapter() == d_meshAdapters[s]);
              // if(currMap->getMeshAdapter() == d_meshAdapters[s])
              {
                if(currMap->isMaster() == false) {
                  currMap->applyStart(nullVec, u, nullVec, 1.0, 0.0);
                  break;
                }
              }
            }
          }//end for m
          break;
        }
      }//end for s

      for(int s = 0; s < numSolvers; s++) {
        if(d_meshIds[s] == (id - 1)) {
          boost::shared_ptr<AMP::Operator::Operator> currOp = d_columnOperator->getOperator(s);
          AMP::LinearAlgebra::Variable::shared_ptr currOpVar = currOp->getInputVariable();
          for(int m = 0; m < numMaps; m++) {
            boost::shared_ptr<AMP::Operator::NodeToNodeMap> currMap = boost::dynamic_pointer_cast<AMP::Operator::NodeToNodeMap>(d_n2nmaps->getOperator(m));
            AMP::LinearAlgebra::Variable::shared_ptr currMapVar = currMap->getInputVariable();
            if((*currMapVar) == (*currOpVar)) {
              AMP_ASSERT(currMap->getMeshAdapter() == d_meshAdapters[s]);
              // if(currMap->getMeshAdapter() == d_meshAdapters[s])
              {
                if(currMap->isMaster() == true) {
                  currMap->applyFinish(nullVec, u, nullVec, 1.0, 0.0);
                  break;
                }
              }
            }
          }//end for m
          break;
        }
      }//end for s

      for(int s = 0; s < numSolvers; s++) {
        if(d_meshIds[s] == id) {
          boost::shared_ptr<AMP::Operator::Operator> currOp = d_columnOperator->getOperator(s);
          AMP::LinearAlgebra::Variable::shared_ptr currOpVar = currOp->getInputVariable();
          for(int m = 0; m < numMaps; m++) {
            boost::shared_ptr<AMP::Operator::NodeToNodeMap> currMap = boost::dynamic_pointer_cast<AMP::Operator::NodeToNodeMap>(d_n2nmaps->getOperator(m));
            AMP::LinearAlgebra::Variable::shared_ptr currMapVar = currMap->getInputVariable();
            if((*currMapVar) == (*currOpVar)) {
              AMP_ASSERT(currMap->getMeshAdapter() == d_meshAdapters[s]);
              // if(currMap->getMeshAdapter() == d_meshAdapters[s]) 
              {
                if(currMap->isMaster() == false) {
                  currMap->applyFinish(nullVec, u, nullVec, 1.0, 0.0);
                  break;
                }
              }
            }
          }//end for m
          break;
        }
      }//end for s

      for(int s = 0; s < numSolvers; s++) {
        if(d_meshIds[s] == id) {
          boost::shared_ptr<AMP::Solver::SolverStrategy> currSolver = d_columnSolver->getSolver(s);
          boost::shared_ptr<AMP::Operator::Operator> currOp = d_columnOperator->getOperator(s);
          if(!d_symmetricCorrection) {
            AMP_ASSERT(currOp == currSolver->getOperator());
          }
          AMP::LinearAlgebra::Variable::shared_ptr inputVar = currOp->getInputVariable();
          AMP::LinearAlgebra::Variable::shared_ptr outputVar = currOp->getOutputVariable();
          AMP::LinearAlgebra::Vector::shared_ptr subFrozenVec = d_frozenVector->subsetVectorForVariable(inputVar);
          AMP::LinearAlgebra::Vector::shared_ptr subUvec = u->subsetVectorForVariable(inputVar);
          AMP::LinearAlgebra::Vector::shared_ptr subFvec = f->subsetVectorForVariable(outputVar);
          AMP::LinearAlgebra::Vector::shared_ptr subFcopyVec = subFvec->cloneVector();
          AMP::LinearAlgebra::Vector::shared_ptr subUcopyVec = subUvec->cloneVector();
          AMP::Mesh::DOFMap::shared_ptr dof_map = d_meshAdapters[s]->getDOFMap(inputVar);
          unsigned int numIds = d_boundaryIds[s].size();
          subUvec->zero();
          if(d_symmetricCorrection) {
            for(unsigned int j = 0; j < numIds; j++) {
              AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator bnd = 
                d_meshAdapters[s]->beginOwnedBoundary( d_boundaryIds[s][j] );
              AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator end_bnd = 
                d_meshAdapters[s]->endOwnedBoundary( d_boundaryIds[s][j] );

              for( ; bnd != end_bnd; ++bnd) {
                std::vector<unsigned int> bndGlobalIds;
                dof_map->getDOFs(*bnd, bndGlobalIds, d_dofIds[s][j]);

                for(unsigned int i = 0; i < bndGlobalIds.size(); i++) {
                  double dVal = subFrozenVec->getLocalValueByGlobalID( bndGlobalIds[i] );
                  double fVal = subFvec->getLocalValueByGlobalID( bndGlobalIds[i] );
                  subUvec->setLocalValueByGlobalID(bndGlobalIds[i], (dVal + fVal));
                }//end i
              }//end bnd
            }//end j
            currOp->apply(nullVec, subUvec, subUcopyVec, 1.0, 0.0);
            subFcopyVec->subtract ( subFvec , subUcopyVec );
          } else {
            subFcopyVec->copyVector ( subFvec );
          }
          for(unsigned int j = 0; j < numIds; j++) {
            AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator bnd = 
              d_meshAdapters[s]->beginOwnedBoundary( d_boundaryIds[s][j] );
            AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator end_bnd = 
              d_meshAdapters[s]->endOwnedBoundary( d_boundaryIds[s][j] );

            for( ; bnd != end_bnd; ++bnd) {
              std::vector<unsigned int> bndGlobalIds;
              dof_map->getDOFs(*bnd, bndGlobalIds, d_dofIds[s][j]);

              for(unsigned int i = 0; i < bndGlobalIds.size(); i++) {
                double dVal = subFrozenVec->getLocalValueByGlobalID( bndGlobalIds[i] );
                double fVal = subFvec->getLocalValueByGlobalID( bndGlobalIds[i] );
                subFcopyVec->setLocalValueByGlobalID(bndGlobalIds[i], (dVal + fVal));
              }//end i
            }//end bnd
          }//end j

          if(d_iDebugPrintInfoLevel > 1) {
            AMP::LinearAlgebra::Vector::shared_ptr opDiag = ((boost::dynamic_pointer_cast<AMP::Operator::LinearOperator>(currOp))->
                getMatrix())->extractDiagonal();

            double opDiagL1 = opDiag->L1Norm();
            double opDiagMax = opDiag->maxNorm();
            double opDiagL2 = opDiag->L2Norm();
            AMP::pout<<"Pellet #"<<id<<" Diag L1 = "<<std::setprecision(15)<<opDiagL1<<std::endl;
            AMP::pout<<"Pellet #"<<id<<" Diag Max = "<<std::setprecision(15)<<opDiagMax<<std::endl;
            AMP::pout<<"Pellet #"<<id<<" Diag L2 = "<<std::setprecision(15)<<opDiagL2<<std::endl;

            AMP::pout<<"Pellet #"<<id<<" frozen L2-norm = "<<std::setprecision(15)<<subFrozenVec->L2Norm()<<std::endl;
            AMP::pout<<"Pellet #"<<id<<" frozen max-norm = "<<std::setprecision(15)<<subFrozenVec->maxNorm()<<std::endl;

            AMP::pout<<"Pellet #"<<id<<" f L2-norm = "<<std::setprecision(15)<<subFvec->L2Norm()<<std::endl;

            AMP::pout<<"Pellet #"<<id<<" fcopy L2-norm = "<<std::setprecision(15)<<subFcopyVec->L2Norm()<<std::endl;
          }

          currSolver->solve(subFcopyVec, subUvec);

          if(d_iDebugPrintInfoLevel > 1) {
            AMP::pout<<"Pellet #"<<id<<" After solve U L2-norm = "<<std::setprecision(15)<<subUvec->L2Norm()<<std::endl;
            AMP::pout<<"Pellet #"<<id<<" After solve U Max-norm = "<<std::setprecision(15)<<subUvec->maxNorm()<<std::endl;
          }

          if(d_iDebugPrintInfoLevel > 1) {
            if(d_symmetricCorrection)
            {
              subFcopyVec->copyVector(subFvec);
              for(unsigned int j = 0; j < numIds; j++) {
                AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator bnd =
                  d_meshAdapters[s]->beginOwnedBoundary( d_boundaryIds[s][j] );
                AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator end_bnd = 
                  d_meshAdapters[s]->endOwnedBoundary( d_boundaryIds[s][j] );

                for( ; bnd != end_bnd; ++bnd) {
                  std::vector<unsigned int> bndGlobalIds;
                  dof_map->getDOFs(*bnd, bndGlobalIds, d_dofIds[s][j]);

                  for(unsigned int i = 0; i < bndGlobalIds.size(); i++) {
                    double dVal = subFrozenVec->getLocalValueByGlobalID( bndGlobalIds[i] );
                    double fVal = subFvec->getLocalValueByGlobalID( bndGlobalIds[i] );
                    subFcopyVec->setLocalValueByGlobalID(bndGlobalIds[i], (dVal + fVal));
                  }//end i
                }//end bnd
              }//end j
            }
            currOp->apply(nullVec, subUvec, subUcopyVec, 1.0, 0.0);
            AMP::LinearAlgebra::Vector::shared_ptr errVec = subFcopyVec->cloneVector();
            errVec->subtract ( subFcopyVec , subUcopyVec );
            double maxErrNorm = errVec->maxNorm();
            double l2ErrNorm = errVec->L2Norm();
            AMP::pout<<"Pellet #"<<id<<" L2Err = "<<std::setprecision(15)<<
              l2ErrNorm<<" maxErr = "<<std::setprecision(15)<<maxErrNorm<<std::endl;
          }

          break;
        } 
      }//end for s
    }//end for id

    if(d_iDebugPrintInfoLevel > 2) {
      AMP::pout << "L2Norm of solution after serial multi-pellet preconditioning solve"<<std::setprecision(15) << u->L2Norm() << std::endl;
    }
  }

  void MultiPelletMechanicsSolver :: solveScan3(boost::shared_ptr<AMP::LinearAlgebra::Vector> f, 
      boost::shared_ptr<AMP::LinearAlgebra::Vector>  u) {
    int numSolvers = d_columnSolver->getNumberOfSolvers();
    AMP::LinearAlgebra::Vector::shared_ptr nullVec;

    for(int s = 0; s < numSolvers; s++) {
      boost::shared_ptr<AMP::Solver::SolverStrategy> currSolver = d_columnSolver->getSolver(s);
      boost::shared_ptr<AMP::Operator::Operator> currOp = d_columnOperator->getOperator(s);
      AMP::LinearAlgebra::Variable::shared_ptr inputVar = currOp->getInputVariable();
      AMP::LinearAlgebra::Variable::shared_ptr outputVar = currOp->getOutputVariable();
      AMP::LinearAlgebra::Vector::shared_ptr subUvec = u->subsetVectorForVariable(inputVar);
      AMP::LinearAlgebra::Vector::shared_ptr subFvec = f->subsetVectorForVariable(outputVar);
      if(d_meshIds[s] == 1) { 
        currSolver->solve(subFvec, subUvec);
      } else {
        AMP::LinearAlgebra::Vector::shared_ptr subFcopyVec = subFvec->cloneVector();
        AMP::LinearAlgebra::Vector::shared_ptr subUcopyVec = subUvec->cloneVector();
        subUvec->zero();
        AMP::Mesh::DOFMap::shared_ptr dof_map = d_meshAdapters[s]->getDOFMap(inputVar);
        unsigned int numIds = d_boundaryIds[s].size();
        for(unsigned int j = 0; j < numIds; j++) {
          AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator bnd = d_meshAdapters[s]->beginOwnedBoundary( d_boundaryIds[s][j] );
          AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator end_bnd = d_meshAdapters[s]->endOwnedBoundary( d_boundaryIds[s][j] );

          for( ; bnd != end_bnd; ++bnd) {
            std::vector<unsigned int> bndGlobalIds;
            dof_map->getDOFs(*bnd, bndGlobalIds, d_dofIds[s][j]);

            for(unsigned int i = 0; i < bndGlobalIds.size(); i++) {
              double fVal = subFvec->getLocalValueByGlobalID( bndGlobalIds[i] );
              subUvec->setLocalValueByGlobalID(bndGlobalIds[i], fVal);
            }//end i
          }//end bnd
        }//end j
        currOp->apply(nullVec, subUvec, subUcopyVec, 1.0, 0.0);
        subFcopyVec->subtract ( subFvec , subUcopyVec );
        for(unsigned int j = 0; j < numIds; j++) {
          AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator bnd = d_meshAdapters[s]->beginOwnedBoundary( d_boundaryIds[s][j] );
          AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator end_bnd = d_meshAdapters[s]->endOwnedBoundary( d_boundaryIds[s][j] );

          for( ; bnd != end_bnd; ++bnd) {
            std::vector<unsigned int> bndGlobalIds;
            dof_map->getDOFs(*bnd, bndGlobalIds, d_dofIds[s][j]);

            for(unsigned int i = 0; i < bndGlobalIds.size(); i++) {
              double fVal = subFvec->getLocalValueByGlobalID( bndGlobalIds[i] );
              subFcopyVec->setLocalValueByGlobalID(bndGlobalIds[i], fVal);
            }//end i
          }//end bnd
        }//end j
        currSolver->solve(subFcopyVec, subUvec);
      }
    }//end for s

    boost::shared_ptr<AMP::LinearAlgebra::Vector> fCopy = f->cloneVector();

    for(int iter = 0; iter < d_iMaxIterations; iter++) {
      std::vector<double> corrections;
      computeScanCorrections(u, corrections);

      fCopy->copyVector(f);

      for(int s = 0; s < numSolvers; s++) {
        if(d_meshIds[s] != 1) {
          boost::shared_ptr<AMP::Operator::Operator> currOp = d_columnOperator->getOperator(s);

          AMP::LinearAlgebra::Variable::shared_ptr inputVar = currOp->getInputVariable();
          AMP::LinearAlgebra::Variable::shared_ptr outputVar = currOp->getOutputVariable();

          AMP::LinearAlgebra::Vector::shared_ptr subFvec = f->subsetVectorForVariable(outputVar);
          AMP::LinearAlgebra::Vector::shared_ptr subFcopyVec = fCopy->subsetVectorForVariable(outputVar);

          AMP::Mesh::DOFMap::shared_ptr dof_map = d_meshAdapters[s]->getDOFMap(inputVar);

          unsigned int numIds = d_boundaryIds[s].size();

          for(unsigned int j = 0; j < numIds; j++) {
            AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator bnd = d_meshAdapters[s]->beginOwnedBoundary( d_boundaryIds[s][j] );
            AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator end_bnd = d_meshAdapters[s]->endOwnedBoundary( d_boundaryIds[s][j] );

            for( ; bnd != end_bnd; ++bnd) {
              std::vector<unsigned int> bndGlobalIds;
              dof_map->getDOFs(*bnd, bndGlobalIds, d_dofIds[s][j]);

              for(unsigned int i = 0; i < bndGlobalIds.size(); i++) {
                double dVal = corrections[s];
                double fVal = subFvec->getLocalValueByGlobalID( bndGlobalIds[i] );
                subFcopyVec->setLocalValueByGlobalID(bndGlobalIds[i], (dVal + fVal));
              }//end i
            }//end bnd
          }//end j
        }//end if
      }//end s      

      for(int s = 0; s < numSolvers; s++) {
        boost::shared_ptr<AMP::Solver::SolverStrategy> currSolver = d_columnSolver->getSolver(s);
        boost::shared_ptr<AMP::Operator::Operator> currOp = d_columnOperator->getOperator(s);
        AMP::LinearAlgebra::Variable::shared_ptr inputVar = currOp->getInputVariable();
        AMP::LinearAlgebra::Variable::shared_ptr outputVar = currOp->getOutputVariable();
        AMP::LinearAlgebra::Vector::shared_ptr subUvec = u->subsetVectorForVariable(inputVar);
        AMP::LinearAlgebra::Vector::shared_ptr subFvec = fCopy->subsetVectorForVariable(outputVar);
        if(d_meshIds[s] == 1) { 
          currSolver->solve(subFvec, subUvec);
        } else {
          AMP::LinearAlgebra::Vector::shared_ptr subFcopyVec = subFvec->cloneVector();
          AMP::LinearAlgebra::Vector::shared_ptr subUcopyVec = subUvec->cloneVector();
          subUvec->zero();
          AMP::Mesh::DOFMap::shared_ptr dof_map = d_meshAdapters[s]->getDOFMap(inputVar);
          unsigned int numIds = d_boundaryIds[s].size();
          for(unsigned int j = 0; j < numIds; j++) {
            AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator bnd = d_meshAdapters[s]->beginOwnedBoundary( d_boundaryIds[s][j] );
            AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator end_bnd = d_meshAdapters[s]->endOwnedBoundary( d_boundaryIds[s][j] );

            for( ; bnd != end_bnd; ++bnd) {
              std::vector<unsigned int> bndGlobalIds;
              dof_map->getDOFs(*bnd, bndGlobalIds, d_dofIds[s][j]);

              for(unsigned int i = 0; i < bndGlobalIds.size(); i++) {
                double fVal = subFvec->getLocalValueByGlobalID( bndGlobalIds[i] );
                subUvec->setLocalValueByGlobalID(bndGlobalIds[i], fVal);
              }//end i
            }//end bnd
          }//end j
          currOp->apply(nullVec, subUvec, subUcopyVec, 1.0, 0.0);
          subFcopyVec->subtract ( subFvec , subUcopyVec );
          for(unsigned int j = 0; j < numIds; j++) {
            AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator bnd = d_meshAdapters[s]->beginOwnedBoundary( d_boundaryIds[s][j] );
            AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator end_bnd = d_meshAdapters[s]->endOwnedBoundary( d_boundaryIds[s][j] );

            for( ; bnd != end_bnd; ++bnd) {
              std::vector<unsigned int> bndGlobalIds;
              dof_map->getDOFs(*bnd, bndGlobalIds, d_dofIds[s][j]);

              for(unsigned int i = 0; i < bndGlobalIds.size(); i++) {
                double fVal = subFvec->getLocalValueByGlobalID( bndGlobalIds[i] );
                subFcopyVec->setLocalValueByGlobalID(bndGlobalIds[i], fVal);
              }//end i
            }//end bnd
          }//end j
          currSolver->solve(subFcopyVec, subUvec);
        }
      }//end for s

    }//end for iter

  }

  void MultiPelletMechanicsSolver :: solveScan2(boost::shared_ptr<AMP::LinearAlgebra::Vector> f, 
      boost::shared_ptr<AMP::LinearAlgebra::Vector>  u) {
    int numSolvers = d_columnSolver->getNumberOfSolvers();

    d_columnSolver->solve(f, u); 

    boost::shared_ptr<AMP::LinearAlgebra::Vector> fCopy = f->cloneVector();

    for(int iter = 0; iter < d_iMaxIterations; iter++) {
      std::vector<double> corrections;
      computeScanCorrections(u, corrections);

      fCopy->copyVector(f);

      for(int s = 0; s < numSolvers; s++) {
        if(d_meshIds[s] != 1) {
          boost::shared_ptr<AMP::Solver::SolverStrategy> currSolver = d_columnSolver->getSolver(s);
          boost::shared_ptr<AMP::Operator::Operator> currOp = d_columnOperator->getOperator(s);
          AMP_ASSERT(currOp == currSolver->getOperator());

          AMP::LinearAlgebra::Variable::shared_ptr inputVar = currOp->getInputVariable();
          AMP::LinearAlgebra::Variable::shared_ptr outputVar = currOp->getOutputVariable();

          AMP::LinearAlgebra::Vector::shared_ptr subFvec = f->subsetVectorForVariable(outputVar);
          AMP::LinearAlgebra::Vector::shared_ptr subFcopyVec = fCopy->subsetVectorForVariable(outputVar);

          AMP::Mesh::DOFMap::shared_ptr dof_map = d_meshAdapters[s]->getDOFMap(inputVar);

          unsigned int numIds = d_boundaryIds[s].size();

          for(unsigned int j = 0; j < numIds; j++) {
            AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator bnd = d_meshAdapters[s]->beginOwnedBoundary( d_boundaryIds[s][j] );
            AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator end_bnd = d_meshAdapters[s]->endOwnedBoundary( d_boundaryIds[s][j] );

            for( ; bnd != end_bnd; ++bnd) {
              std::vector<unsigned int> bndGlobalIds;
              dof_map->getDOFs(*bnd, bndGlobalIds, d_dofIds[s][j]);

              for(unsigned int i = 0; i < bndGlobalIds.size(); i++) {
                double dVal = corrections[s];
                double fVal = subFvec->getLocalValueByGlobalID( bndGlobalIds[i] );
                subFcopyVec->setLocalValueByGlobalID(bndGlobalIds[i], (dVal + fVal));
              }//end i
            }//end bnd
          }//end j
        }//end if
      }//end s      

      d_columnSolver->solve(fCopy, u);
    }//end for iter

  }

  void MultiPelletMechanicsSolver :: solveScan1(boost::shared_ptr<AMP::LinearAlgebra::Vector> f, 
      boost::shared_ptr<AMP::LinearAlgebra::Vector>  u) {
    d_columnSolver->solve(f, u);

    std::vector<double> corrections;
    computeScanCorrections(u, corrections);

    for(unsigned int s = 0; s < d_meshAdapters.size(); s++) {
      boost::shared_ptr<AMP::Solver::SolverStrategy> currSolver = d_columnSolver->getSolver(s);
      boost::shared_ptr<AMP::Operator::Operator> currOp = d_columnOperator->getOperator(s);
      AMP_ASSERT(currOp == currSolver->getOperator());
      AMP::LinearAlgebra::Variable::shared_ptr inputVar = currOp->getInputVariable();
      AMP::LinearAlgebra::Vector::shared_ptr subUvec = u->subsetVectorForVariable(inputVar);
      AMP::Mesh::DOFMap::shared_ptr dof_map = d_meshAdapters[s]->getDOFMap(inputVar);
      AMP::Mesh::MeshManager::Adapter::OwnedNodeIterator nd = d_meshAdapters[s]->beginOwnedNode( );
      AMP::Mesh::MeshManager::Adapter::OwnedNodeIterator end_nd = d_meshAdapters[s]->endOwnedNode( );
      for( ; nd != end_nd; ++nd) {
        std::vector<unsigned int> ndGlobalIds;
        std::vector<unsigned int> dummy(1);
        dummy[0] = 2;
        dof_map->getDOFs(*nd, ndGlobalIds, dummy);
        for(unsigned int i = 0; i < ndGlobalIds.size(); i++) {
          subUvec->addLocalValueByGlobalID( ndGlobalIds[i], corrections[s] );
        }//end i
      }//end nd
    }//end s

  }

  void MultiPelletMechanicsSolver :: computeScanCorrections(boost::shared_ptr<AMP::LinearAlgebra::Vector>  u,
      std::vector<double> & corrections) {

    //d_meshIds must be sorted on each processor. 
    //Each pellet either uses AMP_COMM_WORLD or its own communicator.
    //If a pellet uses its own communicator, then the set of processors used by
    //one pellet is disjoint from the set of processors used by any other
    //pellet. 
    assert((d_meshAdapters.size() == 1u) || (d_meshAdapters.size() == d_totalNumberOfPellets));

    std::vector<double> localNumerator(d_meshAdapters.size());
    std::vector<int> localDenominator(d_meshAdapters.size());
    for(unsigned int s = 0; s < d_meshAdapters.size(); s++) {
      localNumerator[s] = 0.0;
      localDenominator[s] = 0;
      boost::shared_ptr<AMP::Operator::Operator> currOp = d_columnOperator->getOperator(s);
      AMP::LinearAlgebra::Variable::shared_ptr inputVar = currOp->getInputVariable();
      AMP::LinearAlgebra::Vector::shared_ptr subUvec = u->subsetVectorForVariable(inputVar);
      AMP::Mesh::DOFMap::shared_ptr dof_map = d_meshAdapters[s]->getDOFMap(inputVar);
      //2 is the ID for the top surface
      AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator bnd = d_meshAdapters[s]->beginOwnedBoundary( 2 );
      AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator end_bnd = d_meshAdapters[s]->endOwnedBoundary( 2 );
      for( ; bnd != end_bnd; ++bnd) {
        std::vector<unsigned int> bndGlobalIds;
        std::vector<unsigned int> dummy(1);
        //Select Z-direction
        dummy[0] = 2;
        dof_map->getDOFs(*bnd, bndGlobalIds, dummy);
        for(unsigned int i = 0; i < bndGlobalIds.size(); i++) {
          double uVal = subUvec->getLocalValueByGlobalID( bndGlobalIds[i] );
          localNumerator[s] += uVal;
          localDenominator[s]++;
        }//end i
      }//end bnd
    }//end s

    bool isRank0;
    std::vector<double> meanZdisp(d_meshAdapters.size());
    if(d_meshAdapters.size() == 1) {
      AMP_MPI pelletComm = d_meshAdapters[0]->getComm();
      int rankInPelletComm = pelletComm.getRank();
      double meanVal = pelletComm.sumReduce(localNumerator[0]);
      int temp = pelletComm.sumReduce(localDenominator[0]);
      if(!rankInPelletComm) {
        meanZdisp[0] = meanVal/static_cast<double>(temp);
        isRank0 = true;
      } else {
        isRank0 = false;
      }
    } else {
      AMP_MPI globalComm(AMP_COMM_WORLD);
      int rank = globalComm.getRank();
      std::vector<double> meanVal(d_totalNumberOfPellets);
      std::vector<int> temp(d_totalNumberOfPellets);
      globalComm.sumReduce((double*)&(localNumerator[0]),(double*)&(meanVal[0]),d_totalNumberOfPellets);
      globalComm.sumReduce((int*)&(localDenominator[0]),(int*)&(temp[0]),d_totalNumberOfPellets);
      if(!rank) {
        for(size_t s = 0; s < d_totalNumberOfPellets; s++) {
          meanZdisp[s] = meanVal[s]/static_cast<double>(temp[s]);
        }//end s
        isRank0 = true;
      } else {
        isRank0 = false;
      }
    }

    corrections.resize(d_meshAdapters.size());
    if(d_meshAdapters.size() == 1) {
      AMP_MPI newComm;
      AMP_MPI globalComm(AMP_COMM_WORLD);
      if(isRank0) {
        newComm = globalComm.split(1,d_meshIds[0]);
        newComm.sumScan((double*)&(meanZdisp[0]),(double*)&(corrections[0]),1);
        corrections[0] = corrections[0] - meanZdisp[0];
      } else {
        newComm = globalComm.split(0,d_meshIds[0]);
      }
    } else {
      if(isRank0) {
        corrections[0] = 0;
        for(size_t i = 1; i < d_totalNumberOfPellets; i++) {
          corrections[i] = corrections[i - 1] + meanZdisp[i - 1];
        }
      }
    }

    if(d_meshAdapters.size() == 1) {
      AMP_MPI pelletComm = d_meshAdapters[0]->getComm();
      corrections[0] = pelletComm.bcast(corrections[0],0);
    } else {
      AMP_MPI globalComm(AMP_COMM_WORLD);
      globalComm.bcast((double*)&(corrections[0]),d_totalNumberOfPellets,0);
    }

  }

  void MultiPelletMechanicsSolver :: solveNoScan(boost::shared_ptr<AMP::LinearAlgebra::Vector> f, 
      boost::shared_ptr<AMP::LinearAlgebra::Vector>  u) {
    int numSolvers = d_columnSolver->getNumberOfSolvers();

    d_columnSolver->solve(f, u);

    boost::shared_ptr<AMP::LinearAlgebra::Vector> fCopy = f->cloneVector();

    for(int iter = 0; iter < d_iMaxIterations; iter++) {
      AMP::LinearAlgebra::Vector::shared_ptr nullVec;
      d_n2nmaps->apply(nullVec, u, nullVec, 1.0, 0.0);

      fCopy->copyVector(f);

      for(int s = 0; s < numSolvers; s++) {
        if(d_meshIds[s] != 1) {
          boost::shared_ptr<AMP::Solver::SolverStrategy> currSolver = d_columnSolver->getSolver(s);
          boost::shared_ptr<AMP::Operator::Operator> currOp = d_columnOperator->getOperator(s);
          AMP_ASSERT(currOp == currSolver->getOperator());

          AMP::LinearAlgebra::Variable::shared_ptr inputVar = currOp->getInputVariable();
          AMP::LinearAlgebra::Variable::shared_ptr outputVar = currOp->getOutputVariable();

          AMP::LinearAlgebra::Vector::shared_ptr subFrozenVec = d_frozenVector->subsetVectorForVariable(inputVar);
          AMP::LinearAlgebra::Vector::shared_ptr subFvec = f->subsetVectorForVariable(outputVar);
          AMP::LinearAlgebra::Vector::shared_ptr subFcopyVec = fCopy->subsetVectorForVariable(outputVar);

          AMP::Mesh::DOFMap::shared_ptr dof_map = d_meshAdapters[s]->getDOFMap(inputVar);

          unsigned int numIds = d_boundaryIds[s].size();

          for(unsigned int j = 0; j < numIds; j++) {
            AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator bnd = d_meshAdapters[s]->beginOwnedBoundary( d_boundaryIds[s][j] );
            AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator end_bnd = d_meshAdapters[s]->endOwnedBoundary( d_boundaryIds[s][j] );

            for( ; bnd != end_bnd; ++bnd) {
              std::vector<unsigned int> bndGlobalIds;
              dof_map->getDOFs(*bnd, bndGlobalIds, d_dofIds[s][j]);

              for(unsigned int i = 0; i < bndGlobalIds.size(); i++) {
                double dVal = subFrozenVec->getLocalValueByGlobalID( bndGlobalIds[i] );
                double fVal = subFvec->getLocalValueByGlobalID( bndGlobalIds[i] );
                subFcopyVec->setLocalValueByGlobalID(bndGlobalIds[i], (dVal + fVal));
              }//end i
            }//end bnd
          }//end j
        }//end if
      }//end s      

      d_columnSolver->solve(fCopy, u);
    }//end iter
  }

void MultiPelletMechanicsSolver :: reset(boost::shared_ptr<SolverStrategyParameters> params)
{
  d_columnSolver->reset(params);
}

}
}
*/
