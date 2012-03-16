
#include <iostream>
#include <string>
#include <cassert>
#include <fstream>

#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"

#include "discretization/simpleDOF_Manager.h"
#include "vectors/VectorBuilder.h"
#include "matrices/MatrixBuilder.h"

/* AMP files */
#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/PIO.h"

void myTest(AMP::UnitTest *ut)
{
  std::string exeName("testMatSetValues");
  std::string input_file = "input_" + exeName;

  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);
  int npes = globalComm.getSize();
  int rank = globalComm.getRank();

  char outFile[256];
  sprintf(outFile, "outMatSetValues_%d", npes);

  AMP::PIO::logAllNodes(outFile);

  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  AMP_INSIST(input_db->keyExists("Mesh"), "Key ''Mesh'' is missing!");
  boost::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase("Mesh");
  boost::shared_ptr<AMP::Mesh::MeshParameters> meshParams(new AMP::Mesh::MeshParameters(mesh_db));
  meshParams->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));
  AMP::Mesh::Mesh::shared_ptr meshAdapter = AMP::Mesh::Mesh::buildMesh(meshParams);

  AMP::LinearAlgebra::Variable::shared_ptr var(new AMP::LinearAlgebra::Variable("dummy"));

  AMP::Discretization::DOFManager::shared_ptr dofMap = AMP::Discretization::simpleDOFManager::create(
      meshAdapter, AMP::Mesh::Vertex, 1, 3, true); 

  AMP::LinearAlgebra::Vector::shared_ptr vec = AMP::LinearAlgebra::createVector(dofMap, var, true);
  boost::shared_ptr<AMP::LinearAlgebra::Matrix> mat = AMP::LinearAlgebra::createMatrix(vec, vec);

  size_t locSize = vec->getLocalSize();
  size_t globSize = vec->getGlobalSize();
  size_t locStartId = vec->getLocalStartID();

  //  AMP::LinearAlgebra::Vector::shared_ptr inVec = vec->cloneVector();
  //  AMP::LinearAlgebra::Vector::shared_ptr outVec = vec->cloneVector();

  const double MatVal = 123.0;

  /*
     for(int idx = locStartId; idx < (locStartId + locSize); idx++) {
     inVec->setLocalValueByGlobalID(idx, static_cast<double>(idx + 1));
     }//end idx
     */

  AMP::plog<<"Rank = "<<rank<<": locSize = "<<locSize<<" globSize = "<<globSize<<std::endl;
  AMP::plog<<"Rank = "<<rank<<": locStartID = "<<locStartId<<std::endl;

  AMP::Mesh::MeshIterator nd = meshAdapter->getIterator(AMP::Mesh::Vertex, 0);
  //AMP::Mesh::MeshIterator nd = meshAdapter->getIDsetIterator(AMP::Mesh::Vertex, 2, 0 );
  AMP::Mesh::MeshIterator end_nd = nd.end();

  int locNonZeroCnt = 0;
  for( ; nd != end_nd; ++nd) {
    std::vector<size_t> ndDofIds;
    dofMap->getDOFs(nd->globalID(), ndDofIds);
    locNonZeroCnt += ((ndDofIds.size())*(ndDofIds.size()));

    std::vector< AMP::Mesh::MeshElement::shared_ptr > neighbors = nd->getNeighbors();
    for(size_t n = 0; n < neighbors.size(); ++n) {
      std::vector<size_t> nhDofIds;
      dofMap->getDOFs(neighbors[n]->globalID(), nhDofIds);
      locNonZeroCnt += ((ndDofIds.size())*(nhDofIds.size()));
    }//end for n
  }//end nd
  locNonZeroCnt *= 2;

  AMP::plog<<"Rank = "<<rank<<": locNonZeroCnt = "<<locNonZeroCnt<<" (NOT TRUE NZ CNT) "<<std::endl;

  std::vector<int> allNonZeroCnts(npes);
  globalComm.allGather<int>(locNonZeroCnt, &(allNonZeroCnts[0]));

  for(int proc = 0; proc < npes; proc++) {
    for(int locTestCnt = 0; locTestCnt < allNonZeroCnts[proc]; locTestCnt++) {
      mat->zero();

      int rowIdx = -1;
      int colIdx = -1;
      if(rank == proc) {
        nd = meshAdapter->getIterator(AMP::Mesh::Vertex, 0);
        //nd = meshAdapter->getIDsetIterator(AMP::Mesh::Vertex, 2, 0 );
        end_nd = nd.end();
        int cnt = 0;
        for( ; nd != end_nd; ++nd) {
          std::vector<size_t> ndDofIds;
          dofMap->getDOFs(nd->globalID(), ndDofIds);

          for(int r = 0; r < ndDofIds.size(); r++) {
            for(int c = 0; c < ndDofIds.size(); c++) {
              if(cnt == locTestCnt) {
                mat->setValueByGlobalID(ndDofIds[r], ndDofIds[c], MatVal);
                //mat->addValueByGlobalID(ndDofIds[r], ndDofIds[c], MatVal);
                rowIdx = ndDofIds[r];
                colIdx = ndDofIds[c];
              }
              cnt++;
              if(cnt == locTestCnt) {
                mat->setValueByGlobalID(ndDofIds[c], ndDofIds[r], MatVal);
                //mat->addValueByGlobalID(ndDofIds[c], ndDofIds[r], MatVal);
                rowIdx = ndDofIds[c];
                colIdx = ndDofIds[r];
              }
              cnt++;
            }//end c
          }//end r

          std::vector< AMP::Mesh::MeshElement::shared_ptr > neighbors = nd->getNeighbors();

          for(size_t n = 0; n < neighbors.size(); ++n) {
            std::vector<size_t> nhDofIds;
            dofMap->getDOFs(neighbors[n]->globalID(), nhDofIds);

            for(int r = 0; r < ndDofIds.size(); r++) {
              for(int c = 0; c < nhDofIds.size(); c++) {
                if(cnt == locTestCnt) {
                  mat->setValueByGlobalID(ndDofIds[r], nhDofIds[c], MatVal);
                  //mat->addValueByGlobalID(ndDofIds[r], nhDofIds[c], MatVal);
                  rowIdx = ndDofIds[r];
                  colIdx = nhDofIds[c];
                }
                cnt++;
                if(cnt == locTestCnt) {
                  mat->setValueByGlobalID(nhDofIds[c], ndDofIds[r], MatVal);
                  //mat->addValueByGlobalID(nhDofIds[c], ndDofIds[r], MatVal);
                  rowIdx = nhDofIds[c];
                  colIdx = ndDofIds[r];
                }
                cnt++;
              }//end c
            }//end r
          }//end for n
        }//end nd
      }

      mat->makeConsistent();

      //     outVec->setRandomValues();
      //     mat->mult(inVec, outVec);

      globalComm.bcast<int>(&rowIdx, 1, proc);
      globalComm.bcast<int>(&colIdx, 1, proc);

      std::cout<<"Testing: "<<rowIdx<<" : "<<colIdx<<std::endl;

      for(int idx = locStartId; idx < (locStartId + locSize); idx++) {
        std::vector<unsigned int> cols;
        std::vector<double> vals;
        mat->getRowByGlobalID(idx, cols, vals);
        if(idx == rowIdx) {
          bool found = false;
          for(int i = 0; i < cols.size(); i++) {
            if(cols[i] == colIdx) {
              AMP_ASSERT(vals[i] == MatVal);
              found = true;
            } else {
              AMP_ASSERT(vals[i] == 0.0);
            }
          }//end i
          AMP_ASSERT(found);
        } else {
          for(int i = 0; i < cols.size(); i++) {
            AMP_ASSERT(vals[i] == 0.0);
          }//end i
        }
      }//end idx

      /*
         for(int idx = locStartId; idx < (locStartId + locSize); idx++) {
         if(idx == rowIdx) {
         AMP_ASSERT((outVec->getLocalValueByGlobalID(idx)) == (MatVal*(static_cast<double>(colIdx + 1))));
         } else {
         AMP_ASSERT((outVec->getLocalValueByGlobalID(idx)) == 0.0);
         }
         }//end idx
         */

    }//end locTestCnt
  }//end proc

  ut->passes(exeName);
}

int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  AMP::UnitTest ut;

  try {
    myTest(&ut);
  } catch (std::exception &err) {
    std::cout << "ERROR: While testing "<<argv[0] << err.what() << std::endl;
    ut.failure("ERROR: While testing");
  } catch( ... ) {
    std::cout << "ERROR: While testing "<<argv[0] << "An unknown exception was thrown." << std::endl;
    ut.failure("ERROR: While testing");
  }

  ut.report();
  int num_failed = ut.NumFailGlobal();

  AMP::AMPManager::shutdown();
  return num_failed;
}   



