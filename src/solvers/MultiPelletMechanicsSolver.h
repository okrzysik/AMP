
#ifndef included_AMP_MultiPelletMechanicsSolver
#define included_AMP_MultiPelletMechanicsSolver

#include "SolverStrategy.h"
#include "ColumnSolver.h"
#include "operators/ColumnOperator.h"
//#include "operators/map/AsyncMapColumnOperator.h"

namespace AMP {
namespace Solver {

  typedef SolverStrategyParameters MultiPelletMechanicsSolverParameters;

  class MultiPelletMechanicsSolver: public SolverStrategy {
    public:
      MultiPelletMechanicsSolver() {  }

      MultiPelletMechanicsSolver(boost::shared_ptr<MultiPelletMechanicsSolverParameters> params) 
        : SolverStrategy(params) {
          d_useScan = (params->d_db)->getBoolWithDefault("USE_SCAN", false);
          d_scanFastLocalSolve = (params->d_db)->getBoolWithDefault("SCAN_WITH_FAST_LOCAL_SOLVE", false);
          d_solveAfterScan = (params->d_db)->getBoolWithDefault("SOLVE_AFTER_SCAN", false);
          d_useSerial = (params->d_db)->getBoolWithDefault("USE_SERIAL", false);
          d_resetColumnOperator = (params->d_db)->getBoolWithDefault("ResetColumnOperator", false);
          d_symmetricCorrection = (params->d_db)->getBoolWithDefault("SYMMETRIC_CORRECTION", true);
          d_columnOperator = boost::dynamic_pointer_cast<AMP::Operator::ColumnOperator>(d_pOperator);
        }

      ~MultiPelletMechanicsSolver() { }

      void solve(boost::shared_ptr<AMP::LinearAlgebra::Vector> f, boost::shared_ptr<AMP::LinearAlgebra::Vector>  u) {
        if(d_useSerial) {
          solveSerial(f, u);
        } else if(d_useScan) {
          if(d_solveAfterScan) {
            if(d_scanFastLocalSolve) {
              solveScan3(f, u);
            } else {
              solveScan2(f, u);
            }
          } else {
            solveScan1(f, u);
          }
        } else {
          solveNoScan(f, u);
        }
      }

      void solveScan1(boost::shared_ptr<AMP::LinearAlgebra::Vector> f, boost::shared_ptr<AMP::LinearAlgebra::Vector>  u);

      void solveScan2(boost::shared_ptr<AMP::LinearAlgebra::Vector> f, boost::shared_ptr<AMP::LinearAlgebra::Vector>  u);

      void solveScan3(boost::shared_ptr<AMP::LinearAlgebra::Vector> f, boost::shared_ptr<AMP::LinearAlgebra::Vector>  u);

      void computeScanCorrections(boost::shared_ptr<AMP::LinearAlgebra::Vector> u, std::vector<double> & corrections);

      void solveNoScan(boost::shared_ptr<AMP::LinearAlgebra::Vector> f, boost::shared_ptr<AMP::LinearAlgebra::Vector>  u);

      void solveSerial(boost::shared_ptr<AMP::LinearAlgebra::Vector> f, boost::shared_ptr<AMP::LinearAlgebra::Vector>  u);

      void initialize(boost::shared_ptr<SolverStrategyParameters> const parameters) { }

      void reset(boost::shared_ptr<SolverStrategyParameters> params);
      
      void resetOperator(const boost::shared_ptr<AMP::Operator::OperatorParameters> params) {
        if(d_resetColumnOperator) {
          d_columnOperator->reset(params);
        }
        d_columnSolver->resetOperator(params);
      }

      void setFrozenVector(AMP::LinearAlgebra::Vector::shared_ptr vec) { 
        d_frozenVector = vec;
      }

      //void setMaps(boost::shared_ptr<AMP::Operator::AsyncMapColumnOperator> maps) {
      //  d_n2nmaps = maps;
      //}

      void setColumnSolver(boost::shared_ptr<ColumnSolver> colSolver) {
        d_columnSolver = colSolver;
      }

      void setMeshAdapters(std::vector<AMP::Mesh::MeshManager::Adapter::shared_ptr> meshes) {
        d_meshAdapters = meshes;
      }

      void setMeshIds(std::vector<unsigned int> ids) {
        d_meshIds = ids;
      }

      void setTotalNumberOfPellets(size_t num) {
        d_totalNumberOfPellets = num;
      }

      void setBoundaryIds(std::vector<std::vector<short int> > vec) {
        d_boundaryIds = vec ;
      }

      void setDofIds(std::vector<std::vector<std::vector<unsigned int> > > vec) {
        d_dofIds = vec;
      }

    protected:

    private:
      boost::shared_ptr<ColumnSolver> d_columnSolver;
      boost::shared_ptr<AMP::Operator::ColumnOperator> d_columnOperator;
      //boost::shared_ptr<AMP::Operator::AsyncMapColumnOperator>  d_n2nmaps;
      size_t d_totalNumberOfPellets;
      AMP::LinearAlgebra::Vector::shared_ptr d_frozenVector;
      std::vector<AMP::Mesh::Mes::shared_ptr> d_mesh;
      std::vector<unsigned int> d_meshIds;
      std::vector<std::vector<short int> > d_boundaryIds;
      std::vector<std::vector<std::vector<unsigned int> > > d_dofIds;
      bool d_useScan;
      bool d_scanFastLocalSolve;
      bool d_useSerial;
      bool d_solveAfterScan;
      bool d_resetColumnOperator;
      bool d_symmetricCorrection;
  };

}
}

#endif

