
#ifndef included_AMP_WeldOperator
#define included_AMP_WeldOperator

#include "mpi.h"

#include "utils/InputDatabase.h"

#include "boost/shared_ptr.hpp"

#include "operators/Operator.h"
#include "vectors/Vector.h"
#include "vectors/Variable.h"
#include <string>

#ifdef DEBUG_CHECK_ASSERTIONS
#include <cassert>
#endif

namespace AMP {
  namespace Operator {

    class WeldOperator : public Operator
    {

      public :
        WeldOperator (const boost::shared_ptr<OperatorParameters> & params) : 
          Operator (params) { }

        virtual ~WeldOperator() { }

        virtual void reset(const boost::shared_ptr<OperatorParameters>& params) {
          (void) params;
        }

        virtual void apply(const AMP::LinearAlgebra::Vector::shared_ptr &f, 
            const  AMP::LinearAlgebra::Vector::shared_ptr &u, AMP::LinearAlgebra::Vector::shared_ptr  &r,
            const double a = -1.0, const double b = 1.0) {
          AMP::LinearAlgebra::Vector::shared_ptr inVec = u->subsetVectorForVariable(d_inpVar);

          AMP::Mesh::DOFMap::shared_ptr dof_map = d_MeshAdapter->getDOFMap(d_inpVar);

          AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator bnd = d_MeshAdapter->beginOwnedBoundary( d_inputBoundaryId );
          AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator end_bnd = d_MeshAdapter->endOwnedBoundary( d_inputBoundaryId );

          double val;
          int cnt = 0;
          for( ; bnd != end_bnd; ++bnd) {
            std::vector<unsigned int> bndGlobalIds;
            std::vector<unsigned int> singleton(1);
            singleton[0] = 2;
            dof_map->getDOFs(*bnd, bndGlobalIds, singleton);
            val = inVec->getLocalValueByGlobalID(bndGlobalIds[0]);
            cnt++;
          }//end for bnd
          assert((cnt == 0) or (cnt == 1));
          //d_comm is constructed so that rank 0 has the input boundary
//          MPI_Bcast(&val, 1, MPI_DOUBLE, 0, d_comm);
          val = d_comm.bcast(val, 0);
          d_outVec->setToScalar(val);
        }

        unsigned int d_inputBoundaryId;
        AMP::LinearAlgebra::Variable::shared_ptr d_inpVar;
        AMP::LinearAlgebra::Vector::shared_ptr d_outVec;
        AMP::AMP_MPI d_comm;
    };
  }
}

#endif


