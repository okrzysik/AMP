
#ifndef included_AMP_NoteToSegmentConstraintsOperator
#define included_AMP_NoteToSegmentConstraintsOperator

#include "boost/shared_ptr.hpp"
#include "matrices/Matrix.h"
#include "operators/Operator.h"
#include "operators/contact/NodeToSegmentConstraintsOperatorParameters.h"
#include "vectors/Vector.h"

#include "ampmesh/dendro/DendroSearch.h"
#include <vector>

namespace AMP {
  namespace Operator {

    /**
      An abstract base class for representing a linear operator. This class 
      stores the matrix representation of the linear operator. It provides
      an implementation of the apply() function.
      @see Operator
      */
    class NodeToSegmentConstraintsOperator : public Operator 
    {

      public :

        /**
          Constructor. This resets the matrix shared pointer.
          @param [in] params 
          */
        NodeToSegmentConstraintsOperator (const boost::shared_ptr<NodeToSegmentConstraintsOperatorParameters> & params)
          : Operator(params)
        {

          d_DOFmanager = (params->d_DOFmanager);

          d_MasterMeshID = (params->d_MasterMeshID);
          d_SlaveMeshID = (params->d_SlaveMeshID);
          d_MasterBoundaryID = (params->d_MasterBoundaryID);
          d_SlaveBoundaryID = (params->d_SlaveBoundaryID);
        }

        /**
          Destructor
          */
        virtual ~NodeToSegmentConstraintsOperator() { }

        /**
         * This function is useful for re-initializing/updating an operator
         * \param params
         *        parameter object containing parameters to change
         */
        virtual void reset(const boost::shared_ptr<OperatorParameters> & params);

        /**
          The apply function for this operator, A, performs the following operation:
          r = a*A(u) + b*f, if f is not NULL and r = a*A(u), if f is NULL.
          Here, A(u) is simply a Matrix-Vector multiplication.
          @param [in] f auxillary/rhs vector. 
          @param [in] u input vector. 
          @param [out] r residual/output vector. 
          @param [in] a first constant used in the expression: r = a*A(u) + b*f. The default value is -1.
          @param [in] b second constant used in the expression: r = a*A(u) + b*f. The default value is 1.
          */
        virtual void apply(const AMP::LinearAlgebra::Vector::shared_ptr &f, const AMP::LinearAlgebra::Vector::shared_ptr &u,
            AMP::LinearAlgebra::Vector::shared_ptr  &r, const double a = -1.0, const double b = 1.0);

        virtual void applyTranspose(const AMP::LinearAlgebra::Vector::shared_ptr &f, const AMP::LinearAlgebra::Vector::shared_ptr &u,
            AMP::LinearAlgebra::Vector::shared_ptr  &r, const double a = -1.0, const double b = 1.0);

      protected :

      private :

        AMP::Discretization::DOFManager::shared_ptr d_DOFmanager;

        unsigned int d_MasterMeshID;
        unsigned int d_SlaveMeshID;

        int d_MasterBoundaryID;
        int d_SlaveBoundaryID;

        std::vector<int> d_SendCnts;
        std::vector<int> d_SendDisps;
        std::vector<int> d_RecvCnts;
        std::vector<int> d_RecvDisps;
        std::vector<int> d_TransposeSendCnts;
        std::vector<int> d_TransposeSendDisps;
        std::vector<int> d_TransposeRecvCnts;
        std::vector<int> d_TransposeRecvDisps;

        std::vector<size_t> d_SlaveGlobalIDs;
        std::vector<size_t> d_RecvMasterGlobalIDs;
        std::vector<size_t> d_MasterVerticesOwnerRanks;
        std::vector<double> d_MasterShapeFunctionsValues;
    };

  }
}

#endif


