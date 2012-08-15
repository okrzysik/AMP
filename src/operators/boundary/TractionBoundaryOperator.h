
#ifndef included_AMP_TractionBoundaryOperator
#define included_AMP_TractionBoundaryOperator

#include "BoundaryOperator.h"
#include "TractionBoundaryOperatorParameters.h"
#include "elem.h"

namespace AMP {
  namespace Operator {

    class TractionBoundaryOperator : public BoundaryOperator {
      public:

        TractionBoundaryOperator(const boost::shared_ptr<TractionBoundaryOperatorParameters> & params);

        virtual ~TractionBoundaryOperator() { }

        void setTraction(AMP::LinearAlgebra::Vector::shared_ptr vec) {
          d_traction = mySubsetVector(vec, d_inputVar);
        }

        void apply( AMP::LinearAlgebra::Vector::const_shared_ptr, AMP::LinearAlgebra::Vector::const_shared_ptr,
            AMP::LinearAlgebra::Vector::shared_ptr r, const double, const double);

        void addRHScorrection(AMP::LinearAlgebra::Vector::shared_ptr rhs);

      protected :

        AMP::LinearAlgebra::Vector::shared_ptr mySubsetVector(AMP::LinearAlgebra::Vector::shared_ptr vec, 
            AMP::LinearAlgebra::Variable::shared_ptr var);

        void computeCorrection();

        void createCurrentLibMeshElement();

        void destroyCurrentLibMeshElement();

        std::vector<AMP::Mesh::MeshElement> d_currNodes;
        ::Elem* d_currElemPtr;
        AMP::LinearAlgebra::Variable::shared_ptr d_inputVar;
        AMP::LinearAlgebra::Variable::shared_ptr d_outputVar;
        AMP::LinearAlgebra::Vector::shared_ptr d_traction;
        AMP::LinearAlgebra::Vector::shared_ptr d_correction;
        bool d_residualMode;
        short int d_boundaryId;
    };

  }
}

#endif



