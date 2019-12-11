
#ifndef included_AMP_DTK_MultiDofDTKMapOperator
#define included_AMP_DTK_MultiDofDTKMapOperator

#include "AMP/ampmesh/Mesh.h"
#include "AMP/operators/Operator.h"
#include "AMP/operators/map/dtk/DTKMapOperator.h"

#include <string>

namespace AMP {
namespace Operator {

class MultiDofDTKMapOperatorParameters : public OperatorParameters
{
public:
    // Constructor.
    explicit MultiDofDTKMapOperatorParameters( std::shared_ptr<AMP::Database> db )
        : OperatorParameters( db )
    { /* ... */
    }

    AMP::LinearAlgebra::Vector::const_shared_ptr d_SourceVector;
    AMP::LinearAlgebra::Vector::shared_ptr d_TargetVector;

    AMP_MPI d_globalComm;
    AMP::Mesh::Mesh::shared_ptr d_Mesh1;
    AMP::Mesh::Mesh::shared_ptr d_Mesh2;
    int d_BoundaryID1;
    int d_BoundaryID2;
    std::string d_Variable1;
    std::string d_Variable2;
    std::size_t d_StrideOffset1;
    std::size_t d_StrideOffset2;
    std::size_t d_StrideLength1;
    std::size_t d_StrideLength2;
};


class MultiDofDTKMapOperator : public Operator
{
public:
    explicit MultiDofDTKMapOperator( const std::shared_ptr<OperatorParameters> &params );

    void apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                AMP::LinearAlgebra::Vector::shared_ptr r );

private:
    AMP::LinearAlgebra::VS_Comm createCommSelect( AMP_MPI globalComm, bool createOnThisRank );

    std::shared_ptr<MultiDofDTKMapOperatorParameters> d_multiDofDTKMapOpParams;
    std::shared_ptr<AMP::Operator::DTKMapOperator> d_Map12;
    std::shared_ptr<AMP::Operator::DTKMapOperator> d_Map21;
    AMP::LinearAlgebra::Vector::const_shared_ptr d_SourceVectorMap12;
    AMP::LinearAlgebra::Vector::shared_ptr d_TargetVectorMap12;
    AMP::LinearAlgebra::Vector::const_shared_ptr d_SourceVectorMap21;
    AMP::LinearAlgebra::Vector::shared_ptr d_TargetVectorMap21;
};
} // namespace Operator
} // namespace AMP

#endif
