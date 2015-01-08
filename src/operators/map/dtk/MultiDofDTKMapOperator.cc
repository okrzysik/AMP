#include <operators/map/dtk/MultiDofDTKMapOperator.h>

namespace AMP {
namespace Operator {

//---------------------------------------------------------------------------//
// Constructor
MultiDofDTKMapOperator::MultiDofDTKMapOperator( const AMP::shared_ptr<OperatorParameters>& params )
{
    // Get the operator parameters.
    AMP::shared_ptr<MultiDofDTKMapOperatorParameters> multiDofDTKMapOpParams =
        AMP::dynamic_pointer_cast<MultiDofDTKMapOperatorParameters>(params);
    AMP_ASSERT( multiDofDTKMapOpParams );

    AMP::Mesh::Mesh::shared_ptr mesh1 = multiDofDTKMapOpParams->d_Mesh1;
    AMP::Mesh::Mesh::shared_ptr mesh2 = multiDofDTKMapOpParams->d_Mesh2;
    int boundaryID1 = multiDofDTKMapOpParams->d_BoundaryID1;
    int boundaryID2 = multiDofDTKMapOpParams->d_BoundaryID2;
    std::string variable1 = multiDofDTKMapOpParams->d_Variable1;
    std::string variable2 = multiDofDTKMapOpParams->d_Variable2;
    std::size_t strideOffset1 = multiDofDTKMapOpParams->d_StrideOffset1;
    std::size_t strideOffset2 = multiDofDTKMapOpParams->d_StrideOffset2;
    std::size_t strideLength1 = multiDofDTKMapOpParams->d_StrideLength1;
    std::size_t strideLength2 = multiDofDTKMapOpParams->d_StrideLength2;
    AMP::LinearAlgebra::Vector::shared_ptr sourceVector = multiDofDTKMapOpParams->d_SourceVector;
    AMP::LinearAlgebra::Vector::shared_ptr targetVector = multiDofDTKMapOpParams->d_TargetVector;
    
    AMP::Mesh::Mesh::shared_ptr boundaryMesh1 = 
        mesh1->Subset(mesh1->getBoundaryIDIterator(AMP::Mesh::Volume, boundaryID1));
    AMP::Mesh::Mesh::shared_ptr boundaryMesh2 = 
        mesh2->Subset(mesh2->getBoundaryIDIterator(AMP::Mesh::Volume, boundaryID2));

    AMP::shared_ptr<AMP::Database> nullDatabase;
    // Build map 1 -> 2
    d_SourceVectorMap12 = sourceVector
            ->select(AMP::LinearAlgebra::VS_Mesh(boundaryMesh1)                 , "var")
            ->select(AMP::LinearAlgebra::VS_ByVariableName(variable1)           , "var")
            ->select(AMP::LinearAlgebra::VS_Stride(strideOffset1, strideLength1), "var");
    d_TargetVectorMap12 = targetVector
            ->select(AMP::LinearAlgebra::VS_Mesh(boundaryMesh2)                 , "var")
            ->select(AMP::LinearAlgebra::VS_ByVariableName(variable1)           , "var")
            ->select(AMP::LinearAlgebra::VS_Stride(strideOffset1, strideLength1), "var");
    AMP::shared_ptr<AMP::Operator::DTKMapOperatorParameters> map12Params(new AMP::Operator::DTKMapOperatorParameters(nullDatabase));
    map12Params->d_domain_mesh = boundaryMesh1;
    map12Params->d_range_mesh  = boundaryMesh2;
    map12Params->d_domain_dofs = d_SourceVectorMap12->getDOFManager();
    map12Params->d_range_dofs  = d_TargetVectorMap12->getDOFManager();
    d_Map12 = AMP::shared_ptr<AMP::Operator::DTKMapOperator>(new AMP::Operator::DTKMapOperator(map12Params));

    // Build map 2 -> 1
    d_SourceVectorMap21 = sourceVector
            ->select(AMP::LinearAlgebra::VS_Mesh(boundaryMesh2)                 , "var")
            ->select(AMP::LinearAlgebra::VS_ByVariableName(variable2)           , "var")
            ->select(AMP::LinearAlgebra::VS_Stride(strideOffset2, strideLength2), "var");
    d_TargetVectorMap21 = targetVector
            ->select(AMP::LinearAlgebra::VS_Mesh(boundaryMesh1)                 , "var")
            ->select(AMP::LinearAlgebra::VS_ByVariableName(variable2)           , "var")
            ->select(AMP::LinearAlgebra::VS_Stride(strideOffset2, strideLength2), "var");
    AMP::shared_ptr<AMP::Operator::DTKMapOperatorParameters> map21Params(new AMP::Operator::DTKMapOperatorParameters(nullDatabase));
    map21Params->d_domain_mesh = boundaryMesh2;
    map21Params->d_range_mesh  = boundaryMesh1;
    map21Params->d_domain_dofs = d_SourceVectorMap21->getDOFManager();
    map21Params->d_range_dofs  = d_TargetVectorMap21->getDOFManager();
    d_Map21 = AMP::shared_ptr<AMP::Operator::DTKMapOperator>(new AMP::Operator::DTKMapOperator(map21Params));

}

void 
MultiDofDTKMapOperator::
apply( AMP::LinearAlgebra::Vector::const_shared_ptr f, 
			 AMP::LinearAlgebra::Vector::const_shared_ptr u, 
			 AMP::LinearAlgebra::Vector::shared_ptr       r,
			 const double                                 a, 
			 const double                                 b )
{
    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    // QUESTION:  should we apply on u rather than on d_SourceVectorMapXY ?
    //            in that case we would have to perform select again
    d_Map12->apply(nullVec, d_SourceVectorMap12, d_TargetVectorMap12);
    d_Map21->apply(nullVec, d_SourceVectorMap21, d_TargetVectorMap21);
}


}
}


