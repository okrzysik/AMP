#include "AMP/operators/map/dtk/MultiDofDTKMapOperator.h"

namespace AMP::Operator {

//---------------------------------------------------------------------------//
// Constructor
MultiDofDTKMapOperator::MultiDofDTKMapOperator( std::shared_ptr<const OperatorParameters> params )
{
    // Get the operator parameters.
    auto multiDofDTKMapOpParams =
        std::dynamic_pointer_cast<const MultiDofDTKMapOperatorParameters>( params );
    AMP_ASSERT( multiDofDTKMapOpParams );
    d_multiDofDTKMapOpParams =
        std::dynamic_pointer_cast<const MultiDofDTKMapOperatorParameters>( params );

    auto mesh1            = multiDofDTKMapOpParams->d_Mesh1;
    auto mesh2            = multiDofDTKMapOpParams->d_Mesh2;
    int boundaryID1       = multiDofDTKMapOpParams->d_BoundaryID1;
    int boundaryID2       = multiDofDTKMapOpParams->d_BoundaryID2;
    std::string variable1 = multiDofDTKMapOpParams->d_Variable1;
    std::string variable2 = multiDofDTKMapOpParams->d_Variable2;
    size_t strideOffset1  = multiDofDTKMapOpParams->d_StrideOffset1;
    size_t strideOffset2  = multiDofDTKMapOpParams->d_StrideOffset2;
    size_t strideLength1  = multiDofDTKMapOpParams->d_StrideLength1;
    size_t strideLength2  = multiDofDTKMapOpParams->d_StrideLength2;
    auto sourceVector     = multiDofDTKMapOpParams->d_SourceVector;
    auto targetVector     = multiDofDTKMapOpParams->d_TargetVector;

    std::shared_ptr<AMP::Mesh::Mesh> boundaryMesh1_vol, boundaryMesh1_ver;
    std::shared_ptr<AMP::Mesh::Mesh> boundaryMesh2_vol, boundaryMesh2_ver;
    std::shared_ptr<AMP::Discretization::DOFManager> sourceDofManager12, sourceDofManager21;
    std::shared_ptr<AMP::Discretization::DOFManager> targetDofManager12, targetDofManager21;

    if ( mesh1 ) {
        boundaryMesh1_vol =
            mesh1->Subset( mesh1->getBoundaryIDIterator( AMP::Mesh::GeomType::Cell, boundaryID1 ) );
        boundaryMesh1_ver = mesh1->Subset(
            mesh1->getBoundaryIDIterator( AMP::Mesh::GeomType::Vertex, boundaryID1 ) );
        // Build map 1 -> 2
    }

    if ( boundaryMesh1_vol ) {
        d_SourceVectorMap12 =
            sourceVector->select( AMP::LinearAlgebra::VS_Mesh( boundaryMesh1_vol ) )
                ->select( AMP::LinearAlgebra::VS_ByVariableName( variable1 ) )
                ->select( AMP::LinearAlgebra::VS_Stride( strideOffset1, strideLength1 ) );
        sourceDofManager12 = d_SourceVectorMap12->getDOFManager();
    }

    int rank = multiDofDTKMapOpParams->d_globalComm.getRank();
    std::fstream fout;
    std::string fileName = "debug_" + std::to_string( rank );
    fout.open( fileName.c_str(), std::fstream::out );

    if ( mesh2 ) {
        AMP::Mesh::MeshIterator iterator_vol =
            mesh2->getBoundaryIDIterator( AMP::Mesh::GeomType::Cell, boundaryID2 );
        AMP::Mesh::MeshIterator iterator_ver =
            mesh2->getBoundaryIDIterator( AMP::Mesh::GeomType::Vertex, boundaryID2 );
        boundaryMesh2_vol = mesh2->Subset( iterator_vol );
        boundaryMesh2_ver = mesh2->Subset( iterator_ver );
    }

    if ( boundaryMesh2_ver ) {
        d_TargetVectorMap12 =
            targetVector->select( AMP::LinearAlgebra::VS_Mesh( boundaryMesh2_ver ) )
                ->select( AMP::LinearAlgebra::VS_ByVariableName( variable2 ) )
                ->select( AMP::LinearAlgebra::VS_Stride( strideOffset2, strideLength2 ) );
        targetDofManager12 = d_TargetVectorMap12->getDOFManager();
    }

    auto map12Params = std::make_shared<AMP::Operator::DTKMapOperatorParameters>( nullptr );
    map12Params->d_domain_mesh = boundaryMesh1_vol;
    map12Params->d_range_mesh  = boundaryMesh2_ver;
    map12Params->d_domain_dofs = sourceDofManager12;
    map12Params->d_range_dofs  = targetDofManager12;
    map12Params->d_globalComm  = multiDofDTKMapOpParams->d_globalComm;
    d_Map12                    = std::make_shared<AMP::Operator::DTKMapOperator>( map12Params );

    if ( boundaryMesh2_vol ) {
        // Build map 2 -> 1
        d_SourceVectorMap21 =
            sourceVector->select( AMP::LinearAlgebra::VS_Mesh( boundaryMesh2_vol ) )
                ->select( AMP::LinearAlgebra::VS_ByVariableName( variable2 ) )
                ->select( AMP::LinearAlgebra::VS_Stride( strideOffset2, strideLength2 ) );
        sourceDofManager21 = d_SourceVectorMap21->getDOFManager();
    }

    if ( boundaryMesh1_ver ) {
        d_TargetVectorMap21 =
            targetVector->select( AMP::LinearAlgebra::VS_Mesh( boundaryMesh1_ver ) )
                ->select( AMP::LinearAlgebra::VS_ByVariableName( variable1 ) )
                ->select( AMP::LinearAlgebra::VS_Stride( strideOffset1, strideLength1 ) );
        targetDofManager21 = d_TargetVectorMap21->getDOFManager();
    }

    auto map21Params = std::make_shared<AMP::Operator::DTKMapOperatorParameters>( nullptr );
    map21Params->d_domain_mesh = boundaryMesh2_vol;
    map21Params->d_range_mesh  = boundaryMesh1_ver;
    map21Params->d_domain_dofs = sourceDofManager21;
    map21Params->d_range_dofs  = targetDofManager21;
    map21Params->d_globalComm  = multiDofDTKMapOpParams->d_globalComm;
    d_Map21                    = std::make_shared<AMP::Operator::DTKMapOperator>( map21Params );
}


void MultiDofDTKMapOperator::apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                                    AMP::LinearAlgebra::Vector::shared_ptr r )
{

    auto mesh1            = d_multiDofDTKMapOpParams->d_Mesh1;
    auto mesh2            = d_multiDofDTKMapOpParams->d_Mesh2;
    int boundaryID1       = d_multiDofDTKMapOpParams->d_BoundaryID1;
    int boundaryID2       = d_multiDofDTKMapOpParams->d_BoundaryID2;
    std::string variable1 = d_multiDofDTKMapOpParams->d_Variable1;
    std::string variable2 = d_multiDofDTKMapOpParams->d_Variable2;
    size_t strideOffset1  = d_multiDofDTKMapOpParams->d_StrideOffset1;
    size_t strideOffset2  = d_multiDofDTKMapOpParams->d_StrideOffset2;
    size_t strideLength1  = d_multiDofDTKMapOpParams->d_StrideLength1;
    size_t strideLength2  = d_multiDofDTKMapOpParams->d_StrideLength2;

    std::shared_ptr<AMP::Mesh::Mesh> boundaryMesh1;
    std::shared_ptr<AMP::Mesh::Mesh> boundaryMesh2;
    if ( mesh1 ) {
        boundaryMesh1 =
            mesh1->Subset( mesh1->getBoundaryIDIterator( AMP::Mesh::GeomType::Cell, boundaryID1 ) );
    }
    if ( mesh2 ) {
        boundaryMesh2 =
            mesh2->Subset( mesh2->getBoundaryIDIterator( AMP::Mesh::GeomType::Cell, boundaryID2 ) );
    }

    if ( boundaryMesh1 ) {
        // Build map 1 -> 2
        d_SourceVectorMap12 =
            u->select( AMP::LinearAlgebra::VS_Mesh( boundaryMesh1 ) )
                ->select( AMP::LinearAlgebra::VS_ByVariableName( variable1 ) )
                ->select( AMP::LinearAlgebra::VS_Stride( strideOffset1, strideLength1 ) );
    }
    if ( boundaryMesh2 ) {
        // Build map 2 -> 1
        d_SourceVectorMap21 =
            u->select( AMP::LinearAlgebra::VS_Mesh( boundaryMesh2 ) )
                ->select( AMP::LinearAlgebra::VS_ByVariableName( variable2 ) )
                ->select( AMP::LinearAlgebra::VS_Stride( strideOffset2, strideLength2 ) );
    }

    // QUESTION:  should we apply on u rather than on d_SourceVectorMapXY ?
    //            in that case we would have to perform select again
    d_Map12->apply( d_SourceVectorMap12, d_TargetVectorMap12 );
    d_Map21->apply( d_SourceVectorMap21, d_TargetVectorMap21 );

    d_multiDofDTKMapOpParams->d_TargetVector->makeConsistent(
        AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );
}
} // namespace AMP::Operator
