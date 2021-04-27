#include "AMP/operators/map/dtk/MultiDofDTKMapOperator.h"

namespace AMP {
namespace Operator {

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

    AMP::Mesh::Mesh::shared_ptr mesh1 = multiDofDTKMapOpParams->d_Mesh1;
    AMP::Mesh::Mesh::shared_ptr mesh2 = multiDofDTKMapOpParams->d_Mesh2;
    int boundaryID1                   = multiDofDTKMapOpParams->d_BoundaryID1;
    int boundaryID2                   = multiDofDTKMapOpParams->d_BoundaryID2;
    std::string variable1             = multiDofDTKMapOpParams->d_Variable1;
    std::string variable2             = multiDofDTKMapOpParams->d_Variable2;
    std::size_t strideOffset1         = multiDofDTKMapOpParams->d_StrideOffset1;
    std::size_t strideOffset2         = multiDofDTKMapOpParams->d_StrideOffset2;
    std::size_t strideLength1         = multiDofDTKMapOpParams->d_StrideLength1;
    std::size_t strideLength2         = multiDofDTKMapOpParams->d_StrideLength2;
    AMP::LinearAlgebra::Vector::const_shared_ptr sourceVector =
        multiDofDTKMapOpParams->d_SourceVector;
    AMP::LinearAlgebra::Vector::shared_ptr targetVector = multiDofDTKMapOpParams->d_TargetVector;

    AMP::Mesh::Mesh::shared_ptr boundaryMesh1_vol, boundaryMesh1_ver;
    AMP::Mesh::Mesh::shared_ptr boundaryMesh2_vol, boundaryMesh2_ver;
    std::shared_ptr<AMP::Discretization::DOFManager> sourceDofManager12, sourceDofManager21;
    std::shared_ptr<AMP::Discretization::DOFManager> targetDofManager12, targetDofManager21;
    std::shared_ptr<AMP::Database> nullDatabase;

    if ( mesh1 ) {
        boundaryMesh1_vol = mesh1->Subset(
            mesh1->getBoundaryIDIterator( AMP::Mesh::GeomType::Volume, boundaryID1 ) );
        boundaryMesh1_ver = mesh1->Subset(
            mesh1->getBoundaryIDIterator( AMP::Mesh::GeomType::Vertex, boundaryID1 ) );
        // Build map 1 -> 2
    }
    AMP::LinearAlgebra::VS_Comm bndMesh1VolCommSelect =
        createCommSelect( multiDofDTKMapOpParams->d_globalComm, ( boundaryMesh1_vol != nullptr ) );
    AMP::LinearAlgebra::VS_Comm bndMesh1VerCommSelect =
        createCommSelect( multiDofDTKMapOpParams->d_globalComm, ( boundaryMesh1_ver != nullptr ) );
    AMP::LinearAlgebra::Vector::const_shared_ptr commSubsetSourceVec =
        sourceVector->constSelect( bndMesh1VolCommSelect, variable1 );

    if ( boundaryMesh1_vol ) {
        d_SourceVectorMap12 =
            commSubsetSourceVec
                ->constSelect( AMP::LinearAlgebra::VS_Mesh( boundaryMesh1_vol ), "var" )
                ->constSelect( AMP::LinearAlgebra::VS_ByVariableName( variable1 ), "var" )
                ->constSelect( AMP::LinearAlgebra::VS_Stride( strideOffset1, strideLength1 ),
                               "var" );
        sourceDofManager12 = d_SourceVectorMap12->getDOFManager();
    }

    int rank = multiDofDTKMapOpParams->d_globalComm.getRank();
    std::fstream fout;
    std::string fileName = "debug_" + AMP::Utilities::intToString( rank );
    fout.open( fileName.c_str(), std::fstream::out );

    if ( mesh2 ) {
        AMP::Mesh::MeshIterator iterator_vol =
            mesh2->getBoundaryIDIterator( AMP::Mesh::GeomType::Volume, boundaryID2 );
        AMP::Mesh::MeshIterator iterator_ver =
            mesh2->getBoundaryIDIterator( AMP::Mesh::GeomType::Vertex, boundaryID2 );
        boundaryMesh2_vol = mesh2->Subset( iterator_vol );
        boundaryMesh2_ver = mesh2->Subset( iterator_ver );
    }

    AMP::LinearAlgebra::VS_Comm bndMesh2VolCommSelect =
        createCommSelect( multiDofDTKMapOpParams->d_globalComm, ( boundaryMesh2_vol != nullptr ) );
    AMP::LinearAlgebra::VS_Comm bndMesh2VerCommSelect =
        createCommSelect( multiDofDTKMapOpParams->d_globalComm, ( boundaryMesh2_ver != nullptr ) );
    AMP::LinearAlgebra::Vector::shared_ptr commSubsetTargetVec =
        targetVector->select( bndMesh2VerCommSelect, variable2 );

    if ( boundaryMesh2_ver ) {
        d_TargetVectorMap12 =
            commSubsetTargetVec->select( AMP::LinearAlgebra::VS_Mesh( boundaryMesh2_ver ), "var" )
                ->select( AMP::LinearAlgebra::VS_ByVariableName( variable2 ), "var" )
                ->select( AMP::LinearAlgebra::VS_Stride( strideOffset2, strideLength2 ), "var" );
        targetDofManager12 = d_TargetVectorMap12->getDOFManager();
    }

    std::shared_ptr<AMP::Operator::DTKMapOperatorParameters> map12Params(
        new AMP::Operator::DTKMapOperatorParameters( nullDatabase ) );
    map12Params->d_domain_mesh = boundaryMesh1_vol;
    map12Params->d_range_mesh  = boundaryMesh2_ver;
    map12Params->d_domain_dofs = sourceDofManager12;
    map12Params->d_range_dofs  = targetDofManager12;
    map12Params->d_globalComm  = multiDofDTKMapOpParams->d_globalComm;
    d_Map12                    = std::shared_ptr<AMP::Operator::DTKMapOperator>(
        new AMP::Operator::DTKMapOperator( map12Params ) );

    commSubsetSourceVec = sourceVector->constSelect( bndMesh2VolCommSelect, variable2 );
    if ( boundaryMesh2_vol ) {
        // Build map 2 -> 1
        d_SourceVectorMap21 =
            commSubsetSourceVec
                ->constSelect( AMP::LinearAlgebra::VS_Mesh( boundaryMesh2_vol ), "var" )
                ->constSelect( AMP::LinearAlgebra::VS_ByVariableName( variable2 ), "var" )
                ->constSelect( AMP::LinearAlgebra::VS_Stride( strideOffset2, strideLength2 ),
                               "var" );
        sourceDofManager21 = d_SourceVectorMap21->getDOFManager();
    }

    commSubsetTargetVec = targetVector->select( bndMesh1VerCommSelect, variable1 );
    if ( boundaryMesh1_ver ) {
        d_TargetVectorMap21 =
            commSubsetTargetVec->select( AMP::LinearAlgebra::VS_Mesh( boundaryMesh1_ver ), "var" )
                ->select( AMP::LinearAlgebra::VS_ByVariableName( variable1 ), "var" )
                ->select( AMP::LinearAlgebra::VS_Stride( strideOffset1, strideLength1 ), "var" );
        targetDofManager21 = d_TargetVectorMap21->getDOFManager();
    }

    std::shared_ptr<AMP::Operator::DTKMapOperatorParameters> map21Params(
        new AMP::Operator::DTKMapOperatorParameters( nullDatabase ) );
    map21Params->d_domain_mesh = boundaryMesh2_vol;
    map21Params->d_range_mesh  = boundaryMesh1_ver;
    map21Params->d_domain_dofs = sourceDofManager21;
    map21Params->d_range_dofs  = targetDofManager21;
    map21Params->d_globalComm  = multiDofDTKMapOpParams->d_globalComm;
    d_Map21                    = std::shared_ptr<AMP::Operator::DTKMapOperator>(
        new AMP::Operator::DTKMapOperator( map21Params ) );
}

AMP::LinearAlgebra::VS_Comm MultiDofDTKMapOperator::createCommSelect( AMP_MPI globalComm,
                                                                      bool createOnThisRank )
{
    int inComm = createOnThisRank ? 1 : 0;
    // Create a comm spanning the meshes
    AMP::LinearAlgebra::VS_Comm commSelect( AMP_MPI( globalComm.split( inComm ) ) );
    return commSelect;
}
void MultiDofDTKMapOperator::apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                                    AMP::LinearAlgebra::Vector::shared_ptr r )
{

    AMP::Mesh::Mesh::shared_ptr mesh1 = d_multiDofDTKMapOpParams->d_Mesh1;
    AMP::Mesh::Mesh::shared_ptr mesh2 = d_multiDofDTKMapOpParams->d_Mesh2;
    int boundaryID1                   = d_multiDofDTKMapOpParams->d_BoundaryID1;
    int boundaryID2                   = d_multiDofDTKMapOpParams->d_BoundaryID2;
    std::string variable1             = d_multiDofDTKMapOpParams->d_Variable1;
    std::string variable2             = d_multiDofDTKMapOpParams->d_Variable2;
    std::size_t strideOffset1         = d_multiDofDTKMapOpParams->d_StrideOffset1;
    std::size_t strideOffset2         = d_multiDofDTKMapOpParams->d_StrideOffset2;
    std::size_t strideLength1         = d_multiDofDTKMapOpParams->d_StrideLength1;
    std::size_t strideLength2         = d_multiDofDTKMapOpParams->d_StrideLength2;

    AMP::Mesh::Mesh::shared_ptr boundaryMesh1;
    AMP::Mesh::Mesh::shared_ptr boundaryMesh2;
    std::shared_ptr<AMP::Database> nullDatabase;
    if ( mesh1 ) {
        boundaryMesh1 = mesh1->Subset(
            mesh1->getBoundaryIDIterator( AMP::Mesh::GeomType::Volume, boundaryID1 ) );
    }
    AMP::LinearAlgebra::VS_Comm bndMesh1CommSelect =
        createCommSelect( d_multiDofDTKMapOpParams->d_globalComm, ( boundaryMesh1 != nullptr ) );
    AMP::LinearAlgebra::Vector::const_shared_ptr commSubsetUVec =
        u->constSelect( bndMesh1CommSelect, variable1 );

    if ( boundaryMesh1 ) {
        // Build map 1 -> 2
        d_SourceVectorMap12 =
            commSubsetUVec->constSelect( AMP::LinearAlgebra::VS_Mesh( boundaryMesh1 ), "var" )
                ->constSelect( AMP::LinearAlgebra::VS_ByVariableName( variable1 ), "var" )
                ->constSelect( AMP::LinearAlgebra::VS_Stride( strideOffset1, strideLength1 ),
                               "var" );
    }

    if ( mesh2 ) {
        boundaryMesh2 = mesh2->Subset(
            mesh2->getBoundaryIDIterator( AMP::Mesh::GeomType::Volume, boundaryID2 ) );
    }
    AMP::LinearAlgebra::VS_Comm bndMesh2CommSelect =
        createCommSelect( d_multiDofDTKMapOpParams->d_globalComm, ( boundaryMesh2 != nullptr ) );
    commSubsetUVec = u->constSelect( bndMesh2CommSelect, variable2 );

    if ( boundaryMesh2 ) {
        // Build map 2 -> 1
        d_SourceVectorMap21 =
            commSubsetUVec->constSelect( AMP::LinearAlgebra::VS_Mesh( boundaryMesh2 ), "var" )
                ->constSelect( AMP::LinearAlgebra::VS_ByVariableName( variable2 ), "var" )
                ->constSelect( AMP::LinearAlgebra::VS_Stride( strideOffset2, strideLength2 ),
                               "var" );
    }

    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    // QUESTION:  should we apply on u rather than on d_SourceVectorMapXY ?
    //            in that case we would have to perform select again
    d_Map12->apply( d_SourceVectorMap12, d_TargetVectorMap12 );
    d_Map21->apply( d_SourceVectorMap21, d_TargetVectorMap21 );

    d_multiDofDTKMapOpParams->d_TargetVector->makeConsistent(
        AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );
}
} // namespace Operator
} // namespace AMP
