#include "AMP/operators/map/StridedZAxisMap.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/utils/PIO.h"
#include "ProfilerApp.h"


namespace AMP {
namespace Operator {


/************************************************************************
 *  Default constructor                                                  *
 ************************************************************************/
StridedZAxisMap::StridedZAxisMap( std::shared_ptr<const AMP::Operator::OperatorParameters> p )
    : ScalarZAxisMap( p )
{
    auto params = std::dynamic_pointer_cast<const Map3to1to3Parameters>( p );
    AMP_ASSERT( params );

    d_inpDofs   = params->d_db->getWithDefault( "InputDOFsPerObject", 1 );
    d_inpStride = params->d_db->getWithDefault( "InputStride", 0 );
    d_outDofs   = params->d_db->getWithDefault( "OutputDOFsPerObject", 1 );
    d_outStride = params->d_db->getWithDefault( "OutputStride", 0 );
}


/************************************************************************
 *  Destructor                                                           *
 ************************************************************************/
StridedZAxisMap::~StridedZAxisMap() = default;


/************************************************************************
 *  Check if the map type is "StridedZAxis"                              *
 ************************************************************************/
bool StridedZAxisMap::validMapType( const std::string &t ) { return t == "StridedZAxis"; }

void StridedZAxisMap::apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                             AMP::LinearAlgebra::Vector::shared_ptr r )
{

    auto inpVar = getInputVariable();
    auto inpVec = u->constSubsetVectorForVariable( inpVar );
    if ( d_inpDofs != 1 ) {
        auto strided = inpVec->constSelect( AMP::LinearAlgebra::VS_Stride( d_inpStride, d_inpDofs ),
                                            inpVar->getName() );
        AMP_ASSERT( strided );
        AMP::Operator::AsyncMapOperator::apply( strided, r );
    } else {
        AMP::Operator::AsyncMapOperator::apply( inpVec, r );
    }
}

void StridedZAxisMap::setVector( AMP::LinearAlgebra::Vector::shared_ptr result )
{
    auto outVar = getOutputVariable();
    auto outVec = result->subsetVectorForVariable( outVar );
    if ( d_outDofs != 1 ) {
        auto strided = outVec->select( AMP::LinearAlgebra::VS_Stride( d_outStride, d_outDofs ),
                                       outVar->getName() );
        AMP_ASSERT( strided );
        AMP::Operator::Map3to1to3::setVector( strided );
    } else {
        AMP::Operator::Map3to1to3::setVector( outVec );
    }
}


} // namespace Operator
} // namespace AMP
