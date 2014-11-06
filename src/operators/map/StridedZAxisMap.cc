#include "operators/map/StridedZAxisMap.h"
#include "discretization/DOF_Manager.h"
#include "utils/PIO.h"
#include "utils/ProfilerApp.h"


namespace AMP {
namespace Operator {


/************************************************************************
*  Default constructor                                                  *
************************************************************************/
StridedZAxisMap::StridedZAxisMap ( const AMP::shared_ptr<AMP::Operator::OperatorParameters> &p )
    : ScalarZAxisMap ( p )
{
    AMP::shared_ptr <Map3to1to3Parameters>  params = AMP::dynamic_pointer_cast<Map3to1to3Parameters> ( p );
    AMP_ASSERT ( params );

    d_inpDofs   = params->d_db->getIntegerWithDefault ( "InputDOFsPerObject", 1 );
    d_inpStride = params->d_db->getIntegerWithDefault ( "InputStride", 0 );
    d_outDofs   = params->d_db->getIntegerWithDefault ( "OutputDOFsPerObject", 1 );
    d_outStride = params->d_db->getIntegerWithDefault ( "OutputStride", 0 );
}


/************************************************************************
*  De-constructor                                                       *
************************************************************************/
StridedZAxisMap::~StridedZAxisMap ()
{
}

/************************************************************************
*  Check if the map type is "StridedZAxis"                              *
************************************************************************/
bool StridedZAxisMap::validMapType ( const std::string &t )
{
    if ( t == "StridedZAxis" )
        return true;
    return false;
}

void StridedZAxisMap::apply(AMP::LinearAlgebra::Vector::const_shared_ptr f,
        AMP::LinearAlgebra::Vector::const_shared_ptr u, 
        AMP::LinearAlgebra::Vector::shared_ptr r,
        const double a, const double b)
{

    AMP::LinearAlgebra::Variable::shared_ptr inpVar = getInputVariable();
    AMP::LinearAlgebra::Vector::const_shared_ptr  inpPhysics = u->constSubsetVectorForVariable(inpVar);
    AMP::LinearAlgebra::Vector::const_shared_ptr  inpStridedPhysics = inpPhysics->constSelect( AMP::LinearAlgebra::VS_Stride( d_inpStride, d_inpDofs) , inpVar->getName() );

    AMP::Operator::AsyncMapOperator::apply(f, inpStridedPhysics,  r, a, b);

}

void  StridedZAxisMap::setVector ( AMP::LinearAlgebra::Vector::shared_ptr result )
{
    AMP::LinearAlgebra::Variable::shared_ptr outVar = getOutputVariable();
    AMP::LinearAlgebra::Vector::shared_ptr  outPhysics = result->subsetVectorForVariable(outVar);
    AMP::LinearAlgebra::Vector::shared_ptr  outStridedPhysics = outPhysics->select( AMP::LinearAlgebra::VS_Stride( d_outStride, d_outDofs) , outVar->getName() );
 
    AMP::Operator::Map3to1to3::setVector(outStridedPhysics);
}

} // Operator namespace
} // AMP namespace


