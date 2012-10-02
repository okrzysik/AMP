#include "operators/map/SubchannelToCladGPMap.h"
#include "discretization/DOF_Manager.h"
#include "utils/ProfilerApp.h"


namespace AMP {
namespace Operator {




/************************************************************************
*  Default constructor                                                  *
************************************************************************/
SubchannelToCladGPMap::SubchannelToCladGPMap ( const boost::shared_ptr<AMP::Operator::OperatorParameters> &p )
    : SubchannelToCladMap ( p )
{
}


/************************************************************************
*  De-constructor                                                       *
************************************************************************/
SubchannelToCladGPMap::~SubchannelToCladGPMap ()
{
}


/************************************************************************
*  Check if the map type is "SubchannelToCladMap"                       *
************************************************************************/
bool SubchannelToCladGPMap::validMapType ( const std::string &t )
{
    if ( t == "SubchannelToCladGPMap" )
        return true;
    return false;
}


/************************************************************************
*  Fill the return vector for the given subchannel                      *
************************************************************************/    
void SubchannelToCladGPMap::fillReturnVector( AMP::LinearAlgebra::Vector::shared_ptr vec, double range[4], 
    const std::vector<AMP::Mesh::MeshElementID>& ids, const std::vector<double>& z, const std::vector<double>& f )
{
    AMP_ERROR("Not finished\n");
}


} // Operator namespace
} // AMP namespace


