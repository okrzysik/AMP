#if 1
#include "operators/map/CladToSubchannelMap.h"
#include "discretization/DOF_Manager.h"
#include "utils/PIO.h"
#include "utils/ProfilerApp.h"


namespace AMP {
namespace Operator {


/************************************************************************
*  Default constructor                                                  *
************************************************************************/
CladToSubchannelMap::CladToSubchannelMap ( const boost::shared_ptr<AMP::Operator::OperatorParameters> &p )
    : AsyncMapOperator ( p )
{
    boost::shared_ptr<CladToSubchannelMapParameters>  params = boost::dynamic_pointer_cast<CladToSubchannelMapParameters> ( p );
    AMP_ASSERT( params );

    int DofsPerObj = params->d_db->getInteger ( "DOFsPerObject" );
    AMP_INSIST(DofsPerObj==1,"CladToSubchannelMap is currently only designed for 1 DOF per node");

    // Get the iterator over the clad surface
    AMP::Mesh::MeshIterator cladIterator;
    if ( d_mesh1.get() != NULL )
        cladIterator = d_mesh1->getBoundaryIDIterator( AMP::Mesh::Vertex, params->d_BoundaryID1, 0 );
    
    // Get the x-y-z grid for the subchannel mesh
    fillSubchannelGrid();
    
}


/************************************************************************
*  De-constructor                                                       *
************************************************************************/
CladToSubchannelMap::~CladToSubchannelMap ()
{
}


/************************************************************************
*  Check if the map type is "CladToSubchannelMap"                       *
************************************************************************/
bool CladToSubchannelMap::validMapType ( const std::string &t )
{
    if ( t == "CladToSubchannelMap" )
        return true;
    return false;
}


/************************************************************************
*  Function to fill the grid of the subchannel for all processors       *
************************************************************************/
void CladToSubchannelMap::fillSubchannelGrid()
{
    std::set<double> x, y, z;
    if ( d_mesh2.get() != NULL ) {
        AMP::Mesh::MeshIterator it = d_mesh2->getIterator( AMP::Mesh::Vertex, 0 );
        for (size_t i=0; i<it.size(); i++) {
            std::vector<double> coord = it->coord();
            AMP_ASSERT(coord.size()==3);
            x.insert( coord[0] );
            y.insert( coord[1] );
            z.insert( coord[2] );
            ++it;
        }
    }
    d_MapComm.setGather(x);
    d_MapComm.setGather(y);
    d_MapComm.setGather(z);
    double last = 1e300;
    for (std::set<double>::iterator it=x.begin(); it!=x.end(); ++it) {
        if ( Utilities::approx_equal(last,*it,1e-12) )
            x.erase(it);
        else
            last = *it;
    }
    for (std::set<double>::iterator it=y.begin(); it!=y.end(); ++it) {
        if ( Utilities::approx_equal(last,*it,1e-12) )
            y.erase(it);
        else
            last = *it;
    }
    for (std::set<double>::iterator it=z.begin(); it!=z.end(); ++it) {
        if ( Utilities::approx_equal(last,*it,1e-12) )
            z.erase(it);
        else
            last = *it;
    }
    d_x = std::vector<double>(x.begin(),x.end());
    d_y = std::vector<double>(y.begin(),y.end());
    d_z = std::vector<double>(z.begin(),z.end());
    if ( d_mesh2.get() != NULL ) 
        AMP_ASSERT(((d_x.size()-1)*(d_y.size()-1)*(d_z.size()-1))==d_mesh2->numGlobalElements(AMP::Mesh::Volume));

}


/************************************************************************
*  Start the communication                                              *
************************************************************************/
void CladToSubchannelMap::applyStart(const AMP::LinearAlgebra::Vector::shared_ptr &f,
    const AMP::LinearAlgebra::Vector::shared_ptr &u, AMP::LinearAlgebra::Vector::shared_ptr &r,
    const double a, const double b)
{
    AMP_ERROR("Not finished");
}


/************************************************************************
*  Finish the communication                                             *
************************************************************************/
void CladToSubchannelMap::applyFinish(const AMP::LinearAlgebra::Vector::shared_ptr &f,
    const AMP::LinearAlgebra::Vector::shared_ptr &u, AMP::LinearAlgebra::Vector::shared_ptr &r,
    const double a, const double b)
{
    AMP_ERROR("Not finished");
}


/************************************************************************
*  Set the vector                                                       *
************************************************************************/
void  CladToSubchannelMap::setVector( AMP::LinearAlgebra::Vector::shared_ptr &result )
{
    AMP_ERROR("Not finished");
}

} // Operator namespace
} // AMP namespace

#endif
