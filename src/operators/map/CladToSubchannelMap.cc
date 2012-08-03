#if 0
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
    boost::shared_ptr <Map3to1to3Parameters>  params = boost::dynamic_pointer_cast<Map3to1to3Parameters> ( p );
    AMP_ASSERT ( params );

    int DofsPerObj = params->d_db->getInteger ( "DOFsPerObject" );
    AMP_INSIST(DofsPerObj==1,"CladToSubchannelMap is currently only designed for 1 DOF per node");

    // Get the iterator over the clad surface
    AMP::Mesh::MeshIterator cladIterator;
    if ( d_mesh1.get() != NULL )
        cladIterator = d_mesh1->getBoundaryIDIterator( AMP::Mesh::Vertex, params->d_BoundaryID1, 0 );
    
    // Get the x-y grid for the subchannel mesh
    std::set<double> x, y, z;
    if ( d_mesh2.get() != NULL ) {
        AMP::Mesh::MeshIterator it = d_mesh2->getIterator( AMP::Mesh::Vertex, 0 );
        for (size_t i=0; i<it.size(); i++) {
            std::vector<double> coord = it->coord();
            AMP_ASSERT(coor.size()==3);
            x.insert( coord[0] );
            y.insert( coord[1] );
            z.insert( coord[2] );
            ++it;
        }
    }
    d_MapComm.setReduce(x);
    d_MapComm.setReduce(y);
    d_MapComm.setReduce(z);
    double last = 1e400;
    for (std::set<double>::iterator it=x.begin(); it!=x.end(); ++it) {
        if ( approx_equal(last,*it,1e-12) )
            x.erase(it);
        else
            last = *it;
    }
    for (std::set<double>::iterator it=y.begin(); it!=y.end(); ++it) {
        if ( approx_equal(last,*it,1e-12) )
            y.erase(it);
        else
            last = *it;
    }
    for (std::set<double>::iterator it=z.begin(); it!=z.end(); ++it) {
        if ( approx_equal(last,*it,1e-12) )
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
*  De-constructor                                                       *
************************************************************************/
CladToSubchannelMap::~CladToSubchannelMap ()
{
}


/************************************************************************
*  Check if the map type is "ScalarZAxis"                               *
************************************************************************/
bool CladToSubchannelMap::validMapType ( const std::string &t )
{
    if ( t == "ScalarZAxis" )
        return true;
    return false;
}


/************************************************************************
*  buildReturn                                                          *
************************************************************************/
void CladToSubchannelMap::buildReturn ( const AMP::LinearAlgebra::Vector::shared_ptr vec, const AMP::Mesh::Mesh::shared_ptr, 
    const AMP::Mesh::MeshIterator &iterator, const std::map<double,double> &map )
{
    PROFILE_START("buildReturn");

    // Get the endpoints of the map
    AMP_ASSERT(map.size()>1);
    std::map<double,double>::const_iterator lb, ub;
    lb = map.begin();
    ub = map.end(); ub--;
    double z0 = (*lb).first;
    double z1 = (*ub).first;
    double v0 = (*lb).second;
    double v1 = (*ub).second;

    // Loop through the points in the output vector
    const double TOL = 1e-8;
    AMP::Discretization::DOFManager::shared_ptr  DOFs = vec->getDOFManager( );
    AMP::Mesh::MeshIterator cur = iterator.begin();
    AMP::Mesh::MeshIterator end = iterator.end();
    std::vector<size_t> ids;
    size_t dof;
    double pos;
    while ( cur != end ) {

        // Get the current position and DOF
        std::vector<double> x = cur->coord();
        pos = x[2];
        DOFs->getDOFs( cur->globalID(), ids );
        AMP_ASSERT(ids.size()==1);
        dof = ids[0];

        // Check the endpoints
        if ( fabs(pos-z0) <= TOL ) {
            // We are within TOL of the first point
            vec->setValueByGlobalID( dof, v0 );
            cur++;
            continue;
        } else if ( fabs(pos-z1) <= TOL ) {
            // We are within TOL of the last point
            vec->setValueByGlobalID( dof, v1 );
            cur++;
            continue;
        } else if ( pos<z0 || pos>z1 ) {
            // We are outside the bounds of the map
            cur++;
            continue;
        } 

        // Find the first point > the current position
        ub = map.upper_bound( pos );
        if ( ub == map.end() )
           ub--;
        else if ( ub == map.begin() )
           ub++;
        lb = ub;
        lb--;

        // Perform linear interpolation
        double lo = lb->first;
        double hi = ub->first;
        AMP_ASSERT(pos>=lo&&pos<hi);
        double wt = (pos - lo) / (hi - lo);
        double ans = (1.-wt) * lb->second + wt * ub->second;
        vec->setValueByGlobalID ( dof, ans );

        cur++;
    }
    PROFILE_STOP("buildReturn");
}


} // Operator namespace
} // AMP namespace

#endif
