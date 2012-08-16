#include "operators/map/ScalarZAxisMap.h"
#include "discretization/DOF_Manager.h"
#include "utils/PIO.h"
#include "utils/ProfilerApp.h"


namespace AMP {
namespace Operator {


/************************************************************************
*  Default constructor                                                  *
************************************************************************/
ScalarZAxisMap::ScalarZAxisMap ( const boost::shared_ptr<AMP::Operator::OperatorParameters> &p )
    : Map3to1to3 ( p )
{
    boost::shared_ptr <Map3to1to3Parameters>  params = boost::dynamic_pointer_cast<Map3to1to3Parameters> ( p );
    AMP_ASSERT ( params );

    int DofsPerObj = params->d_db->getInteger ( "DOFsPerObject" );
    AMP_INSIST(DofsPerObj==1,"ScalarZAxis is currently only designed for 1 DOF per node");

    // Create the element iterators
    if ( d_mesh1.get() != NULL ) {
        d_srcIterator1 = d_mesh1->getBoundaryIDIterator( AMP::Mesh::Vertex, params->d_BoundaryID1, 0 );
        d_dstIterator1 = d_mesh1->getBoundaryIDIterator( AMP::Mesh::Vertex, params->d_BoundaryID1, 0 );
    }
    if ( d_mesh2.get() != NULL ) {
        d_srcIterator2 = d_mesh2->getBoundaryIDIterator( AMP::Mesh::Vertex, params->d_BoundaryID2, 0 );
        d_dstIterator2 = d_mesh2->getBoundaryIDIterator( AMP::Mesh::Vertex, params->d_BoundaryID2, 0 );
    }
}


/************************************************************************
*  De-constructor                                                       *
************************************************************************/
ScalarZAxisMap::~ScalarZAxisMap ()
{
}


/************************************************************************
*  Check if the map type is "ScalarZAxis"                               *
************************************************************************/
bool ScalarZAxisMap::validMapType ( const std::string &t )
{
    if ( t == "ScalarZAxis" )
        return true;
    return false;
}


/************************************************************************
*  buildMap                                                             *
*  This function constructs the map from the given vector.              *
*  It loops through all values in "cur", storing the value iv using     *
*  the z-position as the 1D key.                                        *
************************************************************************/
std::multimap<double,double>  ScalarZAxisMap::buildMap( AMP::LinearAlgebra::Vector::const_shared_ptr vec, 
    const AMP::Mesh::Mesh::shared_ptr, const AMP::Mesh::MeshIterator &iterator )
{
    PROFILE_START("buildMap");
    std::multimap<double,double> map;
    AMP::Discretization::DOFManager::shared_ptr  dof = vec->getDOFManager( );
    AMP::Mesh::MeshIterator cur = iterator.begin();
    AMP::Mesh::MeshIterator end = iterator.end();
    std::vector<size_t> ids;
    while ( cur != end ) {
        dof->getDOFs( cur->globalID(), ids );
        AMP_ASSERT(ids.size()==1);
        double val = vec->getValueByGlobalID ( ids[0] );
        std::vector<double> x = cur->coord();
        addTo1DMap( map, x[2], val );
        cur++;
    }
    PROFILE_STOP("buildMap");
    return map;
}


/************************************************************************
*  buildReturn                                                          *
************************************************************************/
void ScalarZAxisMap::buildReturn ( const AMP::LinearAlgebra::Vector::shared_ptr vec, const AMP::Mesh::Mesh::shared_ptr, 
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


