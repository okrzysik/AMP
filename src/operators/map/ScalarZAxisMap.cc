#include "operators/map/ScalarZAxisMap.h"
#include "discretization/DOF_Manager.h"
#include "utils/PIO.h"

namespace AMP {
namespace Operator {


/************************************************************************
*  Default constructor                                                  *
************************************************************************/
ScalarZAxisMap::ScalarZAxisMap ( const boost::shared_ptr<AMP::Operator::OperatorParameters> &p )
    : Map3to1to3 ( p )
{
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
void ScalarZAxisMap::buildMap ( const AMP::LinearAlgebra::Vector::shared_ptr v )
{
    AMP::Discretization::DOFManager::shared_ptr  dof = v->getDOFManager( );
    AMP::Mesh::MeshIterator cur = d_BoundaryNodeIterator.begin();
    AMP::Mesh::MeshIterator end = d_BoundaryNodeIterator.end();
    std::vector<unsigned int> ids;
    while ( cur != end )
    {
        dof->getDOFs( *cur, ids );
        AMP_ASSERT(ids.size()==1);
        double val = v->getValueByGlobalID ( ids[0] );
        std::vector<double> x = cur->coord();
        addTo1DMap ( x[2] , val );
        cur++;
    }
}


/************************************************************************
*  buildReturn                                                          *
************************************************************************/
void ScalarZAxisMap::buildReturn ( const AMP::LinearAlgebra::Vector::shared_ptr vec )
{

    // Get the endpoints of the map
    std::multimap<double,double>::iterator lb, ub;
    lb = d_1DMultiMap.begin();
    ub = d_1DMultiMap.end(); ub--;
    double z0 = (*lb).first;
    double z1 = (*ub).first;
    double v0 = (*lb).second;
    double v1 = (*ub).second;

    // Loop through the points in the output vector
    const double TOL = 1e-8;
    AMP::Discretization::DOFManager::shared_ptr  dof = vec->getDOFManager( );
    AMP::Mesh::MeshIterator cur = d_BoundaryNodeIterator.begin();
    AMP::Mesh::MeshIterator end = d_BoundaryNodeIterator.end();
    std::vector<unsigned int> ids;
    while ( cur != end ) {
        // Check the endpoints
        std::vector<double> x = cur->coord();
        double pos = x[2];
        dof->getDOFs( *cur, ids );
        AMP_ASSERT(ids.size()==1);
        if ( fabs(pos-z0) <= TOL ) {
            // We are within TOL of the first point
            vec->setValueByGlobalID( ids[0], v0 );
            cur++;
            continue;
        } else if ( fabs(pos-z1) <= TOL ) {
            // We are within TOL of the last point
            vec->setValueByGlobalID( ids[0], v1 );
            cur++;
            continue;
        } else if ( pos<z0 || pos>z1 ) {
            // We are outside the bounds of the map
            cur++;
            continue;
        } 

        // Find the first point > the current position
        ub = d_1DMultiMap.upper_bound( pos );
        if ( ub == d_1DMultiMap.end() )
            ub--;
        lb = ub--;

        // Perform linear interpolation
        double lo = lb->first;
        double hi = ub->first;
        double wt = (pos - lo) / (hi - lo);
        double ans = (1.-wt) * lb->second + wt * ub->second;
        vec->setValueByGlobalID ( ids[0], ans );

        cur++;
    }
}


} // Operator namespace
} // AMP namespace


