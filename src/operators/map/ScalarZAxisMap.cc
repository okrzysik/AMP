#include "operators/map/ScalarZAxisMap.h"
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
    AMP::Mesh::DOFMap::shared_ptr d = d_MeshAdapter->getDOFMap ( getInputVariable() );
    BNIterator cur = d_BeginNode;
    while ( cur != d_EndNode )
    {
        double val = v->getValueByGlobalID ( d->getGlobalID ( cur->globalID() , 0 ) );
        addTo1DMap ( cur->z() , val );
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
    const double TOL = 1e-12;
    AMP::Mesh::DOFMap::shared_ptr dof = d_MeshAdapter->getDOFMap ( getInputVariable() );
    BNIterator cur = d_BeginNode;
    while ( cur != d_EndNode ) {
        // Check the endpoints
        double pos = cur->z();
        if ( fabs(pos-z0) <= TOL ) {
            // We are within TOL of the first point
            vec->setValueByGlobalID( dof->getGlobalID(cur->globalID(),0), v0 );
            cur++;
            continue;
        } else if ( fabs(pos-z1) <= TOL ) {
            // We are within TOL of the last point
            vec->setValueByGlobalID( dof->getGlobalID(cur->globalID(),0), v1 );
            cur++;
            continue;
        } else if ( pos<z0 || pos>z1 ) {
            // We are outside the bounds of the map
            cur++;
            continue;
        } 

        // Find the first point > the current position
        ub = d_1DMultiMap.upper_bound( cur->z() );
        if ( ub == d_1DMultiMap.end() )
            ub--;
        lb = ub--;

        // Perform linear interpolation
        double lo = lb->first;
        double hi = ub->first;
        double wt = (pos - lo) / (hi - lo);
        double ans = (1.-wt) * lb->second + wt * ub->second;
        //double ans = ( (hi-pos) * lb->second + (pos-lo) * ub->second )/ (hi-lo);
        vec->setValueByGlobalID ( dof->getGlobalID ( cur->globalID() , 0 ) , ans );

        cur++;
    }
}


} // Operator namespace
} // AMP namespace


