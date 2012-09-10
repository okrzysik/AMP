#include "operators/map/ScalarN2GZAxisMap.h"
#include "discretization/DOF_Manager.h"
#include "utils/PIO.h"
#include "utils/ProfilerApp.h"
#include "discretization/simpleDOF_Manager.h"
#include "vectors/VectorBuilder.h"


/* Libmesh files */
#include "fe_type.h"
#include "fe_base.h"
#include "elem.h"
#include "quadrature.h"

#include "enum_order.h"
#include "enum_fe_family.h"
#include "enum_quadrature_type.h"
#include "auto_ptr.h"
#include "string_to_enum.h"

#include "face_quad4.h"
#include "node.h"


namespace AMP {
namespace Operator {


/************************************************************************
*  Default constructor                                                  *
************************************************************************/
ScalarN2GZAxisMap::ScalarN2GZAxisMap ( const boost::shared_ptr<AMP::Operator::OperatorParameters> &p )
    : Map3to1to3 ( p )
{
    boost::shared_ptr <Map3to1to3Parameters>  params = boost::dynamic_pointer_cast<Map3to1to3Parameters> ( p );
    AMP_ASSERT ( params );

    int DofsPerObj = params->d_db->getInteger ( "DOFsPerObject" );
    AMP_INSIST(DofsPerObj==4,"ScalarZAxis is currently only designed for 4 Gp per elem");

    // Create the element iterators
    if ( d_mesh1.get() != NULL ) {
        d_srcIterator1 = d_mesh1->getBoundaryIDIterator( AMP::Mesh::Vertex, params->d_BoundaryID1, 0 );
        d_dstIterator1 = d_mesh1->getBoundaryIDIterator( AMP::Mesh::Face, params->d_BoundaryID1, 0 );
    }
    if ( d_mesh2.get() != NULL ) {
        d_srcIterator2 = d_mesh2->getBoundaryIDIterator( AMP::Mesh::Vertex, params->d_BoundaryID2, 0 );
        d_dstIterator2 = d_mesh2->getBoundaryIDIterator( AMP::Mesh::Face, params->d_BoundaryID2, 0 );
    }
   
    AMP::Mesh::MeshIterator iterator = AMP::Mesh::Mesh::getIterator( AMP::Mesh::Union, d_dstIterator1 , d_dstIterator2 );
    libmeshElements.reinit( iterator );

    // Build the z-coordinates of the return gauss points
    d_z_coord1 = getGaussPoints( d_dstIterator1 );
    d_z_coord2 = getGaussPoints( d_dstIterator2 );
}


/************************************************************************
*  De-constructor                                                       *
************************************************************************/
ScalarN2GZAxisMap::~ScalarN2GZAxisMap ()
{
}

/************************************************************************
*  Check if the map type is "ScalarN2GZAxis"                               *
************************************************************************/
bool ScalarN2GZAxisMap::validMapType ( const std::string &t )
{
    if ( t == "ScalarN2GZAxis" )
        return true;
    return false;
}


/************************************************************************
*  buildMap                                                             *
*  This function constructs the map from the given vector.              *
*  It loops through all values in "cur", storing the value iv using     *
*  the z-position as the 1D key.                                        *
************************************************************************/
std::multimap<double,double>  ScalarN2GZAxisMap::buildMap( AMP::LinearAlgebra::Vector::const_shared_ptr vec, 
    const AMP::Mesh::Mesh::shared_ptr, const AMP::Mesh::MeshIterator &iterator )
{
    PROFILE_START("buildMap");
    std::multimap<double,double> map;
    AMP::Discretization::DOFManager::shared_ptr  dof = vec->getDOFManager( );
    AMP::Mesh::MeshIterator cur = iterator.begin();
    std::vector<size_t> ids;
    for (size_t i=0; i<cur.size(); i++) {
        dof->getDOFs( cur->globalID(), ids );
        AMP_ASSERT(ids.size()==1);
        double val = vec->getValueByGlobalID ( ids[0] );
        std::vector<double> x = cur->coord();
        addTo1DMap( map, x[2], val );
        ++cur;
    }
    PROFILE_STOP("buildMap");
    return map;
}


/************************************************************************
*  Function to build the z-coordinates of the gauss points              *
************************************************************************/
AMP::LinearAlgebra::Vector::const_shared_ptr ScalarN2GZAxisMap::getGaussPoints( const AMP::Mesh::MeshIterator& iterator )
{
    if ( iterator.size()==0 )
        return AMP::LinearAlgebra::Vector::const_shared_ptr();
    if ( iterator==d_dstIterator1 && d_z_coord1.get()!=NULL )
        return d_z_coord1;
    if ( iterator==d_dstIterator2 && d_z_coord2.get()!=NULL )
        return d_z_coord2;
    AMP::Discretization::DOFManager::shared_ptr GpDofMap = 
        AMP::Discretization::simpleDOFManager::create(iterator,4);
    AMP::LinearAlgebra::Variable::shared_ptr var(new AMP::LinearAlgebra::Variable("gauss_z"));
    AMP::LinearAlgebra::Vector::shared_ptr z_pos = AMP::LinearAlgebra::createVector(  GpDofMap, var );
    AMP::Mesh::MeshIterator cur = iterator.begin();
    std::vector<size_t> ids;
    for (size_t i=0; i<cur.size(); i++) {
        // Create the libmesh element
        libMeshEnums::Order feTypeOrder = Utility::string_to_enum<libMeshEnums::Order>("FIRST");
        libMeshEnums::FEFamily feFamily = Utility::string_to_enum<libMeshEnums::FEFamily>("LAGRANGE");
        boost::shared_ptr < ::FEType > d_feType ( new ::FEType(feTypeOrder, feFamily) );
        boost::shared_ptr < ::FEBase > d_fe ( (::FEBase::build(2, (*d_feType))).release() );
        libMeshEnums::Order qruleOrder = Utility::string_to_enum<libMeshEnums::Order>("SECOND");
        boost::shared_ptr < ::QBase > d_qrule ( (::QBase::build("QGAUSS", 2, qruleOrder)).release() );
        d_fe->attach_quadrature_rule( d_qrule.get() );
        d_fe->reinit ( libmeshElements.getElement( cur->globalID() ));
        // Get the current position and DOF
        std::vector<Point> coordinates = d_fe->get_xyz();
        GpDofMap->getDOFs( cur->globalID(), ids );
        for (unsigned int qp = 0; qp < ids.size(); qp++) {
            double pos = coordinates[qp](2);
            z_pos->setLocalValueByGlobalID( ids[qp], pos );
        }
        ++cur;
    }
    return z_pos;
}


/************************************************************************
*  buildReturn                                                          *
************************************************************************/
void ScalarN2GZAxisMap::buildReturn ( const AMP::LinearAlgebra::Vector::shared_ptr vec, const AMP::Mesh::Mesh::shared_ptr, 
    const AMP::Mesh::MeshIterator &iterator, const std::map<double,double> &map )
{
    if ( iterator.size()==0 )
        return;
    PROFILE_START("buildReturn");
    const double TOL = 1e-8;

    // Convert the map to a std::vector
    AMP_ASSERT(map.size()>1);
    std::vector<double> z(map.size()), f(map.size());
    size_t i=0;
    for (std::map<double,double>::const_iterator it=map.begin(); it!=map.end(); it++) {
        z[i] = it->first;
        f[i] = it->second;
        i++;
    }
    for (size_t i=1; i<z.size(); i++)
        AMP_ASSERT(z[i]>(z[i-1]+TOL));
    double z0 = z[0];
    double z1 = z[z.size()-1];
    double v0 = f[0];
    double v1 = f[z.size()-1];

    // Get the coordinates of the gauss points
    AMP::LinearAlgebra::Vector::const_shared_ptr z_pos = getGaussPoints( iterator );
    AMP_ASSERT(z_pos.get()!=NULL);

    // Get the DOF managers
    AMP::Discretization::DOFManager::shared_ptr  DOFs = vec->getDOFManager( );
    AMP::Discretization::DOFManager::shared_ptr  gaussDOFs = z_pos->getDOFManager( );

    // Loop through the points in the output vector
    AMP::Mesh::MeshIterator cur = iterator.begin();
    std::vector<size_t> id1, id2;
    std::vector<double> zi, fi;
    std::vector<bool> keep;
    for (i=0; i<cur.size(); i++) {

        // Get the DOFs
        DOFs->getDOFs( cur->globalID(), id1 );
        gaussDOFs->getDOFs( cur->globalID(), id2 );
        AMP_ASSERT(id1.size()==id2.size());
        size_t N = id1.size();
        zi.resize(N);
        fi.resize(N);
        keep.resize(N);

        // Get the coordinates of the gauss points
        z_pos->getLocalValuesByGlobalID(N,&id2[0],&zi[0]);

        bool keep_all = true;
        for (size_t j=0; j<N; j++) {
            // Check the endpoints
            keep[j] = true;
            if ( fabs(zi[j]-z0) <= TOL ) {
                // We are within TOL of the first point
                zi[j] = v0;
                continue;
            } else if ( fabs(zi[j]-z1) <= TOL ) {
                // We are within TOL of the last point
                zi[j] = v1;
                continue;
            } else if ( zi[j]<z0 || zi[j]>z1 ) {
                // We are outside the bounds of the map
                keep[j] = false;
                keep_all = false;
                continue;
            } 

            // Find the first point > the current position
            size_t k = AMP::Utilities::findfirst(z,zi[j]);
            if ( k==0 ) { k++; }
            if ( k==z.size() ) { k--; }

            // Perform linear interpolation
            double wt = (zi[j]-z[k-1])/(z[k]-z[k-1]);
            fi[j] = (1.0-wt)*f[k-1] + wt*f[k];
        }

        if ( keep_all ) {
            vec->setLocalValuesByGlobalID( N, &id1[0], &fi[0] );
        } else {
            for (size_t j=0; j<N; j++) {
                if ( keep[j] )
                    vec->setLocalValueByGlobalID( id1[j], fi[j] );
            }
        }

        ++cur;
    }
    PROFILE_STOP("buildReturn");
}


} // Operator namespace
} // AMP namespace


