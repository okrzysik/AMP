#include "operators/map/ScalarN2GZAxisMap.h"
#include "ProfilerApp.h"
#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"
#include "utils/PIO.h"
#include "vectors/VectorBuilder.h"


// Libmesh files
DISABLE_WARNINGS
#include "libmesh/auto_ptr.h"
#include "libmesh/elem.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/face_quad4.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_type.h"
#include "libmesh/node.h"
#include "libmesh/quadrature.h"
#include "libmesh/string_to_enum.h"
ENABLE_WARNINGS


namespace AMP {
namespace Operator {


/************************************************************************
*  Default constructor                                                  *
************************************************************************/
ScalarN2GZAxisMap::ScalarN2GZAxisMap( const AMP::shared_ptr<AMP::Operator::OperatorParameters> &p )
    : Map3to1to3( p )
{
    AMP::shared_ptr<Map3to1to3Parameters> params =
        AMP::dynamic_pointer_cast<Map3to1to3Parameters>( p );
    AMP_ASSERT( params );

    int DofsPerObj = params->d_db->getInteger( "DOFsPerObject" );
    AMP_INSIST( DofsPerObj == 4, "ScalarZAxis is currently only designed for 4 Gp per elem" );

    // Create the element iterators
    if ( d_mesh1.get() != nullptr ) {
        d_srcIterator1 =
            d_mesh1->getBoundaryIDIterator( AMP::Mesh::GeomType::Vertex, params->d_BoundaryID1, 0 );
        d_dstIterator1 =
            d_mesh1->getBoundaryIDIterator( AMP::Mesh::GeomType::Face, params->d_BoundaryID1, 0 );
    }
    if ( d_mesh2.get() != nullptr ) {
        d_srcIterator2 =
            d_mesh2->getBoundaryIDIterator( AMP::Mesh::GeomType::Vertex, params->d_BoundaryID2, 0 );
        d_dstIterator2 =
            d_mesh2->getBoundaryIDIterator( AMP::Mesh::GeomType::Face, params->d_BoundaryID2, 0 );
    }

    AMP::Mesh::MeshIterator iterator =
        AMP::Mesh::Mesh::getIterator( AMP::Mesh::SetOP::Union, d_dstIterator1, d_dstIterator2 );
    libmeshElements.reinit( iterator );

    // Build the z-coordinates of the return gauss points
    d_z_coord1 = getGaussPoints( d_dstIterator1 );
    d_z_coord2 = getGaussPoints( d_dstIterator2 );
}


/************************************************************************
*  De-constructor                                                       *
************************************************************************/
ScalarN2GZAxisMap::~ScalarN2GZAxisMap() {}

/************************************************************************
*  Check if the map type is "ScalarN2GZAxis"                               *
************************************************************************/
bool ScalarN2GZAxisMap::validMapType( const std::string &t )
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
std::multimap<double, double>
ScalarN2GZAxisMap::buildMap( AMP::LinearAlgebra::Vector::const_shared_ptr vec,
                             const AMP::Mesh::Mesh::shared_ptr,
                             const AMP::Mesh::MeshIterator &iterator )
{
    if ( iterator.size() == 0 )
        return std::multimap<double, double>();
    PROFILE_START( "buildMap" );
    std::multimap<double, double> map;
    AMP::Discretization::DOFManager::shared_ptr dof = vec->getDOFManager();
    size_t N                                        = iterator.size();
    std::vector<AMP::Mesh::MeshElementID> ids( N );
    std::vector<double> z( N, 0.0 );
    AMP::Mesh::MeshIterator it = iterator.begin();
    for ( size_t i = 0; i < N; ++i, ++it ) {
        ids[i] = it->globalID();
        z[i]   = it->coord( 2 );
    }
    std::vector<size_t> dofs( N );
    dof->getDOFs( ids, dofs );
    AMP_ASSERT( ids.size() == dofs.size() );
    std::vector<double> vals( N );
    vec->getValuesByGlobalID( N, &dofs[0], &vals[0] );
    addTo1DMap( map, z, vals );
    PROFILE_STOP( "buildMap" );
    return map;
}


/************************************************************************
*  Function to build the z-coordinates of the gauss points              *
************************************************************************/
AMP::LinearAlgebra::Vector::const_shared_ptr
ScalarN2GZAxisMap::getGaussPoints( const AMP::Mesh::MeshIterator &iterator )
{
    if ( iterator.size() == 0 )
        return AMP::LinearAlgebra::Vector::const_shared_ptr();
    if ( iterator == d_dstIterator1 && d_z_coord1.get() != nullptr )
        return d_z_coord1;
    if ( iterator == d_dstIterator2 && d_z_coord2.get() != nullptr )
        return d_z_coord2;
    PROFILE_START( "getGaussPoints" );
    AMP::Discretization::DOFManager::shared_ptr GpDofMap =
        AMP::Discretization::simpleDOFManager::create( iterator, 4 );
    AMP::LinearAlgebra::Variable::shared_ptr var( new AMP::LinearAlgebra::Variable( "gauss_z" ) );
    AMP::LinearAlgebra::Vector::shared_ptr z_pos =
        AMP::LinearAlgebra::createVector( GpDofMap, var );
    AMP::Mesh::MeshIterator cur = iterator.begin();
    std::vector<size_t> ids;
    for ( size_t i = 0; i < cur.size(); i++ ) {
        // Create the libmesh element
        libMeshEnums::Order feTypeOrder = Utility::string_to_enum<libMeshEnums::Order>( "FIRST" );
        libMeshEnums::FEFamily feFamily =
            Utility::string_to_enum<libMeshEnums::FEFamily>( "LAGRANGE" );
        AMP::shared_ptr<::FEType> d_feType( new ::FEType( feTypeOrder, feFamily ) );
        AMP::shared_ptr<::FEBase> d_fe( (::FEBase::build( 2, ( *d_feType ) ) ).release() );
        libMeshEnums::Order qruleOrder = Utility::string_to_enum<libMeshEnums::Order>( "SECOND" );
        AMP::shared_ptr<::QBase> d_qrule( (::QBase::build( "QGAUSS", 2, qruleOrder ) ).release() );
        d_fe->attach_quadrature_rule( d_qrule.get() );
        d_fe->reinit( libmeshElements.getElement( cur->globalID() ) );
        // Get the current position and DOF
        std::vector<Point> coordinates = d_fe->get_xyz();
        GpDofMap->getDOFs( cur->globalID(), ids );
        for ( unsigned int qp = 0; qp < ids.size(); qp++ ) {
            double pos = coordinates[qp]( 2 );
            z_pos->setLocalValueByGlobalID( ids[qp], pos );
        }
        ++cur;
    }
    PROFILE_STOP( "getGaussPoints" );
    return z_pos;
}


/************************************************************************
*  buildReturn                                                          *
************************************************************************/
void ScalarN2GZAxisMap::buildReturn( AMP::LinearAlgebra::Vector::shared_ptr vec,
                                     const AMP::Mesh::Mesh::shared_ptr,
                                     const AMP::Mesh::MeshIterator &iterator,
                                     const std::map<double, double> &map )
{
    if ( iterator.size() == 0 )
        return;
    PROFILE_START( "buildReturn" );
    const double TOL = 1e-8;

    // Convert the map to a std::vector
    AMP_ASSERT( map.size() > 1 );
    std::vector<double> z( map.size() ), f( map.size() );
    auto it_map = map.begin();
    for ( size_t i = 0; it_map != map.end(); ++it_map, ++i ) {
        z[i] = it_map->first;
        f[i] = it_map->second;
    }
    for ( size_t i = 1; i < z.size(); i++ )
        AMP_ASSERT( z[i] > ( z[i - 1] + TOL ) );
    double z0 = z[0];
    double z1 = z[z.size() - 1];

    // Get the coordinates of the gauss points
    AMP::LinearAlgebra::Vector::const_shared_ptr z_pos = getGaussPoints( iterator );
    AMP_ASSERT( z_pos.get() != nullptr );

    // Get the DOF managers
    AMP::Discretization::DOFManager::shared_ptr DOFs      = vec->getDOFManager();
    AMP::Discretization::DOFManager::shared_ptr gaussDOFs = z_pos->getDOFManager();

    // Loop through the points in the output vector
    size_t N0 = iterator.size();
    std::vector<size_t> dofs;
    std::vector<double> zi;
    dofs.reserve( N0 );
    zi.reserve( N0 );
    std::vector<size_t> id1, id2;
    std::vector<double> zi2;
    AMP::Mesh::MeshIterator it_mesh = iterator.begin();
    for ( size_t i = 0; i < N0; ++i, ++it_mesh ) {
        // Get the local DOFs
        DOFs->getDOFs( it_mesh->globalID(), id1 );
        gaussDOFs->getDOFs( it_mesh->globalID(), id2 );
        AMP_ASSERT( id1.size() == id2.size() );
        size_t N2 = id1.size();
        // Get the coordinates of the gauss points
        zi2.resize( N2 );
        z_pos->getLocalValuesByGlobalID( N2, &id2[0], &zi2[0] );
        for ( size_t j = 0; j < N2; j++ ) {
            if ( zi2[j] < z0 - TOL || zi2[j] > z1 + TOL ) {
                // We are outside the bounds of the map
                continue;
            }
            dofs.push_back( id1[j] );
            zi.push_back( zi2[j] );
        }
    }
    std::vector<double> fi( zi.size() );
    for ( size_t i = 0; i < zi.size(); i++ ) {
        // Find the first point > the current position
        size_t k = AMP::Utilities::findfirst( z, zi[i] );
        k        = std::max<size_t>( k, 1 );
        k        = std::min<size_t>( k, z.size() - 1 );
        // Perform linear interpolation
        double wt = ( zi[i] - z[k - 1] ) / ( z[k] - z[k - 1] );
        fi[i]     = ( 1.0 - wt ) * f[k - 1] + wt * f[k];
    }
    vec->setLocalValuesByGlobalID( dofs.size(), &dofs[0], &fi[0] );

    PROFILE_STOP( "buildReturn" );
}


} // Operator namespace
} // AMP namespace
