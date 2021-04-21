#include "AMP/ampmesh/Mesh.h"
#include "AMP/ampmesh/Geometry.h"
#include "AMP/ampmesh/MeshElementVectorIterator.h"
#include "AMP/ampmesh/MeshGeometry.h"
#include "AMP/ampmesh/MeshParameters.h"
#include "AMP/ampmesh/MeshUtilities.h"
#include "AMP/ampmesh/MultiMesh.h"
#include "AMP/ampmesh/SubsetMesh.h"
#include "AMP/utils/AMP_MPI.I"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/kdtree.h"
#ifdef USE_AMP_VECTORS
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"
#endif
#ifdef USE_AMP_DISCRETIZATION
#include "AMP/discretization/simpleDOF_Manager.h"
#endif

#include <cmath>


namespace AMP::Mesh {


static_assert( sizeof( MeshID ) == 8, "unexpected size for MeshID" );
static_assert( sizeof( ElementID ) == 8, "unexpected size for ElementID" );
static_assert( sizeof( MeshElementID ) == 16, "unexpected size for MeshElementID" );


static unsigned int nextLocalMeshID = 1;


/********************************************************
 * Constructors                                          *
 ********************************************************/
Mesh::Mesh( const std::shared_ptr<MeshParameters> &params_in )
{
    // Set the base properties
    AMP_ASSERT( sizeof( MeshElementID ) == 16 );
    d_params    = params_in;
    GeomDim     = GeomType::null;
    PhysicalDim = 0;
    d_max_gcw   = 0;
    d_comm      = d_params->comm;
    d_db        = d_params->d_db;
    AMP_INSIST( !d_comm.isNull(), "Communicator in mesh params must be non NULL" );
    setMeshID();
    d_name = d_db->getWithDefault<std::string>( "MeshName", "NULL" );
}
Mesh::Mesh( const Mesh &rhs )
    : d_params( rhs.d_params ),
      d_geometry( nullptr ),
      GeomDim( rhs.GeomDim ),
      PhysicalDim( rhs.PhysicalDim ),
      d_max_gcw( rhs.d_max_gcw ),
      d_comm( rhs.d_comm ),
      d_db( rhs.d_db ),
      d_meshID( rhs.d_meshID ),
      d_name( rhs.d_name ),
      d_box( rhs.d_box ),
      d_box_local( rhs.d_box_local )
{
    if ( rhs.d_geometry )
        d_geometry = rhs.d_geometry->clone();
}


/********************************************************
 * De-constructor                                        *
 ********************************************************/
Mesh::~Mesh() = default;


/********************************************************
 * Function to set the mesh ID                           *
 * This function will create a unique ID for every mesh. *
 * To accomplish this goal, the ID will consist of the   *
 * rank of the root processor (from the global comm),    *
 * and the number of meshes created by that processor.   *
 ********************************************************/
void Mesh::setMeshID()
{
    if ( d_comm.getRank() == 0 ) {
        // Root will create the meshID
        AMP_MPI globalComm( AMP_COMM_WORLD );
        d_meshID = MeshID( globalComm.getRank(), nextLocalMeshID );
        nextLocalMeshID++;
    }
    // Broadcast the meshID to all processors
    d_meshID = d_comm.bcast( d_meshID, 0 );
}


/********************************************************
 * Function to return the meshID composing the mesh      *
 ********************************************************/
std::vector<MeshID> Mesh::getAllMeshIDs() const { return std::vector<MeshID>( 1, d_meshID ); }
std::vector<MeshID> Mesh::getBaseMeshIDs() const { return std::vector<MeshID>( 1, d_meshID ); }
std::vector<MeshID> Mesh::getLocalMeshIDs() const { return std::vector<MeshID>( 1, d_meshID ); }
std::vector<MeshID> Mesh::getLocalBaseMeshIDs() const { return std::vector<MeshID>( 1, d_meshID ); }


/********************************************************
 * Function to return the mesh with the given ID         *
 ********************************************************/
std::shared_ptr<Mesh> Mesh::Subset( MeshID meshID ) const
{
    if ( d_meshID == meshID )
        return std::const_pointer_cast<Mesh>( shared_from_this() );
    else
        return std::shared_ptr<Mesh>();
}


/********************************************************
 * Function to return the mesh with the given name       *
 ********************************************************/
std::shared_ptr<Mesh> Mesh::Subset( std::string name ) const
{
    if ( d_name == name )
        return std::const_pointer_cast<Mesh>( shared_from_this() );
    else
        return std::shared_ptr<Mesh>();
}


/********************************************************
 * Function to subset a mesh using a mesh iterator       *
 ********************************************************/
std::shared_ptr<Mesh> Mesh::Subset( const MeshIterator &iterator, bool isGlobal ) const
{
    if ( isGlobal ) {
        auto N = d_comm.sumReduce( iterator.size() );
        if ( N == 0 )
            return std::shared_ptr<Mesh>();
    } else if ( iterator.size() == 0 ) {
        return std::shared_ptr<Mesh>();
    }
    auto mesh = std::make_shared<SubsetMesh>( shared_from_this(), iterator, isGlobal );
    return mesh;
}


/********************************************************
 * Function to return the element given an ID            *
 ********************************************************/
MeshElement Mesh::getElement( const MeshElementID &elem_id ) const
{
    MeshID mesh_id = elem_id.meshID();
    AMP_INSIST( mesh_id == d_meshID, "mesh id must match the mesh id of the element" );
    auto it = getIterator( elem_id.type() );
    for ( size_t i = 0; i < it.size(); i++, ++it ) {
        if ( it->globalID() == elem_id )
            return *it;
    }
    return MeshElement();
}


/********************************************************
 * Function to return parents of an element              *
 ********************************************************/
std::vector<MeshElement> Mesh::getElementParents( const MeshElement &, const GeomType ) const
{
    AMP_ERROR( "getElementParents is not implimented for the base class" );
    return std::vector<MeshElement>();
}


/********************************************************
 * Return the position vector                            *
 ********************************************************/
#ifdef USE_AMP_VECTORS
AMP::LinearAlgebra::Vector::shared_ptr Mesh::getPositionVector( std::string name,
                                                                const int gcw ) const
{
#ifdef USE_AMP_DISCRETIZATION
    auto DOFs = AMP::Discretization::simpleDOFManager::create(
        std::const_pointer_cast<Mesh>( shared_from_this() ),
        AMP::Mesh::GeomType::Vertex,
        gcw,
        PhysicalDim,
        true );
    auto nodalVariable = std::make_shared<AMP::LinearAlgebra::Variable>( name );
    auto position      = AMP::LinearAlgebra::createVector( DOFs, nodalVariable, true );
    std::vector<size_t> dofs( PhysicalDim );
    for ( const auto &elem : DOFs->getIterator() ) {
        auto id    = elem.globalID();
        auto coord = elem.coord();
        DOFs->getDOFs( id, dofs );
        position->setValuesByGlobalID( dofs.size(), &dofs[0], &coord[0] );
    }
    return position;
#else
    AMP_ERROR( "getPositionVector requires DISCRETIZATION" );
    return AMP::LinearAlgebra::Vector::shared_ptr();
#endif
}
#endif


/********************************************************
 * Check if the element is a member of the mesh          *
 ********************************************************/
bool Mesh::isMember( const MeshElementID &id ) const { return id.meshID() == d_meshID; }


/********************************************************
 * Functions that aren't implimented for the base class  *
 ********************************************************/
std::shared_ptr<Mesh> Mesh::Subset( Mesh & ) const
{
    AMP_ERROR( "Subset is not implimented for the base class" );
    return std::shared_ptr<Mesh>();
}
MeshIterator Mesh::getIterator( const GeomType, const int ) const
{
    AMP_ERROR( "getIterator is not implimented for the base class" );
    return MeshIterator();
}
MeshIterator Mesh::getSurfaceIterator( const GeomType, const int ) const
{
    AMP_ERROR( "getSurfaceIterator is not implimented for the base class" );
    return MeshIterator();
}
std::vector<int> Mesh::getBoundaryIDs() const
{
    AMP_ERROR( "getBoundaryIDs is not implimented for the base class" );
    return std::vector<int>();
}
MeshIterator Mesh::getBoundaryIDIterator( const GeomType, const int, const int ) const
{
    AMP_ERROR( "getBoundaryIDIterator is not implimented for the base class" );
    return MeshIterator();
}
std::vector<int> Mesh::getBlockIDs() const
{
    AMP_ERROR( "getBlockIDs is not implimented for the base class" );
    return std::vector<int>();
}
MeshIterator Mesh::getBlockIDIterator( const GeomType, const int, const int ) const
{
    AMP_ERROR( "getBlockIDIterator is not implimented for the base class" );
    return MeshIterator();
}
size_t Mesh::numLocalElements( const GeomType ) const
{
    AMP_ERROR( "numLocalElements is not implimented for the base class" );
    return 0;
}
size_t Mesh::numGlobalElements( const GeomType ) const
{
    AMP_ERROR( "numGlobalElements is not implimented for the base class" );
    return 0;
}
size_t Mesh::numGhostElements( const GeomType, int ) const
{
    AMP_ERROR( "numGhostElements is not implimented for the base class" );
    return 0;
}


/********************************************************
 * Compare two meshes                                    *
 ********************************************************/
static double getTol( const std::vector<double> &box, size_t N )
{
    int ndim     = box.size();
    size_t N2    = pow( N, 1.0 / ndim );
    double dx[3] = { 0, 0, 0 };
    for ( int d = 0; d < ndim / 2; d++ )
        dx[d] = ( box[2 * d + 1] - box[2 * d] ) / N2;
    return 0.2 * sqrt( dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2] );
}
static inline std::vector<Point> getPoints( MeshIterator it )
{
    std::vector<Point> p( it.size() );
    for ( size_t i = 0; i < p.size(); i++, ++it )
        p[i] = it->centroid();
    return p;
}
int Mesh::compare( const Mesh &a, const Mesh &b )
{
    // Check if the meshes are equal
    if ( a == b )
        return 1;
    // Special case for multimeshes
    auto a2 = dynamic_cast<const MultiMesh *>( &a );
    auto b2 = dynamic_cast<const MultiMesh *>( &b );
    if ( a2 || b2 ) {
        if ( !a2 || !b2 )
            return false;
        auto list1 = a2->getMeshes();
        auto list2 = b2->getMeshes();
        if ( list1.size() != list2.size() )
            return false;
        int result = 1;
        for ( size_t i = 0; i < list1.size(); i++ ) {
            int test = compare( *list1[i], *list2[i] );
            if ( test == 0 )
                return 0;
            result = std::max( result, test );
        }
        return result;
    }
    // Default comparison
    // Perform simple comparisons
    if ( a.GeomDim != b.GeomDim || a.PhysicalDim != b.PhysicalDim ||
         a.d_comm.compare( b.d_comm ) == 0 )
        return 0;
    if ( a.getBoundaryIDs() != b.getBoundaryIDs() || a.getBlockIDs() != b.getBlockIDs() )
        return 0;
    // Compare domains
    size_t N1  = a.numLocalElements( a.GeomDim );
    size_t N2  = b.numLocalElements( b.GeomDim );
    auto box1  = a.getBoundingBox();
    auto box2  = b.getBoundingBox();
    double tol = getTol( box1, std::min( N1, N2 ) );
    for ( size_t i = 0; i < box1.size(); i++ ) {
        if ( fabs( box1[i] - box2[i] ) > tol )
            return 0;
    }
    // Compare the coordinates
    if ( N1 == N2 ) {
        bool test    = true;
        auto nodes_a = getPoints( a.getIterator( GeomType::Vertex ) );
        auto nodes_b = getPoints( b.getIterator( GeomType::Vertex ) );
        auto elems_a = getPoints( a.getIterator( a.GeomDim ) );
        auto elems_b = getPoints( b.getIterator( b.GeomDim ) );
        kdtree tree_a_node( nodes_a );
        kdtree tree_a_elem( elems_a );
        for ( const auto &p : nodes_b ) {
            auto p2 = tree_a_node.find_nearest( p );
            test    = test && ( p - p2 ).norm() < tol * tol;
        }
        for ( const auto &p : elems_b ) {
            auto p2 = tree_a_elem.find_nearest( p );
            test    = test && ( p - p2 ).norm() < tol * tol;
        }
        if ( test )
            return 2;
    }
    // Get the geometries
    auto geom1 = a.getGeometry();
    auto geom2 = b.getGeometry();
    if ( !geom1 ) {
        auto ptr = std::const_pointer_cast<Mesh>( a.shared_from_this() );
        geom1    = std::make_shared<AMP::Geometry::MeshGeometry>( ptr );
    }
    if ( !geom2 ) {
        auto ptr = std::const_pointer_cast<Mesh>( b.shared_from_this() );
        geom2    = std::make_shared<AMP::Geometry::MeshGeometry>( ptr );
    }
    if ( *geom1 == *geom2 )
        return 3;

    AMP_WARNING( "Not finished" );
    return -1;
}


/********************************************************
 * MeshIterator set operations                           *
 ********************************************************/
MeshIterator Mesh::getIterator( SetOP OP, const MeshIterator &A, const MeshIterator &B )
{
    if ( OP == SetOP::Union ) {
        // Perform a union: A U B
        // Get the union using the mesh IDs
        std::set<MeshElementID> union_set;
        for ( auto &elem : A )
            union_set.insert( elem.globalID() );
        for ( auto &elem : B )
            union_set.insert( elem.globalID() );
        std::vector<MeshElementID> union_ids( union_set.begin(), union_set.end() );
        // Create the iterator
        if ( union_ids.size() == A.size() ) {
            return MeshIterator( A.begin() );
        } else if ( union_ids.size() == B.size() ) {
            return MeshIterator( B.begin() );
        } else {
            auto elements = std::make_shared<std::vector<MeshElement>>( union_ids.size() );
            for ( auto &elem : A ) {
                MeshElementID idA = elem.globalID();
                size_t index      = Utilities::findfirst( union_ids, idA );
                if ( index == union_ids.size() ) {
                    index--;
                }
                if ( union_ids[index] == idA )
                    ( *elements )[index] = elem;
            }
            for ( auto &elem : B ) {
                MeshElementID idB = elem.globalID();
                size_t index      = Utilities::findfirst( union_ids, idB );
                if ( index == union_ids.size() ) {
                    index--;
                }
                if ( union_ids[index] == idB )
                    ( *elements )[index] = elem;
            }
            return MultiVectorIterator( elements, 0 );
        }
    } else if ( OP == SetOP::Intersection ) {
        // Perform a intersection: A n B
        // Get the intersection using the mesh IDs
        if ( A.size() == 0 || B.size() == 0 )
            return MeshIterator();
        std::vector<MeshElementID> idA( A.size() );
        auto it = A.begin();
        for ( size_t i = 0; i < A.size(); ++i, ++it )
            idA[i] = it->globalID();
        Utilities::quicksort( idA );
        std::vector<MeshElementID> intersection;
        intersection.reserve( B.size() );
        for ( auto &elem : B ) {
            MeshElementID idB = elem.globalID();
            size_t index      = Utilities::findfirst( idA, idB );
            if ( index == idA.size() ) {
                index--;
            }
            if ( idA[index] == idB )
                intersection.push_back( idB );
        }
        if ( intersection.empty() )
            return MeshIterator();
        // Sort the intersection and check for duplicates
        Utilities::quicksort( intersection );
        for ( size_t i = 1; i < intersection.size(); i++ )
            AMP_ASSERT( intersection[i] != intersection[i - 1] );
        // Create the iterator
        if ( intersection.size() == A.size() ) {
            return MeshIterator( A.begin() );
        } else if ( intersection.size() == B.size() ) {
            return MeshIterator( B.begin() );
        } else {
            auto elements = std::make_shared<std::vector<MeshElement>>( intersection.size() );
            for ( auto &elem : B ) {
                MeshElementID idB = elem.globalID();
                size_t index      = Utilities::findfirst( intersection, idB );
                if ( index == intersection.size() ) {
                    index--;
                }
                if ( intersection[index] == idB )
                    ( *elements )[index] = elem;
            }
            return MultiVectorIterator( elements, 0 );
        }
    } else if ( OP == SetOP::Complement ) {
        // Perform a SetOP::Complement:  A - B
        // Get the compliment using the mesh IDs
        std::set<MeshElementID> compliment_set;
        for ( auto &elem : A )
            compliment_set.insert( elem.globalID() );
        for ( auto &elem : B )
            compliment_set.erase( elem.globalID() );
        std::vector<MeshElementID> compliment( compliment_set.begin(), compliment_set.end() );
        if ( compliment.empty() )
            return MeshIterator();
        // Create the iterator
        if ( compliment.size() == A.size() ) {
            return MeshIterator( A.begin() );
        } else {
            auto elements = std::make_shared<std::vector<MeshElement>>( compliment.size() );
            for ( auto &elem : A ) {
                MeshElementID idA = elem.globalID();
                size_t index      = Utilities::findfirst( compliment, idA );
                if ( index == compliment.size() ) {
                    index--;
                }
                if ( compliment[index] == idA )
                    ( *elements )[index] = elem;
            }
            return MultiVectorIterator( elements, 0 );
        }
    } else {
        AMP_ERROR( "Unknown set operation" );
    }
    return MeshIterator();
}


} // namespace AMP::Mesh


/********************************************************
 * Instantiate communication of MeshElementID            *
 ********************************************************/
namespace AMP {
typedef AMP::Mesh::MeshElementID ID;
template ID AMP_MPI::bcast<ID>( ID, int ) const;
template void AMP_MPI::bcast<ID>( ID *, int, int ) const;
template void AMP_MPI::send<ID>( const ID *, int, int, int ) const;
template MPI_Request AMP_MPI::Isend<ID>( const ID *, int, int, int ) const;
template void AMP_MPI::recv<ID>( ID *, int &, int, const bool, int ) const;
template MPI_Request AMP_MPI::Irecv<ID>( ID *buf, int, int, int ) const;
template std::vector<ID> AMP_MPI::allGather<ID>( const ID & ) const;
template std::vector<ID> AMP_MPI::allGather<ID>( const std::vector<ID> & ) const;
template void AMP_MPI::allGather<ID>( const ID &, ID * ) const;
template int AMP_MPI::allGather<ID>( const ID *, int, ID *, int *, int *, bool ) const;
template void AMP_MPI::setGather<ID>( std::set<ID> & ) const;
template void AMP_MPI::allToAll<ID>( int, const ID *, ID * ) const;
template int
AMP_MPI::allToAll<ID>( const ID *, const int[], const int[], ID *, int *, int *, bool ) const;
} // namespace AMP
