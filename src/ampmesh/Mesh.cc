#include "AMP/ampmesh/Mesh.h"
#include "AMP/ampmesh/MeshElementVectorIterator.h"
#include "AMP/ampmesh/SubsetMesh.h"
#include "AMP/utils/Utilities.h"

#ifdef USE_AMP_VECTORS
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"
#endif
#ifdef USE_AMP_DISCRETIZATION
#include "AMP/discretization/simpleDOF_Manager.h"
#endif

#include <cmath>


namespace AMP {
namespace Mesh {


static_assert( sizeof( MeshID ) == 8, "unexpected size for MeshID" );
static_assert( sizeof( ElementID ) == 8, "unexpected size for ElementID" );
static_assert( sizeof( MeshElementID ) == 16, "unexpected size for MeshElementID" );


static unsigned int nextLocalMeshID = 1;


/********************************************************
 * Constructors                                          *
 ********************************************************/
Mesh::Mesh( const MeshParameters::shared_ptr &params_in )
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
    d_name = d_db->getStringWithDefault( "MeshName", "NULL" );
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
AMP::shared_ptr<Mesh> Mesh::Subset( MeshID meshID ) const
{
    if ( d_meshID == meshID )
        return AMP::const_pointer_cast<Mesh>( shared_from_this() );
    else
        return AMP::shared_ptr<Mesh>();
}


/********************************************************
 * Function to return the mesh with the given name       *
 ********************************************************/
AMP::shared_ptr<Mesh> Mesh::Subset( std::string name ) const
{
    if ( d_name == name )
        return AMP::const_pointer_cast<Mesh>( shared_from_this() );
    else
        return AMP::shared_ptr<Mesh>();
}


/********************************************************
 * Function to subset a mesh using a mesh iterator       *
 ********************************************************/
AMP::shared_ptr<Mesh> Mesh::Subset( const MeshIterator &iterator, bool isGlobal ) const
{
    if ( isGlobal ) {
        auto N = d_comm.sumReduce( iterator.size() );
        if ( N == 0 )
            return AMP::shared_ptr<Mesh>();
    } else if ( !isGlobal && iterator.size() == 0 ) {
        return AMP::shared_ptr<Mesh>();
    }
    auto mesh = AMP::make_shared<SubsetMesh>( shared_from_this(), iterator, isGlobal );
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
        AMP::const_pointer_cast<Mesh>( shared_from_this() ),
        AMP::Mesh::GeomType::Vertex,
        gcw,
        PhysicalDim,
        true );
    auto nodalVariable = AMP::make_shared<AMP::LinearAlgebra::Variable>( name );
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
AMP::shared_ptr<Mesh> Mesh::Subset( Mesh & ) const
{
    AMP_ERROR( "Subset is not implimented for the base class" );
    return AMP::shared_ptr<Mesh>();
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
            auto elements = AMP::make_shared<std::vector<MeshElement>>( union_ids.size() );
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
            auto elements = AMP::make_shared<std::vector<MeshElement>>( intersection.size() );
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
            auto elements = AMP::make_shared<std::vector<MeshElement>>( compliment.size() );
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


} // namespace Mesh
} // namespace AMP
