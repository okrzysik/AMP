#include <set>
#include <vector>

#include "AMP/ampmesh/MeshElementVectorIterator.h"
#include "AMP/ampmesh/MultiIterator.h"
#include "AMP/ampmesh/SubsetMesh.h"
#include "AMP/utils/AMP_MPI.I"
#ifdef USE_AMP_VECTORS
#include "AMP/vectors/Vector.h"
#endif

namespace AMP {
namespace Mesh {


/********************************************************
 * Constructors                                          *
 ********************************************************/
SubsetMesh::SubsetMesh( std::shared_ptr<const Mesh> mesh,
                        const AMP::Mesh::MeshIterator &iterator_in,
                        bool isGlobal )
    : d_parent_mesh( mesh )
{
    if ( isGlobal ) {
        this->d_comm = mesh->getComm();
    } else {
        for ( const auto &elem : iterator_in ) {
            if ( !elem.globalID().is_local() )
                AMP_ERROR(
                    "Subsetting local mesh for iterator with ghost elements is not supported" );
        }
        this->d_comm = AMP_MPI( AMP_COMM_SELF );
    }
    this->setMeshID();
    this->PhysicalDim = mesh->getDim();
    this->d_name      = mesh->getName() + "_subset";
    // Check the iterator
    auto type = GeomType::null;
    for ( const auto &elem : iterator_in ) {
        if ( type == GeomType::null )
            type = elem.elementType();
        if ( type != elem.elementType() )
            AMP_ERROR( "Subset mesh requires all of the elements to be the same type" );
    }
    int type2 = d_comm.minReduce( (int) type );
    if ( type != GeomType::null && type2 != (int) type )
        AMP_ERROR( "Subset mesh requires all of the elements to be the same type" );
    this->GeomDim = static_cast<GeomType>( type2 );
    auto mesh_ids = mesh->getBaseMeshIDs();
    for ( const auto &elem : iterator_in ) {
        MeshID id  = elem.globalID().meshID();
        bool found = false;
        for ( auto &mesh_id : mesh_ids ) {
            if ( id == mesh_id )
                found = true;
        }
        if ( !found )
            AMP_ERROR( "Iterator contains elements not found in parent meshes" );
    }
    // Create a list of all elements of the desired type
    d_max_gcw = d_parent_mesh->getMaxGhostWidth();
    d_elements =
        std::vector<std::vector<std::shared_ptr<std::vector<MeshElement>>>>( (int) GeomDim + 1 );
    for ( int i = 0; i <= static_cast<int>( GeomDim ); i++ ) {
        d_elements[i] = std::vector<std::shared_ptr<std::vector<MeshElement>>>(
            d_max_gcw + 1, std::make_shared<std::vector<MeshElement>>() );
    }
    int gcw = 0;
    while ( true ) {
        auto iterator1 = Mesh::getIterator(
            SetOP::Intersection, iterator_in, mesh->getIterator( GeomDim, gcw ) );
        auto iterator2 = iterator1.begin();
        if ( gcw > 0 )
            iterator2 =
                Mesh::getIterator( SetOP::Complement, iterator1, mesh->getIterator( GeomDim, 0 ) );
        d_elements[(int) GeomDim][gcw] =
            std::make_shared<std::vector<MeshElement>>( iterator2.size() );
        for ( size_t i = 0; i < iterator2.size(); i++ ) {
            d_elements[(int) GeomDim][gcw]->operator[]( i ) = *iterator2;
            ++iterator2;
        }
        AMP::Utilities::quicksort( *( d_elements[(int) GeomDim][gcw] ) );
        if ( iterator1.size() == iterator_in.size() )
            break;
        gcw++;
    }
    // Create a list of all elements that compose the elements of GeomType
    for ( int t = 0; t < (int) GeomDim; t++ ) {
        d_elements[t] = std::vector<std::shared_ptr<std::vector<MeshElement>>>( d_max_gcw + 1 );
        for ( gcw = 0; gcw <= d_max_gcw; gcw++ ) {
            std::set<MeshElement> list;
            for ( const auto &elem : this->getIterator( GeomDim, gcw ) ) {
                auto elements = elem.getElements( (GeomType) t );
                for ( auto &element : elements ) {
                    if ( gcw == 0 ) {
                        if ( element.globalID().is_local() )
                            list.insert( element );
                    } else {
                        bool found = false;
                        for ( int j = 0; j < gcw; j++ ) {
                            size_t index =
                                AMP::Utilities::findfirst( *( d_elements[t][j] ), element );
                            if ( index == d_elements[t][j]->size() ) {
                                index--;
                            }
                            if ( d_elements[t][j]->operator[]( index ) == element ) {
                                found = true;
                                break;
                            }
                        }
                        if ( !found )
                            list.insert( element );
                    }
                }
            }
            d_elements[t][gcw] =
                std::make_shared<std::vector<MeshElement>>( list.begin(), list.end() );
        }
    }
    // For each entity type, we need to check that any ghost elements are owned by somebody
    for ( int t = 0; t < (int) GeomDim; t++ ) {
        // First get a global list of all ghost elements
        std::vector<MeshElementID> ghost_local;
        for ( gcw = 0; gcw <= d_max_gcw; gcw++ ) {
            ghost_local.reserve( ghost_local.size() + d_elements[t][gcw]->size() );
            for ( size_t i = 0; i < d_elements[t][gcw]->size(); i++ ) {
                MeshElementID id = ( *d_elements[t][gcw] )[i].globalID();
                if ( !id.is_local() )
                    ghost_local.push_back( id );
            }
        }
        size_t N_ghost_local  = ghost_local.size();
        size_t N_ghost_global = d_comm.sumReduce( N_ghost_local );
        if ( N_ghost_global == 0 )
            continue;
        std::vector<MeshElementID> ghost_global( N_ghost_global );
        MeshElementID *send_ptr = nullptr;
        if ( N_ghost_local > 0 ) {
            send_ptr = &ghost_local[0];
        }
        MeshElementID *recv_ptr = &ghost_global[0];
        d_comm.allGather( send_ptr, (int) N_ghost_local, recv_ptr );
        AMP::Utilities::unique( ghost_global );
        // For each ghost, check if we own it, and add it to the list if necessary
        MeshID my_mesh_id    = d_parent_mesh->meshID();
        unsigned int my_rank = d_parent_mesh->getComm().getRank();
        bool changed         = false;
        for ( auto &elem : ghost_global ) {
            AMP_ASSERT( my_mesh_id == elem.meshID() );
            if ( elem.owner_rank() == my_rank ) {
                bool found = false;
                for ( size_t j = 0; j < d_elements[t][0]->size(); j++ ) {
                    if ( ( *d_elements[t][0] )[j].globalID() == elem ) {
                        found = true;
                        break;
                    }
                }
                if ( !found ) {
                    MeshElementID tmp = elem;
                    tmp.set_is_local( true ); // We do own this element
                    MeshElement element = d_parent_mesh->getElement( tmp );
                    AMP_ASSERT( element.globalID() != MeshElementID() );
                    AMP_ASSERT( element.globalID() == elem );
                    d_elements[t][0]->push_back( element );
                    changed = true;
                }
            }
        }
        if ( changed )
            AMP::Utilities::quicksort( *d_elements[t][0] );
    }
    // Count the number of elements of each type
    N_global = std::vector<size_t>( (int) GeomDim + 1 );
    for ( int i = 0; i <= (int) GeomDim; i++ )
        N_global[i] = d_elements[i][0]->size();
    d_comm.sumReduce( &N_global[0], (int) N_global.size() );
    for ( int i = 0; i <= (int) GeomDim; i++ )
        AMP_ASSERT( N_global[i] > 0 );
    // Create the bounding box
    d_box = std::vector<double>( 2 * PhysicalDim );
    for ( int j = 0; j < PhysicalDim; j++ ) {
        d_box[2 * j + 0] = 1e100;
        d_box[2 * j + 1] = -1e100;
    }
    for ( const auto &elem : getIterator( GeomType::Vertex, 0 ) ) {
        auto coord = elem.coord();
        for ( int j = 0; j < PhysicalDim; j++ ) {
            if ( coord[j] < d_box[2 * j + 0] )
                d_box[2 * j + 0] = coord[j];
            if ( coord[j] > d_box[2 * j + 1] )
                d_box[2 * j + 1] = coord[j];
        }
    }
    for ( int j = 0; j < PhysicalDim; j++ ) {
        d_box[2 * j + 0] = d_comm.minReduce( d_box[2 * j + 0] );
        d_box[2 * j + 1] = d_comm.maxReduce( d_box[2 * j + 1] );
    }
    // Create the boundary id sets
    auto boundary_ids = d_parent_mesh->getBoundaryIDs();
    std::set<int> new_boundary_ids;
    for ( int t = 0; t <= (int) GeomDim; t++ ) {
        for ( auto &boundary_id : boundary_ids ) {
            for ( gcw = 0; gcw <= d_max_gcw; gcw++ ) {
                if ( gcw > 0 )
                    continue; // Iterators over id sets with ghost values is not supported in
                              // libmesh yet
                auto iterator1 = MultiVectorIterator( d_elements[t][gcw], 0 );
                auto iterator2 =
                    d_parent_mesh->getBoundaryIDIterator( (GeomType) t, boundary_id, gcw );
                auto iterator = Mesh::getIterator( SetOP::Intersection, iterator1, iterator2 );
                std::shared_ptr<std::vector<MeshElement>> elements;
                if ( iterator.size() == 0 ) {
                    elements = std::make_shared<std::vector<MeshElement>>( 0 );
                } else if ( iterator.size() == iterator1.size() ) {
                    elements = d_elements[t][gcw];
                } else {
                    elements = std::make_shared<std::vector<MeshElement>>( iterator.size() );
                    for ( size_t j = 0; j < iterator.size(); j++ ) {
                        elements->operator[]( j ) = *iterator;
                        ++iterator;
                    }
                }
                if ( gcw == 0 ) {
                    size_t global_size = d_comm.sumReduce( elements->size() );
                    if ( global_size == 0 )
                        break;
                }
                map_id_struct map_id;
                map_id.id   = boundary_id;
                map_id.type = (GeomType) t;
                map_id.gcw  = gcw;
                d_boundarySets.insert( std::make_pair( map_id, elements ) );
                new_boundary_ids.insert( boundary_id );
            }
        }
    }
    d_boundaryIdSets = std::vector<int>( new_boundary_ids.begin(), new_boundary_ids.end() );
    int *send_ptr    = nullptr;
    if ( d_boundaryIdSets.size() > 0 )
        send_ptr = &d_boundaryIdSets[0];
    size_t recv_size = d_comm.sumReduce( d_boundaryIdSets.size() );
    if ( recv_size > 0 ) {
        std::vector<int> recv_list( recv_size, 0 );
        d_comm.allGather( &send_ptr[0], (int) d_boundaryIdSets.size(), &recv_list[0] );
        for ( auto &elem : recv_list )
            new_boundary_ids.insert( elem );
        d_boundaryIdSets = std::vector<int>( new_boundary_ids.begin(), new_boundary_ids.end() );
    }
    // Create the surface sets
    d_surface =
        std::vector<std::vector<std::shared_ptr<std::vector<MeshElement>>>>( (int) GeomDim + 1 );
    for ( int t = 0; t <= (int) GeomDim; t++ ) {
        d_surface[t] = std::vector<std::shared_ptr<std::vector<MeshElement>>>( d_max_gcw + 1 );
        for ( gcw = 0; gcw <= d_max_gcw; gcw++ ) {
            auto iterator1 = MultiVectorIterator( d_elements[t][gcw], 0 );
            auto iterator2 = d_parent_mesh->getSurfaceIterator( (GeomType) t, gcw );
            auto iterator  = Mesh::getIterator( SetOP::Intersection, iterator1, iterator2 );
            std::shared_ptr<std::vector<MeshElement>> elements;
            if ( iterator.size() == 0 ) {
                elements = std::make_shared<std::vector<MeshElement>>( 0 );
            } else if ( iterator.size() == iterator1.size() ) {
                elements = d_elements[t][gcw];
            } else {
                elements = std::make_shared<std::vector<MeshElement>>( iterator.size() );
                for ( size_t j = 0; j < iterator.size(); j++ ) {
                    elements->operator[]( j ) = *iterator;
                    ++iterator;
                }
            }
            d_surface[t][gcw] = elements;
        }
    }
    // Create the block id sets
    if ( GeomDim == d_parent_mesh->getGeomType() ) {
        // Currently only elements support the block IDs
        std::vector<int> block_ids = d_parent_mesh->getBlockIDs();
        std::set<int> new_block_ids;
        for ( const auto &elem : getIterator( GeomDim, 0 ) ) {
            for ( auto &block_id : block_ids ) {
                if ( elem.isInBlock( block_id ) )
                    new_block_ids.insert( block_id );
            }
        }
        d_blockIdSets = std::vector<int>( new_block_ids.begin(), new_block_ids.end() );
        send_ptr      = nullptr;
        if ( d_boundaryIdSets.size() > 0 )
            send_ptr = &d_boundaryIdSets[0];
        recv_size = d_comm.sumReduce( d_blockIdSets.size() );
        if ( recv_size > 0 ) {
            std::vector<int> recv_list( recv_size, 0 );
            d_comm.allGather( &send_ptr[0], (int) d_blockIdSets.size(), &recv_list[0] );
            for ( auto &elem : recv_list )
                new_block_ids.insert( elem );
            d_blockIdSets = std::vector<int>( new_block_ids.begin(), new_block_ids.end() );
        }
    }
}


/********************************************************
 * De-constructor                                        *
 ********************************************************/
SubsetMesh::~SubsetMesh() = default;


/********************************************************
 * Copy the mesh                                         *
 ********************************************************/
std::unique_ptr<Mesh> SubsetMesh::clone() const
{
    AMP_ERROR( "clone is not currently supported with SubsetMesh" );
    return std::unique_ptr<Mesh>();
}


/********************************************************
 * Function to return the meshID composing the mesh      *
 ********************************************************/
std::vector<MeshID> SubsetMesh::getAllMeshIDs() const
{
    auto ids = d_parent_mesh->getAllMeshIDs();
    ids.push_back( d_meshID );
    AMP::Utilities::quicksort( ids );
    return ids;
}
std::vector<MeshID> SubsetMesh::getBaseMeshIDs() const { return d_parent_mesh->getBaseMeshIDs(); }
std::vector<MeshID> SubsetMesh::getLocalMeshIDs() const
{
    auto ids = d_parent_mesh->getLocalMeshIDs();
    ids.push_back( d_meshID );
    AMP::Utilities::quicksort( ids );
    return ids;
}
std::vector<MeshID> SubsetMesh::getLocalBaseMeshIDs() const
{
    return d_parent_mesh->getLocalBaseMeshIDs();
}


/********************************************************
 * Function to return the mesh with the given ID         *
 ********************************************************/
std::shared_ptr<Mesh> SubsetMesh::Subset( MeshID meshID ) const
{
    if ( d_meshID == meshID || d_parent_mesh->meshID() == meshID )
        return std::const_pointer_cast<Mesh>( shared_from_this() );
    else
        return std::shared_ptr<Mesh>();
}


/********************************************************
 * Function to return the mesh with the given name       *
 ********************************************************/
std::shared_ptr<Mesh> SubsetMesh::Subset( std::string name ) const
{
    if ( d_name == name || d_parent_mesh->getName() == name )
        return std::const_pointer_cast<Mesh>( shared_from_this() );
    else
        return std::shared_ptr<Mesh>();
}


/********************************************************
 * Mesh iterators                                        *
 ********************************************************/
MeshIterator SubsetMesh::getIterator( const GeomType type, const int gcw ) const
{
    int gcw2   = gcw;
    auto type2 = static_cast<int>( type );
    if ( (int) d_elements.size() <= type2 || gcw < 0 )
        return MeshIterator();
    if ( d_elements[type2].empty() )
        return MeshIterator();
    if ( gcw2 >= (int) d_elements[type2].size() )
        gcw2 = (int) d_elements[type2].size() - 1;
    if ( gcw2 == 0 )
        return MultiVectorIterator( d_elements[type2][0], 0 );
    std::vector<MeshIterator> iterators( gcw2 + 1 );
    for ( int i = 0; i <= gcw2; i++ )
        iterators[i] = MultiVectorIterator( d_elements[type2][i], 0 );
    return MultiIterator( iterators, 0 );
}
MeshIterator SubsetMesh::getSurfaceIterator( const GeomType type, const int gcw ) const
{
    auto type2 = static_cast<int>( type );
    if ( gcw == 0 )
        return MultiVectorIterator( d_surface[type2][0], 0 );
    if ( gcw >= (int) d_surface[type2].size() )
        AMP_ERROR( "Maximum ghost width exceeded" );
    std::vector<MeshIterator> iterators( gcw + 1 );
    for ( int i = 0; i <= gcw; i++ )
        iterators[i] = MultiVectorIterator( d_surface[type2][i], 0 );
    return MultiIterator( iterators, 0 );
}
std::vector<int> SubsetMesh::getBoundaryIDs() const { return d_boundaryIdSets; }
MeshIterator
SubsetMesh::getBoundaryIDIterator( const GeomType type, const int id, const int gcw ) const
{
    std::vector<MeshIterator> iterators;
    iterators.reserve( gcw + 1 );
    for ( int i = 0; i <= gcw; i++ ) {
        map_id_struct map_id;
        map_id.id   = id;
        map_id.type = type;
        map_id.gcw  = i;
        auto map_it = d_boundarySets.find( map_id );
        if ( map_it == d_boundarySets.end() )
            continue;
        iterators.push_back( MultiVectorIterator( map_it->second, 0 ) );
    }
    if ( iterators.empty() )
        return MeshIterator();
    if ( iterators.size() == 1 )
        return iterators[0];
    return MultiIterator( iterators, 0 );
}
std::vector<int> SubsetMesh::getBlockIDs() const { return d_blockIdSets; }
MeshIterator
SubsetMesh::getBlockIDIterator( const GeomType type, const int id, const int gcw ) const
{
    std::vector<MeshIterator> iterators;
    iterators.reserve( gcw + 1 );
    for ( int i = 0; i <= gcw; i++ ) {
        map_id_struct map_id;
        map_id.id   = id;
        map_id.type = type;
        map_id.gcw  = i;
        auto map_it = d_blockSets.find( map_id );
        if ( map_it == d_blockSets.end() )
            continue;
        iterators.push_back( MultiVectorIterator( map_it->second, 0 ) );
    }
    if ( iterators.empty() )
        return MeshIterator();
    if ( iterators.size() == 1 )
        return iterators[0];
    return MultiIterator( iterators, 0 );
}


/********************************************************
 * Check if the element is a member of the mesh          *
 ********************************************************/
bool SubsetMesh::isMember( const MeshElementID &id ) const
{
    if ( !d_parent_mesh->isMember( id ) )
        return false;
    auto type = static_cast<int>( id.type() );
    if ( type >= static_cast<int>( d_elements.size() ) )
        return false;
    for ( auto &elem : d_elements[type] ) {
        if ( elem == nullptr )
            continue;
        for ( auto &element : *elem.get() ) {
            if ( element == id )
                return true;
        }
    }
    return false;
}


/********************************************************
 * Function to return the element given an ID            *
 ********************************************************/
MeshElement SubsetMesh::getElement( const MeshElementID &elem_id ) const
{
    return d_parent_mesh->getElement( elem_id );
}


/********************************************************
 * Function to return parents of an element              *
 ********************************************************/
std::vector<MeshElement> SubsetMesh::getElementParents( const MeshElement &elem,
                                                        const GeomType type ) const
{
    return d_parent_mesh->getElementParents( elem, type );
}


/********************************************************
 * Other functions                                       *
 ********************************************************/
size_t SubsetMesh::numLocalElements( const GeomType type ) const
{
    return d_elements[static_cast<int>( type )][0]->size();
}
size_t SubsetMesh::numGlobalElements( const GeomType type ) const
{
    return N_global[static_cast<int>( type )];
}
size_t SubsetMesh::numGhostElements( const GeomType type, int gcw ) const
{
    AMP_ASSERT( type <= GeomDim );
    if ( gcw == 0 )
        return 0;
    if ( gcw >= (int) d_elements[(int) type].size() )
        AMP_ERROR( "Maximum ghost width exceeded" );
    return d_elements[(int) type][gcw]->size();
}
uint64_t SubsetMesh::positionHash() const { return d_parent_mesh->positionHash(); }
Mesh::Movable SubsetMesh::isMeshMovable() const { return Mesh::Movable::Fixed; }
void SubsetMesh::displaceMesh( const std::vector<double> & )
{
    AMP_ERROR( "displaceMesh by a constant value does not work for subset mesh" );
}
#ifdef USE_AMP_VECTORS
void SubsetMesh::displaceMesh( const AMP::LinearAlgebra::Vector::const_shared_ptr )
{
    AMP_ERROR( "displaceMesh is not implimented for subset mesh" );
}
#endif


/****************************************************************
 * Check if two meshes are equal                                 *
 ****************************************************************/
bool SubsetMesh::operator==( const Mesh &rhs ) const
{
    // Check if &rhs == this
    if ( this == &rhs )
        return true;
    // Check if we can cast to a MultiMesh
    auto mesh = dynamic_cast<const SubsetMesh *>( &rhs );
    if ( !mesh )
        return false;
    // Perform comparison on sub-meshes
    if ( d_parent_mesh != mesh->d_parent_mesh )
        return false;
    AMP_ERROR( "Not finished" );
    return false;
}


} // namespace Mesh
} // namespace AMP
