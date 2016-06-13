#include <set>
#include <vector>

#include "ampmesh/MeshElementVectorIterator.h"
#include "ampmesh/MultiIterator.h"
#include "ampmesh/SubsetMesh.h"
#ifdef USE_AMP_VECTORS
#include "vectors/Vector.h"
#endif

namespace AMP {
namespace Mesh {


/********************************************************
* Constructors                                          *
********************************************************/
SubsetMesh::SubsetMesh( AMP::shared_ptr<const Mesh> mesh,
                        const AMP::Mesh::MeshIterator &iterator_in,
                        bool isGlobal )
{
    this->d_parent_mesh = mesh;
    if ( isGlobal ) {
        this->d_comm = mesh->getComm();
    } else {
        for ( AMP::Mesh::MeshIterator it = iterator_in; it != iterator_in.end(); ++it ) {
            if ( !it->globalID().is_local() )
                AMP_ERROR(
                    "Subsetting local mesh for iterator with ghost elements is not supported" );
        }
        this->d_comm = AMP_MPI( AMP_COMM_SELF );
    }
    this->setMeshID();
    this->PhysicalDim = mesh->getDim();
    this->d_name      = mesh->getName() + "_subset";
    // Check the iterator
    GeomType type         = null;
    MeshIterator iterator = iterator_in.begin();
    for ( size_t i = 0; i < iterator.size(); i++ ) {
        if ( type == null )
            type = iterator->elementType();
        if ( type != iterator->elementType() )
            AMP_ERROR( "Subset mesh requires all of the elements to be the same type" );
        ++iterator;
    }
    int type2 = d_comm.minReduce( (int) type );
    if ( type != null && type2 != (int) type )
        AMP_ERROR( "Subset mesh requires all of the elements to be the same type" );
    this->GeomDim                = (GeomType) type2;
    std::vector<MeshID> mesh_ids = mesh->getBaseMeshIDs();
    iterator                     = iterator_in.begin();
    for ( size_t i = 0; i < iterator.size(); i++ ) {
        MeshID id  = iterator->globalID().meshID();
        bool found = false;
        for ( auto &mesh_id : mesh_ids ) {
            if ( id == mesh_id )
                found = true;
        }
        if ( !found )
            AMP_ERROR( "Iterator contains elements not found in parent meshes" );
        ++iterator;
    }
    // Create a list of all elements of the desired type
    d_max_gcw = d_parent_mesh->getMaxGhostWidth();
    d_elements =
        std::vector<std::vector<AMP::shared_ptr<std::vector<MeshElement>>>>( (int) GeomDim + 1 );
    for ( int i = 0; i <= GeomDim; i++ ) {
        d_elements[i] = std::vector<AMP::shared_ptr<std::vector<MeshElement>>>(
            d_max_gcw + 1,
            AMP::shared_ptr<std::vector<MeshElement>>( new std::vector<MeshElement>() ) );
    }
    int gcw = 0;
    while ( 1 ) {
        MeshIterator iterator1 =
            Mesh::getIterator( Intersection, iterator_in, mesh->getIterator( GeomDim, gcw ) );
        MeshIterator iterator2 = iterator1.begin();
        if ( gcw > 0 )
            iterator2 = Mesh::getIterator( Complement, iterator1, mesh->getIterator( GeomDim, 0 ) );
        d_elements[GeomDim][gcw] = AMP::shared_ptr<std::vector<MeshElement>>(
            new std::vector<MeshElement>( iterator2.size() ) );
        for ( size_t i = 0; i < iterator2.size(); i++ ) {
            d_elements[GeomDim][gcw]->operator[]( i ) = *iterator2;
            ++iterator2;
        }
        AMP::Utilities::quicksort( *( d_elements[GeomDim][gcw] ) );
        if ( iterator1.size() == iterator_in.size() )
            break;
        gcw++;
    }
    // Create a list of all elements that compose the elements of GeomType
    for ( int t = 0; t < (int) GeomDim; t++ ) {
        d_elements[t] = std::vector<AMP::shared_ptr<std::vector<MeshElement>>>( d_max_gcw + 1 );
        for ( int gcw = 0; gcw <= d_max_gcw; gcw++ ) {
            std::set<MeshElement> list;
            iterator = this->getIterator( GeomDim, gcw );
            for ( size_t it = 0; it < iterator.size(); it++ ) {
                std::vector<MeshElement> elements = iterator->getElements( (GeomType) t );
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
                ++iterator;
            }
            d_elements[t][gcw] = AMP::shared_ptr<std::vector<MeshElement>>(
                new std::vector<MeshElement>( list.begin(), list.end() ) );
        }
    }
    // For each entity type, we need to check that any ghost elements are owned by somebody
    for ( int t = 0; t < (int) GeomDim; t++ ) {
        // First get a global list of all ghost elements
        std::vector<MeshElementID> ghost_local;
        for ( int gcw = 0; gcw <= d_max_gcw; gcw++ ) {
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
    for ( int i     = 0; i <= (int) GeomDim; i++ )
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
    iterator = getIterator( Vertex, 0 );
    for ( size_t i = 0; i < iterator.size(); i++ ) {
        std::vector<double> coord = iterator->coord();
        for ( int j = 0; j < PhysicalDim; j++ ) {
            if ( coord[j] < d_box[2 * j + 0] )
                d_box[2 * j + 0] = coord[j];
            if ( coord[j] > d_box[2 * j + 1] )
                d_box[2 * j + 1] = coord[j];
        }
        ++iterator;
    }
    for ( int j = 0; j < PhysicalDim; j++ ) {
        d_box[2 * j + 0] = d_comm.minReduce( d_box[2 * j + 0] );
        d_box[2 * j + 1] = d_comm.maxReduce( d_box[2 * j + 1] );
    }
    // Create the boundary id sets
    std::vector<int> boundary_ids = d_parent_mesh->getBoundaryIDs();
    std::set<int> new_boundary_ids;
    for ( int t = 0; t <= (int) GeomDim; t++ ) {
        for ( auto &boundary_id : boundary_ids ) {
            for ( int gcw = 0; gcw <= d_max_gcw; gcw++ ) {
                if ( gcw > 0 )
                    continue; // Iterators over id sets with ghost values is not supported in
                              // libmesh yet
                MeshIterator iterator1 = MultiVectorIterator( d_elements[t][gcw], 0 );
                MeshIterator iterator2 =
                    d_parent_mesh->getBoundaryIDIterator( (GeomType) t, boundary_id, gcw );
                iterator = Mesh::getIterator( Intersection, iterator1, iterator2 );
                AMP::shared_ptr<std::vector<MeshElement>> elements;
                if ( iterator.size() == 0 ) {
                    elements = AMP::shared_ptr<std::vector<MeshElement>>(
                        new std::vector<MeshElement>( 0 ) );
                } else if ( iterator.size() == iterator1.size() ) {
                    elements = d_elements[t][gcw];
                } else {
                    elements = AMP::shared_ptr<std::vector<MeshElement>>(
                        new std::vector<MeshElement>( iterator.size() ) );
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
                std::pair<map_id_struct, AMP::shared_ptr<std::vector<MeshElement>>> data(
                    map_id, elements );
                d_boundarySets.insert( data );
                new_boundary_ids.insert( boundary_id );
            }
        }
    }
    d_boundaryIdSets = std::vector<int>( new_boundary_ids.begin(), new_boundary_ids.end() );
    int *send_ptr    = nullptr;
    if ( d_boundaryIdSets.size() > 0 )
        send_ptr     = &d_boundaryIdSets[0];
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
        std::vector<std::vector<AMP::shared_ptr<std::vector<MeshElement>>>>( (int) GeomDim + 1 );
    for ( int t = 0; t <= (int) GeomDim; t++ ) {
        d_surface[t] = std::vector<AMP::shared_ptr<std::vector<MeshElement>>>( d_max_gcw + 1 );
        for ( int gcw = 0; gcw <= d_max_gcw; gcw++ ) {
            MeshIterator iterator1 = MultiVectorIterator( d_elements[t][gcw], 0 );
            MeshIterator iterator2 = d_parent_mesh->getSurfaceIterator( (GeomType) t, gcw );
            iterator               = Mesh::getIterator( Intersection, iterator1, iterator2 );
            AMP::shared_ptr<std::vector<MeshElement>> elements;
            if ( iterator.size() == 0 ) {
                elements =
                    AMP::shared_ptr<std::vector<MeshElement>>( new std::vector<MeshElement>( 0 ) );
            } else if ( iterator.size() == iterator1.size() ) {
                elements = d_elements[t][gcw];
            } else {
                elements = AMP::shared_ptr<std::vector<MeshElement>>(
                    new std::vector<MeshElement>( iterator.size() ) );
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
        iterator = getIterator( GeomDim, 0 );
        for ( size_t i = 0; i < iterator.size(); i++ ) {
            for ( auto &block_id : block_ids ) {
                if ( iterator->isInBlock( block_id ) )
                    new_block_ids.insert( block_id );
            }
            ++iterator;
        }
        d_blockIdSets = std::vector<int>( new_block_ids.begin(), new_block_ids.end() );
        send_ptr      = nullptr;
        if ( d_boundaryIdSets.size() > 0 )
            send_ptr = &d_boundaryIdSets[0];
        recv_size    = d_comm.sumReduce( d_blockIdSets.size() );
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
SubsetMesh::~SubsetMesh() {}


/********************************************************
* Copy the mesh                                         *
********************************************************/
AMP::shared_ptr<Mesh> SubsetMesh::copy() const
{
    AMP_ERROR( "Copy is not currently supported with SubsetMesh" );
    return AMP::shared_ptr<Mesh>();
}


/********************************************************
* Function to return the meshID composing the mesh      *
********************************************************/
std::vector<MeshID> SubsetMesh::getAllMeshIDs() const
{
    std::vector<MeshID> ids = d_parent_mesh->getAllMeshIDs();
    ids.push_back( d_meshID );
    AMP::Utilities::quicksort( ids );
    return ids;
}
std::vector<MeshID> SubsetMesh::getBaseMeshIDs() const { return d_parent_mesh->getBaseMeshIDs(); }
std::vector<MeshID> SubsetMesh::getLocalMeshIDs() const
{
    std::vector<MeshID> ids = d_parent_mesh->getLocalMeshIDs();
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
AMP::shared_ptr<Mesh> SubsetMesh::Subset( MeshID meshID ) const
{
    if ( d_meshID == meshID || d_parent_mesh->meshID() == meshID )
        return AMP::const_pointer_cast<Mesh>( shared_from_this() );
    else
        return AMP::shared_ptr<Mesh>();
}


/********************************************************
* Function to return the mesh with the given name       *
********************************************************/
AMP::shared_ptr<Mesh> SubsetMesh::Subset( std::string name ) const
{
    if ( d_name == name || d_parent_mesh->getName() == name )
        return AMP::const_pointer_cast<Mesh>( shared_from_this() );
    else
        return AMP::shared_ptr<Mesh>();
}


/********************************************************
* Mesh iterators                                        *
********************************************************/
MeshIterator SubsetMesh::getIterator( const GeomType type, const int gcw ) const
{
    int gcw2 = gcw;
    if ( gcw2 >= (int) d_elements[type].size() )
        gcw2 = (int) d_elements[type].size() - 1;
    if ( gcw2 == 0 )
        return MultiVectorIterator( d_elements[type][0], 0 );
    std::vector<AMP::shared_ptr<MeshIterator>> iterators( gcw2 + 1 );
    for ( int i = 0; i <= gcw2; i++ )
        iterators[i] =
            AMP::shared_ptr<MeshIterator>( new MultiVectorIterator( d_elements[type][i], 0 ) );
    return MultiIterator( iterators, 0 );
}
MeshIterator SubsetMesh::getSurfaceIterator( const GeomType type, const int gcw ) const
{
    if ( gcw == 0 )
        return MultiVectorIterator( d_surface[type][0], 0 );
    if ( gcw >= (int) d_surface[type].size() )
        AMP_ERROR( "Maximum ghost width exceeded" );
    std::vector<AMP::shared_ptr<MeshIterator>> iterators( gcw + 1 );
    for ( int i = 0; i <= gcw; i++ )
        iterators[i] =
            AMP::shared_ptr<MeshIterator>( new MultiVectorIterator( d_surface[type][i], 0 ) );
    return MultiIterator( iterators, 0 );
}
std::vector<int> SubsetMesh::getBoundaryIDs() const { return d_boundaryIdSets; }
MeshIterator
SubsetMesh::getBoundaryIDIterator( const GeomType type, const int id, const int gcw ) const
{
    std::vector<AMP::shared_ptr<MeshIterator>> iterators;
    iterators.reserve( gcw + 1 );
    for ( int i = 0; i <= gcw; i++ ) {
        map_id_struct map_id;
        map_id.id   = id;
        map_id.type = type;
        map_id.gcw  = i;
        auto map_it = d_boundarySets.find( map_id );
        if ( map_it == d_boundarySets.end() )
            continue;
        iterators.push_back(
            AMP::shared_ptr<MeshIterator>( new MultiVectorIterator( map_it->second, 0 ) ) );
    }
    if ( iterators.empty() )
        return MeshIterator();
    if ( iterators.size() == 1 )
        return *iterators[0];
    return MultiIterator( iterators, 0 );
}
std::vector<int> SubsetMesh::getBlockIDs() const { return d_blockIdSets; }
MeshIterator SubsetMesh::getBlockIDIterator( const GeomType, const int, const int ) const
{
    AMP_ERROR( "getBlockIDIterator is not implimented for SubsetMesh yet" );
    return MeshIterator();
}


/********************************************************
* Check if the element is a member of the mesh          *
********************************************************/
bool SubsetMesh::isMember( const MeshElementID &id ) const
{
    if ( !d_parent_mesh->isMember( id ) )
        return false;
    int type = static_cast<int>( id.type() );
    if ( type >= static_cast<int>( d_elements.size() ) )
        return false;
    for ( auto &elem : d_elements[type] ) {
        if ( elem == nullptr )
            continue;
        const std::vector<MeshElement> &elements = *( elem.get() );
        for ( auto &element : elements ) {
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
    return d_elements[type][0]->size();
}
size_t SubsetMesh::numGlobalElements( const GeomType type ) const { return N_global[type]; }
size_t SubsetMesh::numGhostElements( const GeomType type, int gcw ) const
{
    AMP_ASSERT( type <= GeomDim );
    if ( gcw == 0 )
        return 0;
    if ( gcw >= (int) d_elements[type].size() )
        AMP_ERROR( "Maximum ghost width exceeded" );
    return d_elements[type][gcw]->size();
}
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


} // Mesh namespace
} // AMP namespace
