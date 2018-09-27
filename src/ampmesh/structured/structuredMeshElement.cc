#include "AMP/ampmesh/structured/structuredMeshElement.h"
#include "AMP/ampmesh/MeshElement.h"
#include "AMP/utils/Utilities.h"


namespace AMP {
namespace Mesh {


// Function to evaluate the magnitude of a cross product in 3d
double cross3magnitude( double a[3], double b[3] )
{
    double v[3];
    v[0] = a[1] * b[2] - a[2] * b[1];
    v[1] = a[2] * b[0] - a[0] * b[2];
    v[2] = a[0] * b[1] - a[1] * b[0];
    return std::sqrt( v[0] * v[0] + v[1] * v[1] + v[2] * v[2] );
}
// Function to evaluate the dot produce of a vector anda cross product in 3d ( a . ( b X c ) )
double dot3cross( double a[3], double b[3], double c[3] )
{
    double v[3];
    v[0] = b[1] * c[2] - b[2] * c[1];
    v[1] = b[2] * c[0] - b[0] * c[2];
    v[2] = b[0] * c[1] - b[1] * c[0];
    return a[0] * v[0] + a[1] * v[1] + a[2] * v[2];
}


/********************************************************
 * Constructors                                          *
 ********************************************************/
structuredMeshElement::structuredMeshElement() { reset(); }
void structuredMeshElement::reset()
{
    typeID        = getTypeID();
    element       = nullptr;
    d_index       = BoxMesh::MeshElementIndex();
    d_globalID    = MeshElementID();
    d_meshType    = GeomType::null;
    d_physicalDim = 0;
}
structuredMeshElement::structuredMeshElement( const BoxMesh::MeshElementIndex &index,
                                              const AMP::Mesh::BoxMesh *mesh )
{
    reset( index, mesh );
}
void structuredMeshElement::reset( const BoxMesh::MeshElementIndex &index,
                                   const AMP::Mesh::BoxMesh *mesh )
{
    typeID        = getTypeID();
    d_mesh        = mesh;
    d_meshType    = d_mesh->getGeomType();
    d_physicalDim = d_mesh->getDim();
    AMP_ASSERT( static_cast<int>( d_meshType ) > 0 && static_cast<int>( d_meshType ) <= 3 );
    d_index    = index;
    d_globalID = d_mesh->convert( index );
}
structuredMeshElement::structuredMeshElement( const structuredMeshElement &rhs )
    : MeshElement(), // Note: we never want to call the base copy constructor
      d_meshType( rhs.d_meshType ),
      d_physicalDim( rhs.d_physicalDim ),
      d_index( rhs.d_index ),
      d_mesh( rhs.d_mesh )
{
    typeID     = getTypeID();
    element    = nullptr;
    d_globalID = rhs.d_globalID;
}
structuredMeshElement &structuredMeshElement::operator=( const structuredMeshElement &rhs )
{
    if ( this == &rhs ) // protect against invalid self-assignment
        return *this;
    this->typeID        = getTypeID();
    this->element       = nullptr;
    this->d_globalID    = rhs.d_globalID;
    this->d_meshType    = rhs.d_meshType;
    this->d_physicalDim = rhs.d_physicalDim;
    this->d_index       = rhs.d_index;
    this->d_mesh        = rhs.d_mesh;
    return *this;
}


/****************************************************************
 * De-constructor                                                *
 ****************************************************************/
structuredMeshElement::~structuredMeshElement() = default;


/****************************************************************
 * Function to clone the element                                 *
 ****************************************************************/
MeshElement *structuredMeshElement::clone() const { return new structuredMeshElement( *this ); }


/****************************************************************
 * Return the global rank of the owner rank                      *
 ****************************************************************/
unsigned int structuredMeshElement::globalOwnerRank() const
{
    return d_mesh->getComm().globalRanks()[d_globalID.owner_rank()];
}


/****************************************************************
 * Function to get the elements composing the current element    *
 * We use a Canonical numbering system                           *
 ****************************************************************/
void structuredMeshElement::getElements( const GeomType type,
                                         std::vector<MeshElement> &elements ) const
{
    int N = 0;
    BoxMesh::MeshElementIndex index[12];
    getElementIndex( type, N, index );
    // Create the elements
    elements.clear();
    elements.reserve( N );
    for ( int i = 0; i < N; i++ )
        elements.emplace_back( structuredMeshElement( index[i], d_mesh ) );
}
void structuredMeshElement::getElementsID( const GeomType type,
                                           std::vector<MeshElementID> &ID ) const
{
    int N = 0;
    BoxMesh::MeshElementIndex index[12];
    getElementIndex( type, N, index );
    ID.resize( N );
    for ( int i = 0; i < N; i++ )
        ID[i] = d_mesh->convert( index[i] );
}
void structuredMeshElement::getElementIndex( const GeomType type,
                                             int &N,
                                             BoxMesh::MeshElementIndex *index ) const
{
    AMP_ASSERT( type <= d_meshType );
    N = 0;
    if ( type == d_globalID.type() ) {
        N        = 1;
        index[0] = d_index;
        return;
    }
    const int *ijk = d_index.d_index;
    if ( type == GeomType::Vertex ) {
        // We want to get the verticies composing the elements
        if ( d_globalID.type() == d_meshType ) {
            // We are dealing with a entity of type dim and want the verticies
            if ( d_meshType == GeomType::Edge ) {
                N = 2;
                index[0].reset( GeomType::Vertex, 0, ijk[0] );
                index[1].reset( GeomType::Vertex, 0, ijk[0] + 1 );
            } else if ( d_meshType == GeomType::Face ) {
                N = 4;
                index[0].reset( GeomType::Vertex, 0, ijk[0], ijk[1] );
                index[1].reset( GeomType::Vertex, 0, ijk[0] + 1, ijk[1] );
                index[2].reset( GeomType::Vertex, 0, ijk[0] + 1, ijk[1] + 1 );
                index[3].reset( GeomType::Vertex, 0, ijk[0], ijk[1] + 1 );
            } else if ( d_meshType == GeomType::Volume ) {
                N = 8;
                index[0].reset( GeomType::Vertex, 0, ijk[0], ijk[1], ijk[2] );
                index[1].reset( GeomType::Vertex, 0, ijk[0] + 1, ijk[1], ijk[2] );
                index[2].reset( GeomType::Vertex, 0, ijk[0] + 1, ijk[1] + 1, ijk[2] );
                index[3].reset( GeomType::Vertex, 0, ijk[0], ijk[1] + 1, ijk[2] );
                index[4].reset( GeomType::Vertex, 0, ijk[0], ijk[1], ijk[2] + 1 );
                index[5].reset( GeomType::Vertex, 0, ijk[0] + 1, ijk[1], ijk[2] + 1 );
                index[6].reset( GeomType::Vertex, 0, ijk[0] + 1, ijk[1] + 1, ijk[2] + 1 );
                index[7].reset( GeomType::Vertex, 0, ijk[0], ijk[1] + 1, ijk[2] + 1 );
            } else {
                AMP_ERROR( "Dimension not supported yet" );
            }
        } else if ( d_globalID.type() == GeomType::Edge ) {
            N                                = 2;
            index[0]                         = d_index;
            index[1]                         = d_index;
            index[0].d_type                  = 0;
            index[1].d_type                  = 0;
            index[0].d_side                  = 0;
            index[1].d_side                  = 0;
            index[0].d_index[d_index.d_side] = ijk[d_index.d_side];
            index[1].d_index[d_index.d_side] = ijk[d_index.d_side] + 1;
        } else if ( d_globalID.type() == GeomType::Face ) {
            N = 4;
            if ( d_index.d_side == 0 ) {
                index[0].reset( GeomType::Vertex, 0, ijk[0], ijk[1], ijk[2] );
                index[1].reset( GeomType::Vertex, 0, ijk[0], ijk[1] + 1, ijk[2] );
                index[2].reset( GeomType::Vertex, 0, ijk[0], ijk[1] + 1, ijk[2] + 1 );
                index[3].reset( GeomType::Vertex, 0, ijk[0], ijk[1], ijk[2] + 1 );
            } else if ( d_index.d_side == 1 ) {
                index[0].reset( GeomType::Vertex, 0, ijk[0], ijk[1], ijk[2] );
                index[1].reset( GeomType::Vertex, 0, ijk[0] + 1, ijk[1], ijk[2] );
                index[2].reset( GeomType::Vertex, 0, ijk[0] + 1, ijk[1], ijk[2] + 1 );
                index[3].reset( GeomType::Vertex, 0, ijk[0], ijk[1], ijk[2] + 1 );
            } else if ( d_index.d_side == 2 ) {
                index[0].reset( GeomType::Vertex, 0, ijk[0], ijk[1], ijk[2] );
                index[1].reset( GeomType::Vertex, 0, ijk[0] + 1, ijk[1], ijk[2] );
                index[2].reset( GeomType::Vertex, 0, ijk[0] + 1, ijk[1] + 1, ijk[2] );
                index[3].reset( GeomType::Vertex, 0, ijk[0], ijk[1] + 1, ijk[2] );
            } else {
                AMP_ERROR( "Internal error" );
            }
        } else if ( d_globalID.type() == GeomType::Volume ) {
            AMP_ERROR( "Not ready for dimensions > 3" );
        } else {
            AMP_ERROR( "Not finsihed" );
        }
    } else if ( type == GeomType::Edge ) {
        if ( d_meshType == GeomType::Face ) {
            N = 4;
            index[0].reset( GeomType::Edge, 0, ijk[0], ijk[1], 0 );
            index[1].reset( GeomType::Edge, 1, ijk[0] + 1, ijk[1], 0 );
            index[2].reset( GeomType::Edge, 0, ijk[0], ijk[1] + 1, 0 );
            index[3].reset( GeomType::Edge, 1, ijk[0], ijk[1], 0 );
        } else if ( d_meshType == GeomType::Volume ) {
            if ( d_globalID.type() == GeomType::Face ) {
                N = 4;
                if ( d_index.d_side == 0 ) {
                    // We are dealing with an x-face
                    index[0].reset( GeomType::Edge, 1, ijk[0], ijk[1], ijk[2] );
                    index[1].reset( GeomType::Edge, 2, ijk[0], ijk[1] + 1, ijk[2] );
                    index[2].reset( GeomType::Edge, 1, ijk[0], ijk[1], ijk[2] + 1 );
                    index[3].reset( GeomType::Edge, 2, ijk[0], ijk[1], ijk[2] );
                } else if ( d_index.d_side == 1 ) {
                    // We are dealing with an y-face
                    index[0].reset( GeomType::Edge, 0, ijk[0], ijk[1], ijk[2] );
                    index[1].reset( GeomType::Edge, 2, ijk[0] + 1, ijk[1], ijk[2] );
                    index[2].reset( GeomType::Edge, 0, ijk[0], ijk[1], ijk[2] + 1 );
                    index[3].reset( GeomType::Edge, 2, ijk[0], ijk[1], ijk[2] );
                } else if ( d_index.d_side == 2 ) {
                    // We are dealing with an z-face
                    index[0].reset( GeomType::Edge, 0, ijk[0], ijk[1], ijk[2] );
                    index[1].reset( GeomType::Edge, 1, ijk[0] + 1, ijk[1], ijk[2] );
                    index[2].reset( GeomType::Edge, 0, ijk[0], ijk[1] + 1, ijk[2] );
                    index[3].reset( GeomType::Edge, 1, ijk[0], ijk[1], ijk[2] );
                } else {
                    AMP_ERROR( "Internal error" );
                }
            } else if ( d_globalID.type() == GeomType::Volume ) {
                AMP_ASSERT( d_index.d_side == 0 );
                N = 12;
                index[0].reset( GeomType::Edge, 0, ijk[0], ijk[1], ijk[2] );
                index[1].reset( GeomType::Edge, 1, ijk[0] + 1, ijk[1], ijk[2] );
                index[2].reset( GeomType::Edge, 0, ijk[0], ijk[1] + 1, ijk[2] );
                index[3].reset( GeomType::Edge, 1, ijk[0], ijk[1], ijk[2] );
                index[4].reset( GeomType::Edge, 2, ijk[0], ijk[1], ijk[2] );
                index[5].reset( GeomType::Edge, 2, ijk[0] + 1, ijk[1], ijk[2] );
                index[6].reset( GeomType::Edge, 2, ijk[0] + 1, ijk[1] + 1, ijk[2] );
                index[7].reset( GeomType::Edge, 2, ijk[0], ijk[1] + 1, ijk[2] );
                index[8].reset( GeomType::Edge, 0, ijk[0], ijk[1], ijk[2] + 1 );
                index[9].reset( GeomType::Edge, 1, ijk[0] + 1, ijk[1], ijk[2] + 1 );
                index[10].reset( GeomType::Edge, 0, ijk[0], ijk[1] + 1, ijk[2] + 1 );
                index[11].reset( GeomType::Edge, 1, ijk[0], ijk[1], ijk[2] + 1 );
            } else {
                AMP_ERROR( "Dimensions > 3 are not supported yet" );
            }
        } else {
            AMP_ERROR( "Dimensions > 3 are not supported yet" );
        }
    } else if ( type == GeomType::Face ) {
        if ( d_globalID.type() == GeomType::Volume ) {
            N = 6;
            index[0].reset( GeomType::Face, 1, ijk[0], ijk[1], ijk[2] );
            index[1].reset( GeomType::Face, 0, ijk[0] + 1, ijk[1], ijk[2] );
            index[2].reset( GeomType::Face, 1, ijk[0], ijk[1] + 1, ijk[2] );
            index[3].reset( GeomType::Face, 0, ijk[0], ijk[1], ijk[2] );
            index[4].reset( GeomType::Face, 2, ijk[0], ijk[1], ijk[2] );
            index[5].reset( GeomType::Face, 2, ijk[0], ijk[1], ijk[2] + 1 );
        } else {
            AMP_ERROR( "Dimensions > 3 are not supported yet" );
        }
    } else if ( type == GeomType::Volume ) {
        AMP_ERROR( "Dimensions > 3 are not supported yet" );
    } else {
        AMP_ERROR( "Not finished" );
    }
    // Fix any elements that are beyond a periodic boundary
    for ( int d = 0; d < static_cast<int>( d_meshType ); d++ ) {
        if ( d_mesh->d_isPeriodic[d] ) {
            int size = d_mesh->d_globalSize[d];
            for ( int i = 0; i < N; i++ ) {
                if ( index[i].d_index[d] < 0 )
                    index[i].d_index[d] += size;
                else if ( index[i].d_index[d] >= size )
                    index[i].d_index[d] -= size;
            }
        }
    }
}


/****************************************************************
 * Function to get the neighboring elements                      *
 ****************************************************************/
void structuredMeshElement::getNeighbors( std::vector<MeshElement::shared_ptr> &neighbors ) const
{
    int N = 0;
    BoxMesh::MeshElementIndex index[27];
    getNeighborIndex( N, index );
    // Get the neighbor elements
    neighbors.resize( N );
    bool periodic[3];
    int size[3];
    for ( int d = 0; d < static_cast<int>( d_meshType ); d++ ) {
        periodic[d] = d_mesh->d_isPeriodic[d];
        size[d]     = d_mesh->d_globalSize[d];
    }
    for ( int i = 0; i < N; i++ ) {
        bool in_mesh = true;
        auto &elem   = index[i];
        for ( int d = 0; d < static_cast<int>( d_meshType ); d++ ) {
            if ( periodic[d] ) {
                if ( elem.d_index[d] < 0 )
                    elem.d_index[d] += size[d];
                if ( elem.d_index[d] >= size[d] )
                    elem.d_index[d] -= size[d];
            } else {
                if ( elem.d_index[d] < 0 )
                    in_mesh = false;
                if ( d_globalID.type() == d_meshType ) {
                    if ( elem.d_index[d] >= size[d] )
                        in_mesh = false;
                } else {
                    if ( elem.d_index[d] > size[d] )
                        in_mesh = false;
                }
            }
        }
        if ( in_mesh )
            neighbors[i].reset( new structuredMeshElement( elem, d_mesh ) );
        else if ( d_globalID.type() != GeomType::Vertex )
            neighbors[i].reset();
    }
}
void structuredMeshElement::getNeighborIndex( int &N, BoxMesh::MeshElementIndex *index ) const
{
    const int *ijk = d_index.d_index;
    if ( d_globalID.type() == GeomType::Vertex ) {
        // Get the list of neighbor nodex (there are no null neighbors)
        // The node neighbors are the list of nodes that share any element
        if ( d_meshType == GeomType::Edge ) {
            N = 2;
            index[0].reset( GeomType::Vertex, 0, ijk[0] - 1 );
            index[1].reset( GeomType::Vertex, 0, ijk[0] + 1 );
        } else if ( d_meshType == GeomType::Face ) {
            N = 8;
            index[0].reset( GeomType::Vertex, 0, ijk[0] - 1, ijk[1] - 1 );
            index[1].reset( GeomType::Vertex, 0, ijk[0], ijk[1] - 1 );
            index[2].reset( GeomType::Vertex, 0, ijk[0] + 1, ijk[1] - 1 );
            index[3].reset( GeomType::Vertex, 0, ijk[0] - 1, ijk[1] );
            index[4].reset( GeomType::Vertex, 0, ijk[0] + 1, ijk[1] );
            index[5].reset( GeomType::Vertex, 0, ijk[0] - 1, ijk[1] + 1 );
            index[6].reset( GeomType::Vertex, 0, ijk[0], ijk[1] + 1 );
            index[7].reset( GeomType::Vertex, 0, ijk[0] + 1, ijk[1] + 1 );
        } else if ( d_meshType == GeomType::Volume ) {
            N = 0;
            for ( int k = -1; k <= 1; k++ ) {
                for ( int j = -1; j <= 1; j++ ) {
                    for ( int i = -1; i <= 1; i++ ) {
                        if ( i == 0 && j == 0 && k == 0 )
                            continue;
                        index[N].reset( GeomType::Vertex, 0, ijk[0] + i, ijk[1] + j, ijk[2] + k );
                        N++;
                    }
                }
            }
        } else {
            AMP_ERROR( "Dimension not supported yet" );
        }
    } else if ( d_globalID.type() == GeomType::Edge ) {
        if ( d_meshType == GeomType::Edge ) {
            N = 2;
            index[0].reset( GeomType::Edge, 0, ijk[0] - 1 );
            index[1].reset( GeomType::Edge, 0, ijk[0] + 1 );
        } else {
            // GeomType::Edge neighbors in dimensions > 1 are not supported yet
        }
    } else if ( d_globalID.type() == GeomType::Face ) {
        if ( d_meshType == GeomType::Face ) {
            N = 4;
            index[0].reset( GeomType::Face, 0, ijk[0], ijk[1] - 1 );
            index[1].reset( GeomType::Face, 0, ijk[0] + 1, ijk[0] );
            index[2].reset( GeomType::Face, 0, ijk[0], ijk[1] + 1 );
            index[3].reset( GeomType::Face, 0, ijk[0] - 1, ijk[1] );
        } else {
            // GeomType::Face neighbors in dimensions > 2 are not supported yet
        }
    } else if ( d_globalID.type() == GeomType::Volume ) {
        if ( d_meshType == GeomType::Volume ) {
            N = 6;
            index[0].reset( GeomType::Volume, 0, ijk[0], ijk[1] - 1, ijk[2] );
            index[1].reset( GeomType::Volume, 0, ijk[0] + 1, ijk[1], ijk[2] );
            index[2].reset( GeomType::Volume, 0, ijk[0], ijk[1] + 1, ijk[2] );
            index[3].reset( GeomType::Volume, 0, ijk[0] - 1, ijk[1], ijk[2] );
            index[4].reset( GeomType::Volume, 0, ijk[0], ijk[1], ijk[2] - 1 );
            index[5].reset( GeomType::Volume, 0, ijk[0], ijk[1], ijk[2] + 1 );
        } else {
            // GeomType::Volume neighbors in dimensions > 3 are not supported yet
        }
    } else {
        AMP_ERROR( "Unknown entity type" );
    }
}


/****************************************************************
 * Function to get the parent elements                           *
 ****************************************************************/
std::vector<MeshElement> structuredMeshElement::getParents( GeomType type ) const
{
    AMP_INSIST( static_cast<int>( type ) >= d_index.d_type,
                "We can't return a parent of geometric type < current type" );
    // Get the indicies of the parent elements (ignore boundaries for now)
    std::vector<BoxMesh::MeshElementIndex> index_list;
    const int *ijk = d_index.d_index;
    if ( d_index.d_type == static_cast<int>( type ) ) {
        // We are looking for the current element
        return std::vector<MeshElement>( 1, MeshElement( *this ) );
    } else if ( static_cast<int>( type ) == d_index.d_type + 1 && type == d_mesh->getGeomType() ) {
        // We have an entity that is the geometric type-1 and we want to get the parents of the
        // geometric type of the
        // mesh
        BoxMesh::MeshElementIndex index( type, 0, ijk[0], ijk[1], ijk[2] );
        index_list.emplace_back( index );
        index.d_index[d_index.d_side]--;
        index_list.emplace_back( index );
    } else if ( d_index.d_type == static_cast<int>( GeomType::Vertex ) ) {
        // We want to get the parents of a vertex
        AMP_ASSERT( static_cast<int>( d_meshType ) <= 3 );
        if ( type == d_mesh->getGeomType() ) {
            for ( int i = ijk[0] - 1; i <= ijk[0]; i++ ) {
                for ( int j = ijk[1] - 1; j <= ijk[1]; j++ ) {
                    for ( int k = ijk[2] - 1; k <= ijk[2]; k++ ) {
                        index_list.emplace_back( BoxMesh::MeshElementIndex( type, 0, i, j, k ) );
                    }
                }
            }
        } else if ( type == GeomType::Edge ) {
            for ( int d = 0; d < static_cast<int>( d_meshType ); d++ ) {
                BoxMesh::MeshElementIndex index( type, d, ijk[0], ijk[1], ijk[2] );
                index_list.emplace_back( index );
                index.d_index[d]--;
                index_list.emplace_back( index );
            }
        } else if ( type == GeomType::Face && d_mesh->getGeomType() == GeomType::Volume ) {
            index_list.resize( 12 );
            index_list[0].reset( type, 0, ijk[0], ijk[1] - 1, ijk[2] - 1 );
            index_list[1].reset( type, 0, ijk[0], ijk[1] - 1, ijk[2] );
            index_list[2].reset( type, 0, ijk[0], ijk[1], ijk[2] - 1 );
            index_list[3].reset( type, 0, ijk[0], ijk[1], ijk[2] );
            index_list[4].reset( type, 1, ijk[0] - 1, ijk[1], ijk[2] - 1 );
            index_list[5].reset( type, 1, ijk[0] - 1, ijk[1], ijk[2] );
            index_list[6].reset( type, 1, ijk[0], ijk[1], ijk[2] - 1 );
            index_list[7].reset( type, 1, ijk[0], ijk[1], ijk[2] );
            index_list[8].reset( type, 2, ijk[0] - 1, ijk[1] - 1, ijk[2] );
            index_list[9].reset( type, 2, ijk[0] - 1, ijk[1], ijk[2] );
            index_list[10].reset( type, 2, ijk[0], ijk[1] - 1, ijk[2] );
            index_list[11].reset( type, 2, ijk[0], ijk[1], ijk[2] );
        } else {
            char text[100];
            sprintf( text,
                     "Unknown type: dim=%i, elem_type=%i, type=%i",
                     (int) d_meshType,
                     (int) d_index.d_type,
                     (int) type );
            AMP_ERROR( std::string( text ) );
        }
    } else if ( d_index.d_type == static_cast<int>( GeomType::Edge ) ) {
        // We want to get the parents of an edge
        AMP_ASSERT( static_cast<int>( d_meshType ) <= 3 );
        int i = ijk[0];
        int j = ijk[1];
        int k = ijk[2];
        if ( type == GeomType::Face && d_mesh->getGeomType() == GeomType::Volume ) {
            if ( d_index.d_side == 0 ) {
                index_list.emplace_back( type, 2, i, j - 1, k );
                index_list.emplace_back( type, 2, i, j, k );
                index_list.emplace_back( type, 1, i, j, k - 1 );
                index_list.emplace_back( type, 1, i, j, k );
            } else if ( d_index.d_side == 1 ) {
                index_list.emplace_back( type, 2, i - 1, j, k );
                index_list.emplace_back( type, 2, i, j, k );
                index_list.emplace_back( type, 0, i, j, k - 1 );
                index_list.emplace_back( type, 0, i, j, k );
            } else if ( d_index.d_side == 2 ) {
                index_list.emplace_back( type, 1, i - 1, j, k );
                index_list.emplace_back( type, 1, i, j, k );
                index_list.emplace_back( type, 0, i, j - 1, k );
                index_list.emplace_back( type, 0, i, j, k );
            } else {
                AMP_ERROR( "Internal error" );
            }
        } else if ( type == GeomType::Volume && d_mesh->getGeomType() == GeomType::Volume ) {
            if ( d_index.d_side == 0 ) {
                index_list.emplace_back( type, 0, i, j - 1, k - 1 );
                index_list.emplace_back( type, 0, i, j, k - 1 );
                index_list.emplace_back( type, 0, i, j - 1, k );
                index_list.emplace_back( type, 0, i, j, k );
            } else if ( d_index.d_side == 1 ) {
                index_list.emplace_back( type, 0, i - 1, j, k - 1 );
                index_list.emplace_back( type, 0, i, j, k - 1 );
                index_list.emplace_back( type, 0, i - 1, j, k );
                index_list.emplace_back( type, 0, i, j, k );
            } else if ( d_index.d_side == 2 ) {
                index_list.emplace_back( type, 0, i - 1, j - 1, k );
                index_list.emplace_back( type, 0, i, j - 1, k );
                index_list.emplace_back( type, 0, i - 1, j, k );
                index_list.emplace_back( type, 0, i, j, k );
            } else {
                AMP_ERROR( "Internal error" );
            }
        } else {
            char text[100];
            sprintf( text,
                     "Unknown type: dim=%i, elem_type=%i, type=%i",
                     (int) d_meshType,
                     (int) d_index.d_type,
                     (int) type );
            AMP_ERROR( std::string( text ) );
        }
    } else {
        char text[100];
        sprintf( text,
                 "Case not programmed yet: dim=%i, elem_type=%i, type=%i",
                 (int) d_meshType,
                 (int) d_index.d_type,
                 (int) type );
        AMP_ERROR( std::string( text ) );
    }
    // Get some basic properties from the mesh
    auto meshGeomDim = (int) d_mesh->getGeomType();
    bool periodic[3] = { false, false, false };
    for ( int d = 0; d < meshGeomDim; d++ )
        periodic[d] = d_mesh->d_isPeriodic[d];
    int size[3] = { 1, 1, 1 };
    for ( int d = 0; d < meshGeomDim; d++ )
        size[d] = d_mesh->d_globalSize[d];
    // Remove any elements that are outside the physical domain
    size_t k = 0;
    for ( size_t i = 0; i < index_list.size(); i++ ) {
        bool erase = false;
        if ( static_cast<int>( d_meshType ) < 2 && index_list[i].d_index[1] != 0 )
            erase = true;
        if ( static_cast<int>( d_meshType ) < 3 && index_list[i].d_index[2] != 0 )
            erase = true;
        for ( int d = 0; d < meshGeomDim; d++ ) {
            if ( periodic[d] )
                continue;
            int i_max = size[d];
            if ( (int) type == meshGeomDim ) {
                // The geometric type of the mesh must be within the domain
            } else if ( type == GeomType::Vertex ) {
                // Verticies can exist on all faces
                i_max++;
            } else if ( (int) type == meshGeomDim - 1 ) {
                if ( index_list[i].d_side == d )
                    i_max++;
            } else {
                if ( index_list[i].d_side != d )
                    i_max++;
            }
            if ( index_list[i].d_index[d] < 0 || index_list[i].d_index[d] >= i_max )
                erase = true;
        }
        if ( !erase ) {
            index_list[k] = index_list[i];
            k++;
        }
    }
    index_list.resize( k );
    // Fix any elements that are beyond a periodic boundary
    for ( int d = 0; d < static_cast<int>( d_meshType ); d++ ) {
        if ( d_mesh->d_isPeriodic[d] ) {
            int size = d_mesh->d_globalSize[d];
            for ( auto &elem : index_list ) {
                if ( elem.d_index[d] < 0 )
                    elem.d_index[d] += size;
                else if ( elem.d_index[d] >= size )
                    elem.d_index[d] -= size;
            }
        }
    }
    // Create the elements
    AMP::Utilities::quicksort( index_list );
    std::vector<MeshElement> elements( index_list.size() );
    for ( size_t i = 0; i < index_list.size(); i++ )
        elements[i] = structuredMeshElement( index_list[i], d_mesh );
    return elements;
}


/****************************************************************
 * Functions to get the element volume                           *
 ****************************************************************/
double structuredMeshElement::volume() const
{
    if ( d_globalID.type() == GeomType::Vertex ) {
        AMP_ERROR( "volume is is not defined Nodes" );
    }
    int N = 0;
    BoxMesh::MeshElementIndex nodes[8];
    getElementIndex( GeomType::Vertex, N, nodes );
    if ( d_globalID.type() == GeomType::Edge ) {
        AMP_ASSERT( N == 2 );
        double x[2][3];
        d_mesh->coord( nodes[0], x[0] );
        d_mesh->coord( nodes[1], x[1] );
        double dist2 = 0.0;
        for ( int i = 0; i < d_physicalDim; i++ )
            dist2 += ( x[0][i] - x[1][i] ) * ( x[0][i] - x[1][i] );
        return sqrt( dist2 );
    } else if ( d_globalID.type() == GeomType::Face ) {
        // Use 2x2 quadrature to approximate the surface area. See for example,
        // Y. Zhang, C. Bajaj, G. Xu. Surface Smoothing and Quality Improvement
        // of Quadrilateral/Hexahedral Meshes with Geometric Flow. The special
        // issue of the Journal Communications in Numerical Methods in
        // Engineering (CNME), submitted as an invited paper, 2006.
        // http://www.ices.utexas.edu/~jessica/paper/quadhexgf/quadhex_geomflow_CNM.pdf
        AMP_ASSERT( N == 4 );
        double x[4][3];
        d_mesh->coord( nodes[0], x[0] );
        d_mesh->coord( nodes[1], x[1] );
        d_mesh->coord( nodes[2], x[2] );
        d_mesh->coord( nodes[3], x[3] );
        double AB[3] = { 0, 0, 0 }, AC[3] = { 0, 0, 0 }, AD[3] = { 0, 0, 0 },
               AC_AB_AD[3] = { 0, 0, 0 };
        for ( int i = 0; i < d_physicalDim; i++ ) {
            AC[i]       = x[2][i] - x[0][i];     // Vector pointing from A to C
            AB[i]       = x[1][i] - x[0][i];     // Vector pointing from A to B
            AD[i]       = x[3][i] - x[0][i];     // Vector pointing from A to D
            AC_AB_AD[i] = AC[i] - AB[i] - AD[i]; // The diagonal vector minus the side vectors
        }
        if ( AC_AB_AD[0] == 0 && AC_AB_AD[1] == 0 && AC_AB_AD[2] == 0 ) {
            // The points are co-planar
            return cross3magnitude( AB, AD );
        } else {
            const double q[2] = { 0.5 - std::sqrt( 3.0 ) / 6.0, 0.5 + std::sqrt( 3.0 ) / 6.0 };
            double vol        = 0.0;
            double v1[3], v2[3];
            for ( auto &elem : q ) {
                for ( auto &q_j : q ) {
                    v1[0] = AB[0] + elem * AC_AB_AD[0];
                    v1[1] = AB[1] + elem * AC_AB_AD[1];
                    v1[2] = AB[2] + elem * AC_AB_AD[2];
                    v2[0] = AD[0] + q_j * AC_AB_AD[0];
                    v2[1] = AD[1] + q_j * AC_AB_AD[1];
                    v2[2] = AD[2] + q_j * AC_AB_AD[2];
                    vol += cross3magnitude( v1, v2 );
                }
            }
            return 0.25 * vol;
        }
    } else if ( d_globalID.type() == GeomType::Volume ) {
        // Compute the volume of the tri-linear hex by splitting it
        // into 6 sub-pyramids and applying the formula in:
        // "Calculation of the GeomType::Volume of a General Hexahedron
        // for Flow Predictions", AIAA Journal v.23, no.6, 1984, p.954-
        static const unsigned char sub_pyr[6][4] = {
            { 0, 3, 2, 1 }, { 6, 7, 4, 5 }, { 0, 1, 5, 4 },
            { 3, 7, 6, 2 }, { 0, 4, 7, 3 }, { 1, 2, 6, 5 }
        };
        // The centroid is a convenient point to use
        // for the apex of all the pyramids.
        std::vector<double> R = this->centroid();
        AMP_ASSERT( N == 8 );
        double x[8][3];
        for ( int i = 0; i < 8; i++ )
            d_mesh->coord( nodes[i], x[i] );
        int pyr_base[4];
        double vol = 0.0;
        // Compute the volume using 6 sub-pyramids
        for ( auto &elem : sub_pyr ) {
            // Set the nodes of the pyramid base
            for ( unsigned int i = 0; i < 4; ++i )
                pyr_base[i] = elem[i];
            // Compute diff vectors
            double a[3], b[3], c[3], d[3], e[3];
            for ( int i = 0; i < 3; i++ ) {
                a[i] = x[pyr_base[0]][i] - R[i];
                b[i] = x[pyr_base[1]][i] - x[pyr_base[3]][i];
                c[i] = x[pyr_base[2]][i] - x[pyr_base[0]][i];
                d[i] = x[pyr_base[3]][i] - x[pyr_base[0]][i];
                e[i] = x[pyr_base[1]][i] - x[pyr_base[0]][i];
            }
            // Compute pyramid volume
            double sub_vol =
                ( 1.0 / 6.0 ) * dot3cross( a, b, c ) + ( 1.0 / 12.0 ) * dot3cross( c, d, e );
            AMP_ASSERT( sub_vol > 0.0 );
            vol += sub_vol;
        }
        return vol;
    }
    AMP_ERROR( "Internal error" );
    return 0.0;
}


/****************************************************************
 * Misc functions                                                *
 ****************************************************************/
bool structuredMeshElement::containsPoint( const std::vector<double> &, double ) const
{
    AMP_ERROR( "Not finsihed" );
    return false;
}
bool structuredMeshElement::isOnSurface() const
{
    bool on_surface = false;
    const int *ijk  = d_index.d_index;
    for ( int d = 0; d < static_cast<int>( d_meshType ); d++ ) {
        if ( d_mesh->d_isPeriodic[d] )
            continue;
        auto size = (int) d_mesh->d_globalSize[d];
        if ( d_globalID.type() == d_mesh->GeomDim ) {
            // We are dealing with the highest level geometric entity
            if ( ijk[d] == 0 || ijk[d] == size - 1 )
                on_surface = true;
        } else if ( d_globalID.type() == GeomType::Vertex ) {
            // We are dealing with a vertex
            if ( ijk[d] == 0 || ijk[d] == size )
                on_surface = true;
        } else if ( d_globalID.type() == GeomType::Edge ) {
            // We are dealing with a vertex
            if ( ( ijk[d] == 0 || ijk[d] == size ) && d_index.d_side != d )
                on_surface = true;
        } else if ( d_globalID.type() == GeomType::Face ) {
            // We are dealing with a vertex
            if ( ( ijk[d] == 0 || ijk[d] == size ) && d_index.d_side == d )
                on_surface = true;
        } else {
            AMP_ERROR( "Internal error (dim>3?)" );
        }
    }
    return on_surface;
}
bool structuredMeshElement::isOnBoundary( int id ) const
{
    return d_mesh->isOnBoundary( d_index, id );
}
bool structuredMeshElement::isInBlock( int id ) const
{
    if ( id == 0 )
        return true;
    return false;
}


} // namespace Mesh
} // namespace AMP
