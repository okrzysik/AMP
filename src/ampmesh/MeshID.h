#ifndef included_AMP_MeshID
#define included_AMP_MeshID

#include "AMP/utils/shared_ptr.h"
#include <vector>


namespace AMP {
namespace Mesh {


//! Enumeration for basic mesh-based quantities
enum class GeomType : unsigned char { Vertex = 0, Edge = 1, Face = 2, Volume = 3, null = 0xFF };


struct MeshElementID;


/**
 * \struct MeshID
 * \brief  A structure used to identify the mesh
 * \details  This structure provides a unique id that can be used to identify a mesh.
 */
struct MeshID {

public:
    // Constructors used to initialize key values
    constexpr MeshID() : data( 0xFFFFFFFFFFFFFFFF ) {}
    constexpr MeshID( unsigned int root, unsigned int local_id ):
        data( ( static_cast<uint64_t>( root ) << 32 ) + local_id ) {}
    constexpr MeshID( uint64_t id ) : data( id ) {}
    constexpr uint64_t getData() const { return data; }
    // Overload key operators
    constexpr bool operator==( const MeshID &rhs ) const { return data == rhs.data; }
    constexpr bool operator!=( const MeshID &rhs ) const { return data != rhs.data; }
    constexpr bool operator>=( const MeshID &rhs ) const { return data >= rhs.data; }
    constexpr bool operator<=( const MeshID &rhs ) const { return data <= rhs.data; }
    constexpr bool operator>( const MeshID &rhs ) const { return data > rhs.data; }
    constexpr bool operator<( const MeshID &rhs ) const { return data < rhs.data; }

private:
    // We will store the data as a 64-bit data type
    // We do not want the user to get direct access to the data
    uint64_t data;
    friend struct MeshElementID;
};


/**
 * \struct MeshElementID
 * \brief  A structure used to identify the mesh element
 * \details  This structure provides a unique id that can be used to identify a mesh element.
 */
struct MeshElementID {

public:
    // Constructors used to initialize key values
    constexpr MeshElementID(): meshId( 0xFFFFFFFFFFFFFFFF ), elemId( 0x000000FFFFFFFFFF ) {}
    constexpr explicit MeshElementID( bool isLocal,
                            GeomType type_id,
                            unsigned int local_ID,
                            unsigned int owner_rank_id,
                            MeshID mesh_ID ):
        meshId( mesh_ID ), elemId( 0 )
    {
        // Set the bit for is_local
        unsigned int tmp = 0x00000000;
        if ( isLocal )
            tmp = 0x80000000;
        // Add the owner_rank
        tmp += ( 0x007FFFFF & owner_rank_id ) << 8;
        // Add the type_id
        char type = (char) type_id;
        tmp += ( (unsigned char) type );
        // Combine the above data with the local_ID
        elemId = ( ( (uint64_t) tmp ) << 32 ) + ( (uint64_t) local_ID );
    }
    constexpr explicit MeshElementID( MeshID mesh_ID, uint64_t elem_id ) : meshId( mesh_ID ), elemId( elem_id ) {}
    constexpr void resetElemID( uint64_t elem_id ) { elemId = elem_id; }
    // Overload key operators
    constexpr bool operator==( const MeshElementID &rhs ) const
    {
        uint64_t test = ( elemId ^ rhs.elemId ) << 1;
        return meshId == rhs.meshId && test == 0;
    }
    constexpr bool operator!=( const MeshElementID &rhs ) const
    {
        uint64_t test = ( elemId ^ rhs.elemId ) << 1;
        return meshId != rhs.meshId || test != 0;
    }
    constexpr bool operator>=( const MeshElementID &rhs ) const
    {
        // Sort by meshID first
        if ( meshId < rhs.meshId )
            return false;
        else if ( meshId > rhs.meshId )
            return true;
        // Sort by the remaining data
        uint64_t d1 = elemId << 1;
        uint64_t d2 = rhs.elemId << 1;
        return d1 >= d2;
    }
    constexpr bool operator>( const MeshElementID &rhs ) const
    {
        // Sort by meshID first
        if ( meshId < rhs.meshId )
            return false;
        else if ( meshId > rhs.meshId )
            return true;
        // Sort by the remaining data
        uint64_t d1 = elemId << 1;
        uint64_t d2 = rhs.elemId << 1;
        return d1 > d2;
    }
    constexpr bool operator<( const MeshElementID &rhs ) const
    {
        // Sort by meshID first
        if ( meshId > rhs.meshId )
            return false;
        else if ( meshId < rhs.meshId )
            return true;
        // Sort by the remaining data
        uint64_t d1 = elemId << 1;
        uint64_t d2 = rhs.elemId << 1;
        return d1 < d2;
    }
    constexpr bool operator<=( const MeshElementID &rhs ) const
    {
        // Sort by meshID first
        if ( meshId > rhs.meshId )
            return false;
        else if ( meshId < rhs.meshId )
            return true;
        // Sort by the remaining data
        uint64_t d1 = elemId << 1;
        uint64_t d2 = rhs.elemId << 1;
        return d1 <= d2;
    }
    // Access local data
    constexpr bool is_local() const { return ( elemId >> 63 ) != 0; }
    constexpr GeomType type() const
    {
        unsigned char tmp = (unsigned char) ( ( elemId >> 32 ) & 0x00FF );
        return (GeomType) tmp;
    }
    constexpr unsigned int local_id() const
    {
        unsigned int tmp = (unsigned int) ( elemId & 0x00000000FFFFFFFF );
        return tmp;
    }
    constexpr unsigned int owner_rank() const
    {
        unsigned int tmp = elemId >> 32;
        tmp              = ( tmp >> 8 ) & 0x007FFFFF;
        return tmp;
    }
    constexpr MeshID meshID() const { return meshId; }
    constexpr uint64_t elemID() const { return elemId; }
    constexpr void set_is_local( bool isLocal )
    {
        if ( isLocal )
            elemId |= 0x8000000000000000;
        else
            elemId &= 0x7FFFFFFFFFFFFFFF;
    }
    constexpr bool isNull() const
    {
        return meshId == 0xFFFFFFFFFFFFFFFF || elemId == 0x000000FFFFFFFFFF;
    }

private:
    // We will store the data as a 128-bit data type
    // The first 64 bits refer to the meshID
    // The next  1 bits refer to if the element is local
    // The next 23 bits refer to the processor id
    // The next  8 bits refer to the element type
    // The next 32 bits refer to the local id
    MeshID meshId;
    uint64_t elemId;
};


} // namespace Mesh
} // namespace AMP


#endif
