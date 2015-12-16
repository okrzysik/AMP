#ifndef included_AMP_MeshID
#define included_AMP_MeshID

#include "utils/shared_ptr.h"
#include <vector>


namespace AMP {
namespace Mesh {


//! Enumeration for basic mesh-based quantities
enum GeomType { Vertex = 0, Edge = 1, Face = 2, Volume = 3, null = 0xFF };


//! Typedef for a unsigned 64-bit integer
typedef unsigned long long int uint64;


struct MeshElementID;


/**
 * \struct MeshID
 * \brief  A structure used to identify the mesh
 * \details  This structure provides a unique id that can be used to identify a mesh.
 */
struct MeshID {

public:
    // Constructors used to initialize key values
    MeshID()
    {
        data = 0xFFFFFFFFFFFFFFFF; // set the mesh id to -1
    }
    MeshID( unsigned int root, unsigned int local_id )
    {
        // Use the first 32-bits for the rank of the root processor and the second 32-bits for the
        // local id
        data = root;
        data = ( data << 32 ) + local_id;
    }
    MeshID( uint64 id ) { data = id; }
    uint64 getData() const { return data; }
    // Overload key operators
    inline bool operator==( const MeshID &rhs ) const { return data == rhs.data; }
    inline bool operator!=( const MeshID &rhs ) const { return data != rhs.data; }
    inline bool operator>=( const MeshID &rhs ) const { return data >= rhs.data; }
    inline bool operator<=( const MeshID &rhs ) const { return data <= rhs.data; }
    inline bool operator>( const MeshID &rhs ) const { return data > rhs.data; }
    inline bool operator<( const MeshID &rhs ) const { return data < rhs.data; }
private:
    // We will store the data as a 64-bit data type
    // We do not want the user to get direct access to the data
    uint64 data;
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
    MeshElementID()
    {
        data[0] = 0xFFFFFFFFFFFFFFFF; // set the mesh id to -1
        data[1] = 0x000000FFFFFFFFFF; // set is_local to false, type_id to 0xFF, owner_rank to 0,
                                      // and local_id to -1
    }
    MeshElementID( bool isLocal,
                   GeomType type_id,
                   unsigned int local_ID,
                   unsigned int owner_rank_id,
                   MeshID mesh_ID )
    {
        // Copy the meshID
        data[0] = mesh_ID.data;
        // memcpy(data,(uint64*)&mesh_ID,sizeof(uint64));
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
        data[1] = ( ( (uint64) tmp ) << 32 ) + ( (uint64) local_ID );
    }
    // Overload key operators
    inline bool operator==( const MeshElementID &rhs ) const
    {
        uint64 test = ( data[1] ^ rhs.data[1] ) << 1;
        return data[0] == rhs.data[0] && test == 0;
    }
    inline bool operator!=( const MeshElementID &rhs ) const
    {
        uint64 test = ( data[1] ^ rhs.data[1] ) << 1;
        return data[0] != rhs.data[0] || test != 0;
    }
    inline bool operator>=( const MeshElementID &rhs ) const
    {
        // Sort by meshID first
        if ( data[0] < rhs.data[0] )
            return false;
        else if ( data[0] > rhs.data[0] )
            return true;
        // Sort by the remaining data
        uint64 d1 = data[1] << 1;
        uint64 d2 = rhs.data[1] << 1;
        return d1 >= d2;
    }
    inline bool operator>( const MeshElementID &rhs ) const
    {
        // Sort by meshID first
        if ( data[0] < rhs.data[0] )
            return false;
        else if ( data[0] > rhs.data[0] )
            return true;
        // Sort by the remaining data
        uint64 d1 = data[1] << 1;
        uint64 d2 = rhs.data[1] << 1;
        return d1 > d2;
    }
    inline bool operator<( const MeshElementID &rhs ) const
    {
        // Sort by meshID first
        if ( data[0] > rhs.data[0] )
            return false;
        else if ( data[0] < rhs.data[0] )
            return true;
        // Sort by the remaining data
        uint64 d1 = data[1] << 1;
        uint64 d2 = rhs.data[1] << 1;
        return d1 < d2;
    }
    inline bool operator<=( const MeshElementID &rhs ) const
    {
        // Sort by meshID first
        if ( data[0] > rhs.data[0] )
            return false;
        else if ( data[0] < rhs.data[0] )
            return true;
        // Sort by the remaining data
        uint64 d1 = data[1] << 1;
        uint64 d2 = rhs.data[1] << 1;
        return d1 <= d2;
    }
    // Access local data
    inline bool is_local() const { return ( data[1] >> 63 ) != 0; }
    inline GeomType type() const
    {
        unsigned char tmp = (unsigned char) ( ( data[1] >> 32 ) & 0x00FF );
        return (GeomType) tmp;
    }
    inline unsigned int local_id() const
    {
        unsigned int tmp = (unsigned int) ( data[1] & 0x00000000FFFFFFFF );
        return tmp;
    }
    inline unsigned int owner_rank() const
    {
        unsigned int tmp = data[1] >> 32;
        tmp              = ( tmp >> 8 ) & 0x007FFFFF;
        return tmp;
    }
    inline MeshID meshID() const { return MeshID( data[0] ); }
    inline void set_is_local( bool isLocal )
    {
        if ( isLocal )
            data[1] |= 0x8000000000000000;
        else
            data[1] &= 0x7FFFFFFFFFFFFFFF;
    }

private:
    // We will store the data as a 128-bit data type
    // The first 64 bits refer to the meshID
    // The next  1 bits refer to if the element is local
    // The next 23 bits refer to the processor id
    // The next  8 bits refer to the element type
    // The next 32 bits refer to the local id
    uint64 data[2];
};
}
}

#endif
