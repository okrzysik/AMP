#ifndef included_AMP_MeshElementID
#define included_AMP_MeshElementID

#include <vector>
#include <boost/shared_ptr.hpp>


namespace AMP {
namespace Mesh {


//! Enumeration for basic mesh-based quantities
enum GeomType { Vertex=0, Edge=1, Face=2, Volume=3, null=-1 };


typedef unsigned long long int    uint64;


/**
 * \struct MeshElementID
 * \brief  A structure used to identify the mesh element
 * \details  This structure provides a unique id that can be used to identify a mesh element.
 */
struct MeshElementID{

public:
    // Constructors used to initialize key values
	MeshElementID() {
        data[0] = 0xFFFFFFFFFFFFFFFF;       // set the mesh id to -1
        data[1] = 0x00000000FFFFFFFF;       // set is_local to false, local id to 0, owner_rank to 0, and local_id to -1
    }
	MeshElementID(bool isLocal, GeomType type_id, unsigned int local_ID, unsigned int owner_rank_id, size_t mesh_ID) {
        data[0] = mesh_ID;
        // Set the bit for is_local
        unsigned int tmp = 0x00000000;
        if ( isLocal )
            tmp = 0x80000000;
        // Add the owner_rank
        tmp += (0x007FFFFF&owner_rank_id) << 8;
        // Add the type_id
        char type = (char) type_id;
        tmp += ((unsigned char) type);
        // Combine the above data with the local_ID
        data[1] = (((uint64)tmp)<<32) + ((uint64)local_ID);
    }
    // Overload key operators
    inline bool operator== (const MeshElementID& rhs ) const {
        uint64 test = (data[1]^rhs.data[1])<<1;
        return data[0]==rhs.data[0] && test==0;
    }
    inline bool operator!= (const MeshElementID& rhs ) const {
        uint64 test = (data[1]^rhs.data[1])<<1;
        return data[0]!=rhs.data[0] || test!=0;
    }
    inline bool operator>= (const MeshElementID& rhs ) const {
        // Sort by meshID first
        if ( data[0] < rhs.data[0] )
            return false;
        else if ( data[0] > rhs.data[0] )
            return true;
        // Sort by the remaining data
        uint64 d1 = data[1]<<1;
        uint64 d2 = rhs.data[1]<<1;
        return d1>=d2;
    }
    inline bool operator> (const MeshElementID& rhs ) const {
        // Sort by meshID first
        if ( data[0] < rhs.data[0] )
            return false;
        else if ( data[0] > rhs.data[0] )
            return true;
        // Sort by the remaining data
        uint64 d1 = data[1]<<1;
        uint64 d2 = rhs.data[1]<<1;
        return d1>d2;
    }
    inline bool operator< (const MeshElementID& rhs ) const {
        // Sort by meshID first
        if ( data[0] > rhs.data[0] )
            return false;
        else if ( data[0] < rhs.data[0] )
            return true;
        // Sort by the remaining data
        uint64 d1 = data[1]<<1;
        uint64 d2 = rhs.data[1]<<1;
        return d1<d2;
    }
    inline bool operator<= (const MeshElementID& rhs ) const {
        // Sort by meshID first
        if ( data[0] > rhs.data[0] )
            return false;
        else if ( data[0] < rhs.data[0] )
            return true;
        // Sort by the remaining data
        uint64 d1 = data[1]<<1;
        uint64 d2 = rhs.data[1]<<1;
        return d1<=d2;
    }
    // Access local data
    inline bool is_local() { 
        return data[1]>>63;
    }
    inline GeomType type() { 
        char tmp = (unsigned char) ((data[1]>>32)&0x00FF);
        return (GeomType) tmp; 
    }
    inline unsigned int local_id() { 
        int tmp = (int) (data[1]&0x00000000FFFFFFFF);
        return tmp;
    }
    inline unsigned int owner_rank() { 
        unsigned int tmp = (unsigned int) ((data[1]>>40)&0x8FFFFFFF);
        return tmp; 
    }
    inline size_t meshID() { 
        return (size_t) data[0];
    }
    inline void set_is_local(bool isLocal) { 
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

