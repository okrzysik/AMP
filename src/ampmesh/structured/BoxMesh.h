#ifndef included_AMP_BoxMesh
#define included_AMP_BoxMesh

#include "ampmesh/Mesh.h"
#include "ampmesh/MeshID.h"
#include "ampmesh/MeshIterator.h"
#include "utils/AMP_MPI.h"

#ifdef USE_AMP_VECTORS
    namespace AMP {
    namespace LinearAlgebra {
        class Vector;
    }
    }
#endif

#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>
#include <vector>
#include <map>


namespace AMP {
namespace Mesh {

class structuredMeshElement;


/**
 * \class Mesh
 * \brief A class used to represent a logically rectangular box mesh
 * \details  This class provides routines for creating and managing a logically 
 *    rectangular mesh domain.  The mesh is described by the number of elements 
 *    in each direction and may be periodic along any given direction.  
 *    The database may specify some simple options to generate meshes:
 *    Generator - "cube", "sphere", "cyliner", "tube"
 *    Size - ndim array with the number of intervals in each direction.
 *           [nx,ny,nz] for box, [nr,nphi,nz] for cylinder, and [nr,nphi,ntheta] for sphere.
 *    Range - Array specifying the physical size of the mesh.
 *           cube: [ x-min  x-max  y-min  y-max  z-min  z-max ]
 *           sphere: [ r ]
 *           cyliner: [ r  z-min  z-max ]
 *           tube: [ r-min  r-max  z-min  z-max ]
 *    Periodic: Are any dimensions periodic (only applies to cubes)
 */
class BoxMesh: public AMP::Mesh::Mesh
{
public:


    /**
     * \brief Read in mesh files, partition domain, and prepare environment for simulation
     * \details  For trivial parallelsim, this method reads in the meshes on each processor.  Each
     * processor contains a piece of each mesh.  For massive parallelism, each mesh is on its own
     * communicator.  As such, some math libraries must be initialized accordingly.
     * \param params  Parameters for constructing a mesh from an input database
     */
    BoxMesh ( const MeshParameters::shared_ptr &params );


    /**
     * \brief   Estimate the number of elements in the mesh 
     * \details  This function will estimate the number of elements in the mesh. 
     *   This is used so that we can properly balance the meshes across multiple processors.
     *   Ideally this should be both an accurate estimate and very fast.  It should not require
     *   any communication and should not have to actually load a mesh.
     * \param params Parameters for constructing a mesh from an input database
     */
    static size_t estimateMeshSize( const MeshParameters::shared_ptr &params );


    //! Deconstructor
     ~BoxMesh ();


    /* Return the number of local element of the given type
     * \param type   Geometric type
     */
    virtual size_t  numLocalElements( const GeomType type ) const;


    /* Return the global number of elements of the given type
     * Note: depending on the mesh this routine may require global communication across the mesh.
     * \param type   Geometric type
     */
    virtual size_t  numGlobalElements( const GeomType type ) const;


    /* Return the number of ghost elements of the given type on the current processor
     * \param type   Geometric type
     */
    virtual size_t  numGhostElements( const GeomType type, const int gcw ) const;


    /**
     * \brief    Return an MeshIterator over the given geometric objects
     * \details  Return an MeshIterator over the given geometric objects
     * \param type   Geometric type to iterate over
     * \param gcw    Desired ghost cell width
     */
    virtual MeshIterator getIterator ( const GeomType type, const int gcw=0 ) const;


    /**
     * \brief    Return an MeshIterator over the given geometric objects on the surface
     * \details  Return an MeshIterator over the given geometric objects on the surface
     * \param type   Geometric type to iterate over
     * \param gcw    Desired ghost cell width
     */
    virtual MeshIterator getSurfaceIterator ( const GeomType type, const int gcw=0 ) const;


    /**
     * \brief    Return the list of all boundary ID sets in the mesh
     * \details  Return the list of all boundary ID sets in the mesh
     * Note: depending on the mesh this routine may require global communication across the mesh.
     */
    virtual std::vector<int> getBoundaryIDs ( ) const;


    /**
     * \brief    Return an MeshIterator over the given geometric objects on the given boundary ID set
     * \details  Return an MeshIterator over the given geometric objects on the given boundary ID set
     * \param type   Geometric type to iterate over
     * \param id     Boundary id for the elements (example: sideset id)
     * \param gcw    Desired ghost cell width
     */
    virtual MeshIterator getBoundaryIDIterator ( const GeomType type, const int id, const int gcw=0 ) const;

    /**
     * \brief    Return the list of all boundary ID sets in the mesh
     * \details  Return the list of all boundary ID sets in the mesh
     * Note: depending on the mesh this routine may require global communication across the mesh.
     */
    virtual std::vector<int> getBlockIDs ( ) const;


    /**
     * \brief    Return an MeshIterator over the given geometric objects on the given block ID set
     * \details  Return an MeshIterator over the given geometric objects on the given block ID set
     * \param type   Geometric type to iterate over
     * \param id     Block id for the elements (example: block id in cubit, subdomain in libmesh)
     * \param gcw    Desired ghost cell width
     */
    virtual MeshIterator getBlockIDIterator ( const GeomType type, const int id, const int gcw=0 ) const;


    /**
     * \brief    Return an MeshIterator constructed through a set operation of two other MeshIterators.
     * \details  Return an MeshIterator constructed through a set operation of two other MeshIterators.
     * \param OP Set operation to perform.
     *           Union - Perform a union of the iterators ( A U B )
     *           Intersection - Perform an intersection of the iterators ( A n B )
     *           Complement - Perform a compliment of the iterators ( A - B )
     * \param A  Pointer to MeshIterator A
     * \param B  Pointer to MeshIterator B
     */
    static MeshIterator getIterator ( SetOP OP, const MeshIterator &A, const MeshIterator &B);
 

    /**
     * \brief    Return a mesh element given it's id.
     * \details  This function queries the mesh to get an element given the mesh id.
     *    This function is only required to return an element if the id is local.
     *    Ideally, this should be done in O(1) time, but the implimentation is up to 
     *    the underlying mesh.  The base class provides a basic implimentation, but 
     *    uses mesh iterators and requires O(N) time on the number of elements in the mesh.
     * \param id    Mesh element id we are requesting.
     */
    virtual MeshElement getElement ( const MeshElementID &id ) const;


    /**
     * \brief    Displace the entire mesh
     * \details  This function will displace the entire mesh by a scalar value.
     *   This function is a blocking call for the mesh communicator, and requires
     *   the same value on all processors.  The displacement vector should be the 
     *   size of the physical dimension.
     * \param x  Displacement vector
     */
    virtual void displaceMesh( std::vector<double> x );


#ifdef USE_AMP_VECTORS
    /**
     * \brief    Displace the entire mesh
     * \details  This function will displace the entire mesh by displacing
     *   each node by the values provided in the vector.  This function is 
     *   a blocking call for the mesh communicator
     * \param x  Displacement vector.  Must have N DOFs per node where N 
     *           is the physical dimension of the mesh.
     */
    virtual void displaceMesh ( boost::shared_ptr<const AMP::LinearAlgebra::Vector> x );
#endif


protected:

    // Structure to uniquely identify an element
    struct MeshElementIndex{
        unsigned char type;     // Mesh element type
        unsigned char side;     // Are we dealing with x, y, or z faces/edges
        int index[3];           // Global x, y, z index (may be negitive with periodic boundaries)
        MeshElementIndex() {    // Empty constructor
            type = 0;
            side = 0;
            index[0] = 0;
            index[1] = 0;
            index[2] = 0;
        };
        MeshElementIndex(unsigned char type_in, unsigned char side_in, int x, int y=0, int z=0) {
            type = type_in;
            side = side_in;
            index[0] = x;
            index[1] = y;
            index[2] = z;
        };
        inline bool operator== (const MeshElementIndex& rhs ) const {
            return type==rhs.type && side==rhs.side && index[0]==rhs.index[0] 
                && index[1]==rhs.index[1] && index[2]==rhs.index[2];
        }
        inline bool operator!= (const MeshElementIndex& rhs ) const {
            return type!=rhs.type || side!=rhs.side || index[0]!=rhs.index[0] 
                || index[1]!=rhs.index[1] || index[2]!=rhs.index[2];
        }
        inline bool operator> (const MeshElementIndex& rhs ) const {
            if ( type < rhs.type ) { return false; }
            else if ( type > rhs.type ) { return true; }
            if ( side < rhs.side ) { return false; }
            else if ( side > rhs.side ) { return true; }
            for (int i=2; i>=0; i--) {
                if ( index[i] < rhs.index[i] ) { return false; }
                else if ( index[i] > rhs.index[i] ) { return true; }
            }
            return false;
        }
        inline bool operator>= (const MeshElementIndex& rhs ) const {
            return this->operator>(rhs) || this->operator==(rhs);
        }
        inline bool operator< (const MeshElementIndex& rhs ) const {
            return !this->operator>(rhs) && !this->operator==(rhs);
        }
        inline bool operator<= (const MeshElementIndex& rhs ) const {
            return !this->operator>(rhs);
        }
    };

    // Function to initialize the mesh data once the mesh has been created
    void initialize();

    // Helper function to return the indices of the local block owned by the given processor
    std::vector<int> getLocalBlock(unsigned int rank) const;

    // Helper function to return the block and owning rank of the given MeshElementIndex
    std::vector<int> getOwnerBlock(const MeshElementIndex index, unsigned int &rank) const;

    // Helper function to fill the node data for a uniform cartesian mesh
    static void fillCartesianNodes(int dim, const int* globalSize, const double *range, 
        const std::vector<MeshElementIndex> &index, std::vector<double> *coord);

    // Internal data
    bool d_isPeriodic[3];                   // Which directions are periodic
    int d_size[3];                          // The size of the logical domain in each direction
    std::vector<int> d_localSize[3];        // A vector indicating the local sizes for the sub-boxes 
    std::vector<MeshElementIndex> d_index;  // The indicies of the nodes we are storing
    std::vector<double> d_coord[3];         // The coordinates of the nodes

    // Basic mesh data
    typedef boost::shared_ptr<std::vector<MeshElementIndex> >   ElementIndexList;
    std::vector<ElementIndexList>   d_elements[4];
    size_t N_global[4];

    // Boundary and id set data
    std::vector<ElementIndexList>  d_surface_list[4];
    std::vector<int>               d_ids;
    std::map<std::pair<int,GeomType>,std::vector<ElementIndexList> >  d_id_list;

    // Friend functions to access protected functions    
    friend class structuredMeshElement;

};


} // Mesh namespace
} // AMP namespace

#endif

