#ifndef included_AMP_BoxMesh
#define included_AMP_BoxMesh

#include "ampmesh/Mesh.h"
#include "ampmesh/MeshID.h"
#include "ampmesh/MeshIterator.h"

#ifdef USE_AMP_VECTORS
namespace AMP {
namespace LinearAlgebra {
class Vector;
}
}
#endif

#include "utils/shared_ptr.h"
#include <map>
#include <vector>
#include <array>


namespace AMP {
namespace Mesh {

class structuredMeshElement;
class structuredMeshIterator;


/**
 * \class BoxMesh
 * \brief A class used to represent a logically rectangular box mesh
 * \details  This class provides routines for creating and managing a logically
 *    rectangular mesh domain.  The mesh is described by the number of elements
 *    in each direction and may be periodic along any given direction.
 *    The database may specify some simple options to generate meshes:
\verbatim
   MeshName - Name of the mesh
   dim - Dimension of the mesh
   Generator - "cube", "circle, "cylinder", "tube", "sphere", "shell"
   Size - ndim array with the number of intervals in each direction.
          nx, ny, nz are the number of intervals in each direction, nr is the number of intervals
          in the radius, and nphi is the number of intervals in the asmuthal direction.
          cube (2d) - [ nx, ny      ]
          cube (3d) - [ nx, ny, nz  ]
          circle    - [ nr          ]
          cylinder  - [ nr, nz      ]
          tube      - [ nr, nphi, nz]
          sphere    - [ nr          ]
          shell     - [ nr, nphi    ]
   Range - Array specifying the physical size of the mesh.
          cube (2d) - [ x-min, x-max, y-min, y-max, z-min, z-max ]
          cube (3d) - [ x-min, x-max, y-min, y-max               ]
          circle    - [ r                                        ]
          cylinder  - [ r,     z_min, z_max                      ]
          tube      - [ r_min, r_max, z_min, z_max               ]
          sphere    - [ r                                        ]
          shell     - [ r_min, r_max                             ]
   Periodic - Are any dimensions periodic (optional)
          cube (2d) - [ x_dir, y_dir ]
          cube (3d) - [ x_dir, y_dir, z_dir ]
          circle    - Not supported
          cylinder  - [ z_dir ]
          tube      - [ z_dir ]
          sphere    - Not supported
          shell     - Not supported
   GCW -  The maximum ghost cell width to support (optional, default is 1)
   x_offset - Offset in x-direction (optional)
   y_offset - Offset in y-direction (optional)
   z_offset - Offset in z-direction (optional)
\endverbatim
 */
class BoxMesh : public AMP::Mesh::Mesh
{
public:
    /**
     * \class Box
     * \brief Structure to identify a logical box
     * \details  This class contains a logical box
     */
    class Box
    {
    public:
        /**
         * \brief   Default constructor
         * \details  The default constructor
         * \param ifirst  First x-coordinate
         * \param ilast   Last x-coordinate
         * \param jfirst  First y-coordinate
         * \param jlast   Last y-coordinate
         * \param kfirst  First z-coordinate
         * \param klast   Last z-coordinate
         */
        inline explicit Box(
            int ifirst, int ilast, int jfirst = 0, int jlast = 0, int kfirst = 0, int klast = 0 );
        inline Box(); //!< Empty constructor
        int first[3]; //!< Starting element
        int last[3];  //!< Ending element
    private:
    };

    /**
     * \class MeshElementIndex
     * \brief Structure to uniquely identify an element
     * \details  This class help convert between logical coordinates and the mesh element of
     * interest
     */
    class MeshElementIndex
    {
    public:
        //! Empty constructor
        inline MeshElementIndex();
        /**
         * \brief   Default constructor
         * \details  The default constructor
         * \param type  Element type
         * \param side  Side of the parent element (ignored if it is the parent or vertex)
         * \param x     Logical coordinate of the element
         * \param y     Logical coordinate of the element
         * \param x     Logical coordinate of the element
         */
        inline explicit MeshElementIndex(
            GeomType type, unsigned char side, int x, int y = 0, int z = 0 );
        inline void reset(
            GeomType type, unsigned char side, int x, int y = 0, int z = 0 );
        inline bool isNull() const { return d_side==255; }
        inline bool operator==( const MeshElementIndex &rhs ) const; //!< Operator ==
        inline bool operator!=( const MeshElementIndex &rhs ) const; //!< Operator !=
        inline bool operator>( const MeshElementIndex &rhs ) const;  //!< Operator >
        inline bool operator>=( const MeshElementIndex &rhs ) const; //!< Operator >=
        inline bool operator<( const MeshElementIndex &rhs ) const;  //!< Operator <
        inline bool operator<=( const MeshElementIndex &rhs ) const; //!< Operator <=
        inline int index( int d ) const { return d_index[d]; }
        inline int& index( int d ) { return d_index[d]; }
        inline GeomType type() const { return static_cast<GeomType>(d_type); }
        inline unsigned char side() const { return d_side; }
        static inline size_t numElements( const MeshElementIndex& first, const MeshElementIndex& last );
    private:
        unsigned char d_type; //!<  Mesh element type
        unsigned char d_side; //!<  Are we dealing with x, y, or z faces/edges
        int d_index[3];       //!<  Global x, y, z index (may be negitive with periodic boundaries)
        friend class structuredMeshElement;
    };

public:
    /**
     * \brief Read in mesh files, partition domain, and prepare environment for simulation
     * \details  For trivial parallelsim, this method reads in the meshes on each processor.  Each
     * processor contains a piece of each mesh.  For massive parallelism, each mesh is on its own
     * communicator.  As such, some math libraries must be initialized accordingly.
     * \param params  Parameters for constructing a mesh from an input database
     */
    static AMP::shared_ptr<BoxMesh> generate( MeshParameters::shared_ptr params );


    //! Virtual function to copy the mesh (allows use to proply copy the derived class)
    virtual AMP::shared_ptr<Mesh> copy() const override = 0;


    /**
     * \brief   Estimate the number of elements in the mesh
     * \details  This function will estimate the number of elements in the mesh.
     *   This is used so that we can properly balance the meshes across multiple processors.
     *   Ideally this should be both an accurate estimate and very fast.  It should not require
     *   any communication and should not have to actually load a mesh.
     * \param params Parameters for constructing a mesh from an input database
     */
    static size_t estimateMeshSize( const MeshParameters::shared_ptr &params );

    /**
     * \brief   Estimate the number of elements in the mesh
     * \details  This function will estimate the number of elements in the mesh.
     *   This is used so that we can properly balance the meshes across multiple processors.
     *   Ideally this should be both an accurate estimate and very fast.  It should not require
     *   any communication and should not have to actually load a mesh.
     * \param params Parameters for constructing a mesh from an input database
     */
    static std::vector<size_t> estimateLogicalMeshSize( const MeshParameters::shared_ptr &params );


    /**
     * \brief   Return the maximum number of processors that can be used with the mesh
     * \details  This function will return the maximum number of processors that can
     *   be used with the mesh.
     * \param params Parameters for constructing a mesh from an input database
     */
    static size_t maxProcs( const MeshParameters::shared_ptr &params );


    //! Deconstructor
    virtual ~BoxMesh();


    /* Return the number of local element of the given type
     * \param type   Geometric type
     */
    virtual size_t numLocalElements( const GeomType type ) const override final;


    /* Return the global number of elements of the given type
     * Note: depending on the mesh this routine may require global communication across the mesh.
     * \param type   Geometric type
     */
    virtual size_t numGlobalElements( const GeomType type ) const override final;


    /* Return the number of ghost elements of the given type on the current processor
     * \param type   Geometric type
     */
    virtual size_t numGhostElements( const GeomType type, const int gcw ) const override final;


    /**
     * \brief    Return an MeshIterator over the given geometric objects
     * \details  Return an MeshIterator over the given geometric objects
     * \param type   Geometric type to iterate over
     * \param gcw    Desired ghost cell width
     */
    virtual MeshIterator getIterator( const GeomType type, const int gcw = 0 ) const override final;


    /**
     * \brief    Return an MeshIterator over the given geometric objects on the surface
     * \details  Return an MeshIterator over the given geometric objects on the surface
     * \param type   Geometric type to iterate over
     * \param gcw    Desired ghost cell width
     */
    virtual MeshIterator getSurfaceIterator( const GeomType type, const int gcw = 0 ) const override final;


    /**
     * \brief    Return the list of all boundary ID sets in the mesh
     * \details  Return the list of all boundary ID sets in the mesh
     * Note: depending on the mesh this routine may require global communication across the mesh.
     */
    virtual std::vector<int> getBoundaryIDs() const override final;


    /**
     * \brief    Return an MeshIterator over the given geometric objects on the given boundary ID
     * set
     * \details  Return an MeshIterator over the given geometric objects on the given boundary ID
     * set
     * \param type   Geometric type to iterate over
     * \param id     Boundary id for the elements (example: sideset id)
     * \param gcw    Desired ghost cell width
     */
    virtual MeshIterator
    getBoundaryIDIterator( const GeomType type, const int id, const int gcw = 0 ) const override final;

    /**
     * \brief    Return the list of all boundary ID sets in the mesh
     * \details  Return the list of all boundary ID sets in the mesh
     * Note: depending on the mesh this routine may require global communication across the mesh.
     */
    virtual std::vector<int> getBlockIDs() const override final;


    /**
     * \brief    Return an MeshIterator over the given geometric objects on the given block ID set
     * \details  Return an MeshIterator over the given geometric objects on the given block ID set
     * \param type   Geometric type to iterate over
     * \param id     Block id for the elements (example: block id in cubit, subdomain in libmesh)
     * \param gcw    Desired ghost cell width
     */
    virtual MeshIterator
    getBlockIDIterator( const GeomType type, const int id, const int gcw = 0 ) const override final;


    /**
     * \brief    Return an MeshIterator constructed through a set operation of two other
     * MeshIterators.
     * \details  Return an MeshIterator constructed through a set operation of two other
     * MeshIterators.
     * \param OP Set operation to perform.
     *           Union - Perform a union of the iterators ( A U B )
     *           Intersection - Perform an intersection of the iterators ( A n B )
     *           Complement - Perform a compliment of the iterators ( A - B )
     * \param A  Pointer to MeshIterator A
     * \param B  Pointer to MeshIterator B
     */
    static MeshIterator getIterator( SetOP OP, const MeshIterator &A, const MeshIterator &B );


    /**
     * \brief    Return a mesh element given it's id.
     * \details  This function queries the mesh to get an element given the mesh id.
     *    This function is only required to return an element if the id is local.
     *    Ideally, this should be done in O(1) time, but the implimentation is up to
     *    the underlying mesh.  The base class provides a basic implimentation, but
     *    uses mesh iterators and requires O(N) time on the number of elements in the mesh.
     * \param id    Mesh element id we are requesting.
     */
    virtual MeshElement getElement( const MeshElementID &id ) const override final;


    /**
     * \brief    Return the parent elements of the given mesh element
     * \details  This function queries the mesh to get an element given the mesh id,
     *    then returns the parent elements that have the element as a child
     * \param elem  Mesh element of interest
     * \param type  Element type of the parents requested
     */
    virtual std::vector<MeshElement> getElementParents( const MeshElement &elem,
                                                        const GeomType type ) const override final;


    //! Return the global logical box
    inline Box getGlobalBox( int gcw = 0 ) const;


    //! Return the local logical box
    inline Box getLocalBox( int gcw = 0 ) const;


    //! Return the bool vector indicating which directions are periodic
    inline std::vector<bool> periodic() const;

    //! Return the size of the mesh
    inline std::vector<size_t> size() const;


    /**
     * \brief    Return a mesh element given it's id.
     * \details  This function queries the mesh to get an element given the mesh id.
     *    This function is only required to return an element if the id is local.
     *    Ideally, this should be done in O(1) time, but the implimentation is up to
     *    the underlying mesh.  The base class provides a basic implimentation, but
     *    uses mesh iterators and requires O(N) time on the number of elements in the mesh.
     * \param index    Mesh element index we are requesting.
     */
    MeshElement getElement( const MeshElementIndex &index ) const;

    /**
     * \brief    Return a mesh element's coordinates given it's id.
     * \details  This function queries the mesh to get an element's coordinates given the mesh id.
     *    Ideally, this should be done in O(1) time, but the implimentation is up to
     *    the underlying mesh.  
     * \param[in] index     Mesh element index we are requesting.
     * \param[out] pos      Mesh element coordinates
     */
    virtual void coord( const MeshElementIndex &index, double *pos ) const = 0;

    //! Check if the element is on the given boundary id
    virtual bool isOnBoundary( const MeshElementIndex &index, int id ) const;


public: // BoxMesh specific functionality

    /**
     * \brief    Return the logical coordinates
     * \details  This function queries the mesh to get the logical coordinates in [0,1]
     *     from the physical coordinates.  Not all meshes support this functionallity.
     * \param[in] x         Physical coordinates
     * @return              Returns the logical coordinates
     */
    virtual std::array<double,3> physicalToLogical( const double *x ) const = 0;

    /**
     * \brief    Return the element containing the point
     * \details  This function queries the mesh to get the element index given a logical point
     * \param[in] x         Logical coordinates in [0,1]
     * @return              Returns the element containing the logical point.
     *                      Note: it will return a null index (isNull) if no element of
     *                      the given type contains the point.
     */
    MeshElementIndex getElementFromLogical( const std::array<double,3>& x, GeomType type ) const;

    /**
     * \brief    Return the element containing the point
     * \details  This function queries the mesh to get the element index given a physical point.
     *    This functionallity requires physicalToLogical which may not be supported by all meshes.
     * \param[in] x         Physical coordinates
     * @return              Returns the element containing the physical point.
     *                      Note: it will return a null index (isNull) if no element of
     *                      the given type contains the point.
     */
    MeshElementIndex getElementFromPhysical( const double *x, GeomType type ) const;

    //! Convert the MeshElementIndex to the MeshElementID
    inline MeshElementID convert( const MeshElementIndex& id ) const;

    //! Convert the MeshElementID to the MeshElementIndex
    inline MeshElementIndex convert( const MeshElementID& id ) const;

protected:

    // Constructor
    BoxMesh( MeshParameters::shared_ptr );
    BoxMesh( const BoxMesh& );

    // Function to initialize the mesh data once the logical mesh info has been created
    void initialize();

    // Function to finalize the mesh data once the coordinates have been set
    void finalize();

    // Helper function to return the indices of the local block owned by the given processor
    inline std::array<int,6> getLocalBlock( unsigned int rank ) const;

    // Helper functions to identify the iterator blocks
    typedef std::vector<std::pair<MeshElementIndex,MeshElementIndex>> ElementBlocks;
    ElementBlocks getIteratorRange( std::array<int,6> range, const GeomType type, const int gcw ) const;
    static ElementBlocks intersect( const ElementBlocks& v1, const ElementBlocks& v2 );

    // Helper function to create an iterator from an ElementBlocks list
    inline MeshIterator createIterator( const ElementBlocks& list ) const;

    // Helper function to fill the node data for a uniform cartesian mesh
    static void fillCartesianNodes( int dim,
                                    const int *globalSize,
                                    const double *range,
                                    const std::vector<MeshElementIndex> &index,
                                    std::vector<double> *coord );

    // Internal data
    std::array<bool,3> d_isPeriodic;        // Which directions are periodic
    std::array<int,3>  d_globalSize;        // The size of the logical domain in each direction
    std::array<int,3>  d_numBlocks;         // The number of local box in each direction
    std::array<int,6>  d_surfaceId;         // For each surface which id is it part of (if any)
    std::array<bool,6> d_onSurface;         // For each surface which id is it part of (if any)
    ElementBlocks d_globalSurfaceList[6][4]; // List of logical surface elements for each surface/type

    // Friend functions to access protected functions
    friend class structuredMeshElement;
    friend class structuredMeshIterator;

private:
    BoxMesh(); // Private empty constructor

    
};


} // Mesh namespace
} // AMP namespace


#include "ampmesh/structured/BoxMesh.inline.h"


#endif
