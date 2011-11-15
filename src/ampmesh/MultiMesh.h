#ifndef included_AMP_MultiMesh
#define included_AMP_MultiMesh

#include "ampmesh/Mesh.h"


namespace AMP {
namespace Mesh {



/**
 * \class MultiMesh
 * \brief A concrete mesh class for a multi-mesh
 *
 * \details  This class provides routines for creating and accessing multimeshes.
 *   The concept of a multimesh is a mesh that is composed of multiple meshes.
 *   This takes care of the need for a mesh manager while allowing all operations
 *   on a given mesh to apply to multiple meshes.  Note: all meshes within a multimesh
 *   are stored in a flat array.  This applies when we have a multimesh that may contain
 *   other multimeshes that may (or may not) overlap.
 */
class MultiMesh: public Mesh
{
public:

    /**
     * \param params Parameters for constructing a mesh from an input database
     * \brief Read in mesh files, partition domain, and prepare environment for simulation
     * \details  For trivial parallelsim, this method reads in the meshes on each processor.  Each
     * processor contains a piece of each mesh.  For massive parallelism, each mesh is on its own
     * communicator.  As such, some math libraries must be initialized accordingly.
     */
    MultiMesh ( const MeshParameters::shared_ptr &params );

    //! Deconstructor
     ~MultiMesh ();

    //! Function to copy the mesh (allows use to proply copy the derived class)
    Mesh copy() const;

    /**
     * \brief   Estimate the number of elements in the mesh 
     * \details  This function will estimate the number of elements in the mesh. 
     *   This is used so that we can properly balance the meshes across multiple processors.
     *   Ideally this should be both an accurate estimate and very fast.  It should not require
     *   any communication and should not have to actually load a mesh.
     * \param params Parameters for constructing a mesh from an input database
     */
    static size_t estimateMeshSize( const MeshParameters::shared_ptr &params );

    /* Return the number of local element of the given type
     * \param type   Geometric type
     */
    virtual size_t  numLocalElements( const GeomType type ) const;


    /* Return the global number of elements of the given type
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
     * \param gcw    Desired ghost 
     */
    virtual MeshIterator getIterator ( const GeomType type, const int gcw=0 );


    /**
     * \brief    Return the list of all ID sets in the mesh
     * \details  Return the list of all ID sets in the mesh
     */
    virtual std::vector<int> getIDSets ( );


    /**
     * \brief    Return an MeshIterator over the given geometric objects on the given ID set
     * \details  Return an MeshIterator over the given geometric objects on the given ID set
     * \param type   Geometric type to iterate over
     * \param type   id for the elements (example: nodeset id)
     */
    virtual MeshIterator getIDsetIterator ( const GeomType type, const int id, const int gcw=0 );




private:

    //! Empty constructor for a mesh
    MultiMesh ( ) {};

    //! Function to create the databases for the meshes within the multimesh
    static std::vector<boost::shared_ptr<AMP::Database> >  createDatabases(boost::shared_ptr<AMP::Database> database);

    //! A list of all meshes in the multimesh
    std::vector<AMP::Mesh::Mesh::shared_ptr> d_meshes;

};

} // Mesh namespace
} // AMP namespace

#endif

