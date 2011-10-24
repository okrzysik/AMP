#ifndef included_AMP_MultiMesh
#define included_AMP_MultiMesh


namespace AMP { 
namespace Mesh {


/**
 * \class MultiMesh
 * \brief A class used to abstract away mesh management from an application
 *
 * \details  This class provides routines for reading, accessing and writing meshes
 * from an input database.  Also, it handles initialization of the chosen
 * mesh database library.  The database format is of the following:
 *\code
    NumberOfMeshes = n
    Mesh_1 {
      Filename = "filename"
      MeshName = "name of mesh used in rest of file"
      MeshGroup = "name of mesh group"
      x_offset = x
      y_offset = y
      z_offset = z
      NumberOfElements = number of elements
    }
      .
      .
      .
    Mesh_n {
      Filename = "filename"
      MeshName = "name of mesh used in rest of file"
      MeshGroup = "name of mesh group"
      x_offset = x
      y_offset = y
      z_offset = z
      NumberOfElements = number of elements
    }
    \endcode
 */

class MultiMesh: public Mesh
{
public:

    /**
     * \fn Mesh ( const MeshParameters::shared_ptr &params )
     * \param params Parameters for constructing a mesh from an input database
     * \brief Read in mesh files, partition domain, and prepare environment for simulation
     * \details  For trivial parallelsim, this method reads in the meshes on each processor.  Each
     * processor contains a piece of each mesh.  For massive parallelism, each mesh is on its own
     * communicator.  As such, some math libraries must be initialized accordingly.
     */
    MultiMesh ( const MeshParameters::shared_ptr &params );

    /**
     * \brief    Return an MultiMeshIterator over the meshes
     * \details  Return an MultiMeshIterator over the meshes
     */
    MultiMeshIterator getIterator ( );

protected:

    //!  Empty constructor for a mesh
    MultiMesh ( );

    //! The mesh parameters
    const MeshParameters::shared_ptr &params;

};



/**
 * \class MultiMeshIterator
 * \brief A iterator class used to iterate over meshes in a multimesh
 * \details  This class provides routines for iterating over meshes in a multimesh.
 */
class MultiMeshIterator: public iterator
{
public:

    /**
     * \brief MultiMeshIterator constructor
     * \details  This constructor will construct a new MultiMeshIterator from an existing mesh iterator.
     * \param mit Existing MeshIterator
     */
    MultiMeshIterator(const MultiMeshIterator& mit);

    //! Increment
    virtual MultiMeshIterator& operator++();
    
    //! Increment
    virtual MultiMeshIterator operator++(int);

    //! Check if two iterators are equal
    virtual bool operator==(const MultiMeshIterator& rhs);

    //! Check if two iterators are not equal
    virtual bool operator!=(const MultiMeshIterator& rhs);
    
    //! Dereference the iterator
    virtual Mesh& operator*();

protected:

    //!  Empty constructor for a MultiMeshIterator
    MultiMeshIterator ( );

};


}
}

#endif

