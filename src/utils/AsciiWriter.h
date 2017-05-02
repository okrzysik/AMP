#ifndef included_AMP_AsciiWriter
#define included_AMP_AsciiWriter

#include <map>
#include <set>
#include <sstream>
#include <string.h>
#include <vector>

#include "utils/Writer.h"

#ifdef USE_AMP_MESH
#include "ampmesh/Mesh.h"
#endif
#ifdef USE_AMP_VECTORS
#include "vectors/Vector.h"
#endif
#ifdef USE_AMP_MATRICES
#include "matrices/Matrix.h"
#endif


namespace AMP {
namespace Utilities {


/**
 * \class AsciiWriter
 * \brief A class used to abstract away reading/writing files.
 * \details  This class provides routines for reading, accessing and writing vectors
 *    and matrices using a simple ASCII format.
 *    Note: this format only supports a decomposition of 1, and will write all data using rank0
 */
class AsciiWriter : public AMP::Utilities::Writer
{
public:
    //!  Default constructor
    AsciiWriter();

    //!  Default destructor
    virtual ~AsciiWriter();

    //!  Function to return the file extension
    virtual std::string getExtension();

    //!  Function to read a file
    virtual void readFile( const std::string &fname );

    //!  Function to write a file
    virtual void writeFile( const std::string &fname, size_t iteration_count, double time=0 );

#ifdef USE_AMP_MESH
    /**
     * \brief    Function to register a mesh
     * \details  This function will register a mesh with the writer.
     *           Note: if mesh is a MultiMesh, it will register all sub meshes.
     * \param mesh  The mesh to register
     * \param level How many sub meshes do we want?
     *              0: Only register the local base meshes (advanced users only)
     *              1: Register current mesh only (default)
     *              2: Register all meshes (do not seperate for the ranks)
     *              3: Register all mesh pieces including the individual ranks

     * \param path  The directory path for the mesh.  Default is an empty string.
     */
    virtual void registerMesh( AMP::Mesh::Mesh::shared_ptr mesh,
                               int level        = 1,
                               std::string path = std::string() );
#endif

#if defined( USE_AMP_VECTORS ) && defined( USE_AMP_MESH )
    /**
     * \brief    Function to register a vector
     * \details  This function will register a vector with the writer and register it with the given
     * mesh.
     *     This version of registerVector allows the data to be "stored" on the mesh for
     * visualization
     *     or mesh-based operations.
     * \param vec   The vector we want to write
     * \param mesh  The mesh we want to write the vector over.
     *              Note: for many writers the vector must completely cover the mesh.
     *              Note: mesh does not have to be previously registered with registerMesh.
     * \param type  The entity type we want to save (vertex, face, cell, etc.)
     *              Note: some writers only supports writing one entity type.  If the vector
     *              spans multiple entity type (eg cell+vertex) the user should register
     *              the vector multiple times (one for each entity type).
     * \param name  Optional name for the vector.
     */
    virtual void registerVector( AMP::LinearAlgebra::Vector::shared_ptr vec,
                                 AMP::Mesh::Mesh::shared_ptr mesh,
                                 AMP::Mesh::GeomType type,
                                 const std::string &name = "" );
#endif

#if defined( USE_AMP_VECTORS )
    /**
     * \brief    Function to register a vector
     * \details  This function will register a vector with the writer.
     *     This version of registerVector only stores the raw data.  It is not associated with a
     * mesh.
     * \param vec   The vector we want to write
     */
    virtual void registerVector( AMP::LinearAlgebra::Vector::shared_ptr vec );
#endif

#ifdef USE_AMP_MATRICES
    /**
     * \brief    Function to register a matrix
     * \details  This function will register a matrix with the writer.
     *     This version of registerMatrix only stores the raw data.  It is not associated with a
     * mesh.
     * \param mat   The matrix we want to write
     */
    virtual void registerMatrix( AMP::LinearAlgebra::Matrix::shared_ptr mat );
#endif


    // A typedef for ids
    typedef std::pair<unsigned int, unsigned int> global_id;

private:
// List of all vectors that have been registered
#ifdef USE_AMP_VECTORS
    std::map<global_id, AMP::LinearAlgebra::Vector::shared_ptr> d_vectors;
#endif

// List of all matrices that have been registered
#ifdef USE_AMP_MATRICES
    std::map<global_id, AMP::LinearAlgebra::Matrix::shared_ptr> d_matrices;
#endif

    // Helper functions
    static global_id getID( AMP_MPI local_comm, AMP_MPI global_comm );
#ifdef USE_AMP_VECTORS
    static AMP::LinearAlgebra::Vector::const_shared_ptr sendVecToRoot(
        AMP::LinearAlgebra::Vector::const_shared_ptr src_vec, int vec_root, AMP_MPI d_comm );
#endif
#ifdef USE_AMP_MATRICES
    static void sendRowToRoot( AMP::LinearAlgebra::Matrix::const_shared_ptr mat,
                               AMP_MPI d_comm,
                               int row,
                               std::vector<size_t> &col,
                               std::vector<double> &data );
#endif
};
}
}

#endif
