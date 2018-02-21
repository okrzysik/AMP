#ifndef included_AMP_Writer
#define included_AMP_Writer

#include <map>
#include <set>
#include <sstream>
#include <string.h>
#include <vector>

#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/shared_ptr.h"

#ifdef USE_AMP_MESH
#include "AMP/ampmesh/Mesh.h"
#endif
#ifdef USE_AMP_VECTORS
#include "AMP/vectors/Vector.h"
#endif
#ifdef USE_AMP_MATRICES
#include "AMP/matrices/Matrix.h"
#endif


namespace AMP {
namespace Utilities {


/**
 * \class Writer
 * \brief A class used to abstract away reading/writing files.
 * \details  This class provides routines for reading, accessing and writing meshes and vectors.
 *    The writers can be used to generate files for visualization or interfacing with other codes.
 */
class Writer
{
public:
    //!  Convenience typedef
    typedef AMP::shared_ptr<AMP::Utilities::Writer> shared_ptr;

    /**
     * \brief   Function to build a writer
     * \details This function will build a default writer for use.
     * \param type   Writer type:
     *               "None"  - An empty writer will be created
     *               "Silo"  - A silo writer will be created if silo is configured,
     *                        otherwise an empty writer will be created.
     *               "Ascii" - A simple ascii writer
     */
    static AMP::shared_ptr<AMP::Utilities::Writer> buildWriter( const std::string &type );

    /**
     * \brief   Function to build a writer
     * \details This function will build a default writer for use.
     * \param db   Input database for the writer
     */
    static AMP::shared_ptr<AMP::Utilities::Writer> buildWriter( AMP::shared_ptr<AMP::Database> db );

    //!  Default constructor
    Writer();

    //!  Default destructor
    virtual ~Writer();

    //!  Function to return the file extension
    virtual std::string getExtension() = 0;

    /**
     * \brief   Function to set the file decomposition
     * \details This function will set the method used for file IO.  When writing files,
     *    there are different decompositions that affect the performance and usability
     *    of the output files.  By default, this writer will generate a single file.
     * \param decomposition   Decomposition method to use:
     *             1:  This will write all of the data to a single file.
     *                 Note that this requires a serial write and will have the worst performance
     *             2:  Each processor will write a separate file and a separate
     *                 summary file will be written.  Note that this will have better performance
     *                 at large scale, but will write many files simultaneously.
     */
    virtual void setDecomposition( int decomposition );

    //!  Function to read a file
    virtual void readFile( const std::string &fname ) = 0;

    //!  Function to write a file
    virtual void writeFile( const std::string &fname, size_t iteration_count, double time = 0 ) = 0;

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
                               std::string path = std::string() ) = 0;
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
                                 const std::string &name = "" ) = 0;
#endif

#if defined( USE_AMP_VECTORS )
    /**
     * \brief    Function to register a vector
     * \details  This function will register a vector with the writer.
     *     This version of registerVector only stores the raw data.  It is not associated with a
     * mesh.
     * \param vec   The vector we want to write
     */
    virtual void registerVector( AMP::LinearAlgebra::Vector::shared_ptr vec ) = 0;
#endif

#ifdef USE_AMP_MATRICES
    /**
     * \brief    Function to register a matrix
     * \details  This function will register a matrix with the writer.
     *     This version of registerMatrix only stores the raw data.  It is not associated with a
     * mesh.
     * \param mat   The matrix we want to write
     */
    virtual void registerMatrix( AMP::LinearAlgebra::Matrix::shared_ptr mat ) = 0;
#endif


protected:  // Protected member functions

    // Given a filename, strip the directory information and create the directories if needed
    void createDirectories( const std::string& filename );


protected:  // Internal data

    // The comm of the writer
    AMP_MPI d_comm;

    // The decomposition to use
    int d_decomposition;
};
} // namespace Utilities
} // namespace AMP

#endif
