#ifndef included_AMP_Writer
#define included_AMP_Writer

#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"

#include <memory>
#include <string>
#include <vector>


// Declare some classes
namespace AMP::Mesh {
class Mesh;
enum class GeomType : uint8_t;
} // namespace AMP::Mesh
namespace AMP::LinearAlgebra {
class Vector;
class Matrix;
} // namespace AMP::LinearAlgebra


namespace AMP::Utilities {


/**
 * \class Writer
 * \brief A class used to abstract away reading/writing files.
 * \details  This class provides routines for reading, accessing and writing meshes and vectors.
 *    The writers can be used to generate files for visualization or interfacing with other codes.
 */
class Writer
{
public:
    struct WriterProperties {
        std::string type;            // Writer type: Silo, HDF5, Ascii
        std::string extension;       // The primary file extension for the writer
        bool registerMesh;           // Does the writer support registering a mesh
        bool registerVector;         // Does the writer support registering a vector
        bool registerVectorWithMesh; // Does the writer support registering a vector with a mesh
        bool registerMatrix;         // Does the writer support registering a matrix
        WriterProperties();
    };

public:
    /**
     * \brief   Function to build a writer
     * \details This function will build a default writer for use.
     * \param type   Writer type:
     *               "None"  - An empty writer will be created
     *               "Silo"  - A silo writer will be created if silo is configured,
     *                        otherwise an empty writer will be created.
     *               "Ascii" - A simple ascii writer
     *               "HDF5"  - A simple HDF5 writer
     * \param comm   Communicator to use
     */
    static std::shared_ptr<AMP::Utilities::Writer> buildWriter( std::string type,
                                                                AMP_MPI comm = AMP_COMM_WORLD );

    /**
     * \brief   Function to build a writer
     * \details This function will build a default writer for use.
     * \param db   Input database for the writer
     */
    static std::shared_ptr<AMP::Utilities::Writer> buildWriter( std::shared_ptr<AMP::Database> db );


public:
    //!  Default constructor
    Writer();

    //!  Default destructor
    virtual ~Writer();

    //! Delete copy constructor
    Writer( const Writer & ) = delete;

    //! Function to get the writer properties
    virtual WriterProperties getProperties() const = 0;

    //!  Function to return the file extension
    std::string getExtension() const;

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
    virtual void writeFile( const std::string &fname, size_t iteration, double time = 0 ) = 0;

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
    virtual void registerMesh( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                               int level               = 1,
                               const std::string &path = std::string() );

    /**
     * \brief    Function to register a vector
     * \details  This function will register a vector with the writer and register it with the given
     * mesh.
     *     This version of registerVector allows the data to be "stored" on the mesh for
     * visualization
     *     or mesh-based operations.
     * \param vec   The vector we want to write
     * \param mesh  The mesh we want to write the vector over.
     *              Note: any writers require the vector to completely cover the mesh.
     *              Note: mesh does not have to be previously registered with registerMesh.
     * \param type  The entity type we want to save (vertex, face, cell, etc.)
     *              Note: some writers only supports writing one entity type.
     *              If the vector spans multiple entity type (eg cell+vertex)  the user should
     *              register the vector multiple times (one for each entity type).
     * \param name  Optional name for the vector.
     */
    virtual void registerVector( std::shared_ptr<AMP::LinearAlgebra::Vector> vec,
                                 std::shared_ptr<AMP::Mesh::Mesh> mesh,
                                 AMP::Mesh::GeomType type,
                                 const std::string &name = "" );

    /**
     * \brief    Function to register a vector
     * \details  This function will register a vector with the writer.
     *     This version of registerVector only stores the raw data.
     *     It is not associated with a mesh.
     * \param vec   The vector we want to write
     * \param name  Optional name for the vector.
     */
    virtual void registerVector( std::shared_ptr<AMP::LinearAlgebra::Vector> vec,
                                 const std::string &name = "" );

    /**
     * \brief    Function to register a matrix
     * \details  This function will register a matrix with the writer.
     *     This version of registerMatrix only stores the raw data..
     *     It is not associated with a mesh.
     * \param mat   The matrix we want to write
     * \param name  Optional name for the vector.
     */
    virtual void registerMatrix( std::shared_ptr<AMP::LinearAlgebra::Matrix> mat,
                                 const std::string &name = "" );


protected: // Protected member functions
    // Given a filename, strip the directory information and create the directories if needed
    void createDirectories( const std::string &filename );


protected: // Internal data
    // The comm of the writer
    AMP_MPI d_comm;

    // The decomposition to use
    int d_decomposition;
};


} // namespace AMP::Utilities

#endif
