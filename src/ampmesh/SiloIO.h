#ifndef included_AMP_SiloIO
#define included_AMP_SiloIO

#include <map>
#include <set>
#include <sstream>
#include <string.h>
#include <vector>

#include "AMP/ampmesh/Mesh.h"
#include "AMP/utils/Writer.h"

#ifdef USE_AMP_VECTORS
#include "AMP/vectors/Vector.h"
#endif
#ifdef USE_AMP_MATRICES
#include "AMP/matrices/Matrix.h"
#endif
#ifdef USE_EXT_SILO
#include <silo.h>
#endif


namespace AMP {
namespace Mesh {


/**
 * \class SiloIO
 * \brief A class used to abstract away reading/writing files for visualization
 * \details  This class provides routines for reading, accessing and writing meshes and vectors
 * using silo.
 */
class SiloIO : public AMP::Utilities::Writer
{
public:
    //!  Default constructor
    SiloIO();

    //!  Default destructor
    virtual ~SiloIO();

    //!  Function to return the file extension
    virtual std::string getExtension() override;

    //!  Function to read a file
    virtual void readFile( const std::string &fname ) override;

    //!  Function to write a file
    virtual void
    writeFile( const std::string &fname, size_t iteration_count, double time = 0 ) override;

    /**
     * \brief    Function to register a mesh
     * \details  This function will register a mesh with the silo writer.
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
                               std::string path = std::string() ) override;

#ifdef USE_AMP_VECTORS
    /**
     * \brief    Function to register a vector
     * \details  This function will register a vector with the silo writer.
     * \param vec   The vector we want to write
     * \param mesh  The mesh we want to write the vector over.
     *              Note: the vector must completely cover the mesh (silo limitiation).
     *              Note: mesh does not have to be previously registered with registerMesh.
     * \param type  The entity type we want to save (vertex, face, cell, etc.)
     *              Note: silo only supports writing one entity type.  If the vector
     *              spans multiple entity type (eg cell+vertex) the user should register
     *              the vector multiple times (one for each entity type).
     * \param name  Optional name for the vector.
     */
    virtual void registerVector( AMP::LinearAlgebra::Vector::shared_ptr vec,
                                 AMP::Mesh::Mesh::shared_ptr mesh,
                                 AMP::Mesh::GeomType type,
                                 const std::string &name = "" ) override;

    /**
     * \brief    Function to register a vector
     * \details  This function will register a vector with the writer.
     *     This version of registerVector only stores the raw data.  It is not associated with a
     * mesh.
     * \param vec   The vector we want to write
     */
    virtual void registerVector( AMP::LinearAlgebra::Vector::shared_ptr vec ) override;
#endif

#ifdef USE_AMP_MATRICES
    /**
     * \brief    Function to register a matrix
     * \details  This function will register a matrix with the writer.
     *     This version of registerMatrix only stores the raw data.  It is not associated with a
     * mesh.
     * \param mat   The matrix we want to write
     */
    virtual void registerMatrix( AMP::LinearAlgebra::Matrix::shared_ptr mat ) override;
#endif


private:
    // Structure used to hold data for the silo meshes
    struct siloBaseMeshData {
        AMP::Mesh::MeshID id;             // Unique ID to identify the mesh
        AMP::Mesh::Mesh::shared_ptr mesh; // Pointer to the mesh
        int rank;                         // Rank of the current processor on the mesh
        int ownerRank;                    // Global rank of the processor that "owns" the mesh
        std::string meshName;             // Name of the mesh in silo
        std::string path;                 // Path to the mesh in silo
        std::string file;                 // File that will contain the mesh
        std::vector<std::string> varName; // Names of variables for each mesh
        std::vector<AMP::Mesh::GeomType> varType; // Types of variables for each mesh
        std::vector<int> varSize;                 // Number of unknowns per point
#ifdef USE_AMP_VECTORS
        std::vector<AMP::LinearAlgebra::Vector::shared_ptr> vec; // Vectors for each mesh
#endif
        // Function to count the number of bytes needed to pack the data
        size_t size() const;
        // Function to pack the data to a byte array (note: some info may be lost)
        void pack( char * ) const;
        // Function to unpack the data from a byte array (note: some info may be lost)
        static siloBaseMeshData unpack( const char * );
    };

    // Structure used to hold data for the silo multimeshes
    struct siloMultiMeshData {
        AMP::Mesh::MeshID id;                 // Unique ID to identify the mesh
        AMP::Mesh::Mesh::shared_ptr mesh;     // Pointer to the mesh
        int ownerRank;                        // Global rank of the processor that "owns" the mesh
                                              // (usually rank 0 on the mesh comm)
        std::string name;                     // Name of the multimesh in silo
        std::vector<siloBaseMeshData> meshes; // Base mesh info needed to construct the mesh data
        std::vector<std::string> varName;     // Vectors for each mesh
        // Function to count the number of bytes needed to pack the data
        size_t size() const;
        // Function to pack the data to a byte array
        void pack( char * ) const;
        // Function to unpack the data from a byte array
        static siloMultiMeshData unpack( const char * );
        // Constructors
        siloMultiMeshData() : id(), ownerRank( -1 ) {}
        siloMultiMeshData( const siloMultiMeshData & );
        siloMultiMeshData &operator=( const siloMultiMeshData & );
        // Destructor
        ~siloMultiMeshData() {}
    };

    // Function to syncronize multimesh data
    void syncMultiMeshData( std::map<AMP::Mesh::MeshID, siloMultiMeshData> &data,
                            int root = -1 ) const;

    // Function to syncronize variable lists
    void syncVariableList( std::set<std::string> &data, int root = -1 ) const;

// Function to write a single mesh
#ifdef USE_EXT_SILO
    void writeMesh( DBfile *file, const siloBaseMeshData &data, int cycle, double time );
#endif

    // Function to determine which base mesh ids to register a vector with
    std::vector<AMP::Mesh::MeshID> getMeshIDs( AMP::Mesh::Mesh::shared_ptr mesh );

    // Function to write the summary file (the file should already be created, ready to reopen)
    // This function requires global communication
    void writeSummary( std::string filename, int cycle, double time );

    // The dimension
    int d_dim;

    // List of all meshes and thier ids
    std::map<AMP::Mesh::MeshID, siloBaseMeshData> d_baseMeshes;
    std::map<AMP::Mesh::MeshID, siloMultiMeshData> d_multiMeshes;

    // List of all variables
    std::set<std::string> d_varNames;

// List of all vectors that have been registered
#ifdef USE_AMP_VECTORS
    std::vector<AMP::LinearAlgebra::Vector::shared_ptr> d_vectors;
#endif
};
} // namespace Mesh
} // namespace AMP

#endif
