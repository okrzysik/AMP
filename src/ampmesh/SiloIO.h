#ifndef included_AMP_SiloIO
#define included_AMP_SiloIO

#include <string.h>
#include <sstream>
#include <map>

#ifdef USE_SILO
    #include <silo.h>
#endif

#ifdef USE_AMP_VECTORS
    #include "vectors/Vector.h"
#endif


namespace AMP { 
namespace Mesh {


/**
 * \class SiloIO
 * \brief A class used to abstract away reading/writing files for visualization
 * \details  This class provides routines for reading, accessing and writing meshes and vectors
 * using silo.
 */
class SiloIO 
{
public:

    //!  Convenience typedef
    typedef boost::shared_ptr<AMP::Mesh::SiloIO>  shared_ptr;

    //!  Default constructor
    SiloIO();

    //!  Function to return the file extension
    std::string getExtension() { return "silo"; }

    //!  Function to read a file
    void  readFile( const std::string &fname );

    //!  Function to write a file
    void  writeFile( const std::string &fname, size_t iteration_count );

    //!  Function to register a mesh (Note: if mesh is a MultiMesh, it will register all sub meshes)
    void registerMesh( AMP::Mesh::Mesh::shared_ptr mesh );

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
    //!  Function to register a vector
    void registerVector( AMP::LinearAlgebra::Vector::shared_ptr vec, AMP::Mesh::Mesh::shared_ptr mesh,
        AMP::Mesh::GeomType type, const std::string &name = "" );
#endif

private:
    
    // Structure used to hold data for the silo meshes
    struct siloBaseMeshData {
        AMP::Mesh::MeshID               id;         // Unique ID to identify the mesh
        AMP::Mesh::Mesh::shared_ptr     mesh;       // Pointer to the mesh
        int                             rank;       // Rank of the current processor on the mesh (used for name mangling)
        std::string                     meshName;   // Name of the mesh in silo
        std::string                     path;       // Path to the mesh in silo
        std::string                     file;       // File that will contain the mesh
        std::vector<std::string>        varName;    // List of the names of variables associated with each mesh
        std::vector<AMP::Mesh::GeomType> varType;   // List of the types of variables associated with each mesh
        std::vector<int> varSize;                   // Number of unknowns per point
        std::vector<AMP::LinearAlgebra::Vector::shared_ptr> vec; // List of the vectors associated with each mesh
        // Function to count the number of bytes needed to pack the data (note: some info may be lost)
        size_t size();
        // Function to pack the data to a byte array (note: some info may be lost)
        void pack( char* );
        // Function to unpack the data from a byte array (note: some info may be lost)
        static siloBaseMeshData unpack( char* );
    };

    // Structure used to hold data for the silo multimeshes
    struct siloMultiMeshData {
        AMP::Mesh::MeshID               id;         // Unique ID to identify the mesh
        std::string                     name;       // Name of the multimesh in silo
        std::vector<siloBaseMeshData>   meshes;     // Base mesh info needed to construct the mesh data
        // Function to count the number of bytes needed to pack the data
        size_t size();
        // Function to pack the data to a byte array
        void pack( char* );
        // Function to unpack the data from a byte array
        static siloMultiMeshData unpack( char* );
    };
    
    // Function to syncronize multimesh data across all processors
    void syncMultiMeshData( std::map<AMP::Mesh::MeshID,siloMultiMeshData> &data ) const;

    // Function to syncronize variable lists across all processors
    void syncVariableList( std::set<std::string> &data ) const;

    // Function to write a single mesh
    void writeMesh( DBfile *file, const siloBaseMeshData &data );

    // Function to write the summary file (the file should already be created, ready to reopen)
    // This function requires global communication
    void writeSummary( std::string filename );

    // The comm of the writer
    AMP_MPI d_comm;    

    // List of all meshes and thier ids
    std::map<AMP::Mesh::MeshID,siloBaseMeshData>            d_baseMeshes;
    std::map<AMP::Mesh::MeshID,AMP::Mesh::Mesh::shared_ptr> d_multiMeshes;

    // List of all variables
    std::set<std::string>   d_varNames;
};


}
}

#endif


