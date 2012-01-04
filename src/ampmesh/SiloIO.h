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
    //!  Function to register a vector
    void registerVector( AMP::LinearAlgebra::Vector::shared_ptr vec, const std::string &s = "" );
#endif

private:
    
    // The comm of the writer
    AMP_MPI d_comm;    

    // List of all meshes and thier ids
    std::map<AMP::Mesh::MeshID,AMP::Mesh::Mesh::shared_ptr>  d_baseMeshes;
    std::map<AMP::Mesh::MeshID,AMP::Mesh::Mesh::shared_ptr>  d_multiMeshes;

    // Function to write a single mesh
    std::string  writeMesh( DBfile *file, AMP::Mesh::Mesh::shared_ptr mesh );

    // Structure used to hold data for the silo multimeshes
    struct siloMultiMeshData {
        AMP::Mesh::MeshID           id;         // Unique ID to identify the mesh
        std::string                 name;       // Name of the multimesh in silo
        std::vector<std::string>    paths;      // Paths to the meshes
        // Function to count the number of bytes needed to pack the data
        size_t size();
        // Function to pack the data to a byte array
        void pack( char* );
        // Function to unpack the data from a byte array
        static siloMultiMeshData unpack( char* );
    };
    
    // Function to syncronize multimesh data across all processors
    void syncMultiMeshData( std::map<AMP::Mesh::MeshID,siloMultiMeshData> &data );
};


}
}

#endif


