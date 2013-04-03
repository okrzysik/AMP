#ifndef included_AMP_Writer
#define included_AMP_Writer

#include <string.h>
#include <sstream>
#include <vector>
#include <map>
#include <set>

#include "boost/smart_ptr/shared_ptr.hpp"

#include "utils/AMP_MPI.h"

#ifdef USE_AMP_MESH
    #include "ampmesh/Mesh.h"
#endif

#ifdef USE_AMP_VECTORS
    #include "vectors/Vector.h"
#endif


namespace AMP { 
namespace Utilities {


/**
 * \class Writer
 * \brief A class used to abstract away reading/writing files for visualization
 * \details  This class provides routines for reading, accessing and writing meshes and vectors.
 */
class Writer
{
public:

    //!  Convenience typedef
    typedef boost::shared_ptr<AMP::Utilities::Writer>  shared_ptr;

    /**
     * \brief   Function to build a writer
     * \details This function will build a default writer for use.
     * \param type   Writer type:
     *               "None" - An empty writer will be created
     *               "Silo" - A silo writer will be created if silo is configured, 
     *                        otherwise an empty writer will be created.
     */
    static boost::shared_ptr<AMP::Utilities::Writer> buildWriter( std::string type );

    //!  Default constructor
    Writer();

    //!  Default destructor
    virtual ~Writer();

    //!  Function to return the file extension
    virtual std::string getExtension()=0;

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
    virtual void  readFile( const std::string &fname )=0;

    //!  Function to write a file
    virtual void  writeFile( const std::string &fname, size_t iteration_count )=0;

#ifdef USE_AMP_MESH
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
    virtual void registerMesh( AMP::Mesh::Mesh::shared_ptr mesh, int level=1, std::string path=std::string() )=0;
#endif

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
    virtual void registerVector( AMP::LinearAlgebra::Vector::shared_ptr vec, AMP::Mesh::Mesh::shared_ptr mesh,
        AMP::Mesh::GeomType type, const std::string &name = "" )=0;
#endif

protected:

    // The comm of the writer
    AMP_MPI d_comm;

    // The decomposition to use
    int d_decomposition;


};


}
}

#endif


