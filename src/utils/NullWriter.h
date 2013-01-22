#ifndef included_AMP_NullWriter
#define included_AMP_NullWriter

#include "utils/Writer.h"


namespace AMP { 
namespace Utilities {


/**
 * \class NullWriter
 * \brief A class used to abstract away reading/writing files for visualization
 * \details  This class provides default routines that do nothing.
 */
class NullWriter: public AMP::Utilities::Writer
{
public:

    //!  Default constructor
    NullWriter() {}

    //!  Default destructor
    virtual ~NullWriter() {}

    // Inherited functions
    virtual std::string getExtension() { return ""; }
    virtual void  readFile( const std::string &fname ) {};
    virtual void  writeFile( const std::string &fname, size_t iteration_count ) {}
    #ifdef USE_AMP_MESH
        virtual void registerMesh( AMP::Mesh::Mesh::shared_ptr mesh, int level=1, std::string path=std::string() ) {}
    #endif
    #ifdef USE_AMP_VECTORS
        virtual void registerVector( AMP::LinearAlgebra::Vector::shared_ptr vec, AMP::Mesh::Mesh::shared_ptr mesh,
            AMP::Mesh::GeomType type, const std::string &name = "" ) {}
    #endif
};


}
}

#endif


