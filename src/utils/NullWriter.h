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
class NullWriter : public AMP::Utilities::Writer
{
public:
    //!  Default constructor
    NullWriter() {}

    //!  Default destructor
    virtual ~NullWriter() {}

    // Inherited functions
    virtual std::string getExtension() { return ""; }
    virtual void readFile( const std::string & ){};
    virtual void writeFile( const std::string &, size_t, double = 0 ) {}
#ifdef USE_AMP_MESH
    virtual void
    registerMesh( AMP::Mesh::Mesh::shared_ptr, int level = 1, std::string path = std::string() )
    {
        NULL_USE( level );
        NULL_USE( path );
    }
#endif
#ifdef USE_AMP_VECTORS
    virtual void registerVector( AMP::LinearAlgebra::Vector::shared_ptr,
                                 AMP::Mesh::Mesh::shared_ptr,
                                 AMP::Mesh::GeomType,
                                 const std::string &name = "" )
    {
        NULL_USE( name );
    }
    virtual void registerVector( AMP::LinearAlgebra::Vector::shared_ptr ) {}
#endif
#ifdef USE_AMP_MATRICES
    virtual void registerMatrix( AMP::LinearAlgebra::Matrix::shared_ptr ) {}
#endif
};
}
}

#endif
