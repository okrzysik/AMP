#ifndef included_AMP_Vector_VectorFactory
#define included_AMP_Vector_VectorFactory

#include "AMP/utils/FactoryStrategy.hpp"


namespace AMP::IO {
class RestartManager;
}


namespace AMP::LinearAlgebra {

class Vector;

template<class VECTOR>
std::unique_ptr<VECTOR> createVectorFromRestart( int64_t fid, AMP::IO::RestartManager *manager )
{
    return std::make_unique<VECTOR>( fid, manager );
}


//! Vector factory class
class VectorFactory
{
public:
    //! get a singleton instance of the factory
    static VectorFactory &getFactory()
    {
        static VectorFactory singletonInstance;
        return singletonInstance;
    }

    //! Create the mesh from the restart file
    static std::shared_ptr<Vector> create( int64_t fid, AMP::IO::RestartManager *manager );

    //! Register a mesh with the factory
    template<class VECTOR>
    static void registerVector( const std::string &name )
    {
        FactoryStrategy<Vector, int64_t, AMP::IO::RestartManager *>::registerFactory(
            name, createVectorFromRestart<VECTOR> );
    }
};


} // namespace AMP::LinearAlgebra

#endif
