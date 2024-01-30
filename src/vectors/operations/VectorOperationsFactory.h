#ifndef included_AMP_Vector_VectorOperationsFactory
#define included_AMP_Vector_VectorOperationsFactory

#include "AMP/utils/FactoryStrategy.hpp"


namespace AMP::IO {
class RestartManager;
}


namespace AMP::LinearAlgebra {

class VectorOperations;

template<class VECTOROPERATIONS>
std::unique_ptr<VECTOROPERATIONS>
createVectorOperationsFromRestart( int64_t fid, AMP::IO::RestartManager *manager )
{
    return std::make_unique<VECTOROPERATIONS>( fid, manager );
}


//! Vector factory class
class VectorOperationsFactory
{
public:
    //! get a singleton instance of the factory
    static VectorOperationsFactory &getFactory()
    {
        static VectorOperationsFactory singletonInstance;
        return singletonInstance;
    }

    //! Create the vector from the restart file
    static std::shared_ptr<VectorOperations> create( int64_t fid,
                                                     AMP::IO::RestartManager *manager );

    //! Register a vector with the factory
    template<class VECTOROPERATIONS>
    static void registerVector( const std::string &name )
    {
        FactoryStrategy<VectorOperations, int64_t, AMP::IO::RestartManager *>::registerFactory(
            name, createVectorOperationsFromRestart<VECTOROPERATIONS> );
    }
};


} // namespace AMP::LinearAlgebra

#endif
