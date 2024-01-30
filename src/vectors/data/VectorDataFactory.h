#ifndef included_AMP_Vector_VectorDataFactory
#define included_AMP_Vector_VectorDataFactory

#include "AMP/utils/FactoryStrategy.hpp"


namespace AMP::IO {
class RestartManager;
}


namespace AMP::LinearAlgebra {

class VectorData;

template<class VECTORDATA>
std::unique_ptr<VECTORDATA> createVectorDataFromRestart( int64_t fid,
                                                         AMP::IO::RestartManager *manager )
{
    return std::make_unique<VECTORDATA>( fid, manager );
}


//! VectorData factory class
class VectorDataFactory
{
public:
    //! get a singleton instance of the factory
    static VectorDataFactory &getFactory()
    {
        static VectorDataFactory singletonInstance;
        return singletonInstance;
    }

    //! Create the vector data from the restart file
    static std::shared_ptr<VectorData> create( int64_t fid, AMP::IO::RestartManager *manager );

    //! Register vector data with the factory
    template<class VECTORDATA>
    static void registerVectorData( const std::string &name )
    {
        FactoryStrategy<VectorData, int64_t, AMP::IO::RestartManager *>::registerFactory(
            name, createVectorDataFromRestart<VECTORDATA> );
    }
};


} // namespace AMP::LinearAlgebra

#endif
