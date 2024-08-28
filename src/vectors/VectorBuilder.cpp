#include "AMP/vectors/VectorBuilder.h"
#include "AMP/discretization/MultiDOF_Manager.h"
#include "AMP/utils/memory.h"
#include "AMP/vectors/MultiVariable.h"
#include "AMP/vectors/MultiVector.h"
#include "AMP/vectors/VectorBuilder.hpp"

namespace AMP::LinearAlgebra {

/********************************************************
 *  Explicit instantiations                              *
 ********************************************************/

#define INSTANTIATE_VECTOR( TYPE )                                                                 \
    template Vector::shared_ptr                                                                    \
    createVector<TYPE, VectorOperationsDefault<TYPE>, VectorDataDefault<TYPE>>(                    \
        std::shared_ptr<AMP::Discretization::DOFManager>, std::shared_ptr<Variable>, bool );       \
    template Vector::shared_ptr createSimpleVector<TYPE>( size_t, const std::string & );           \
    template Vector::shared_ptr createSimpleVector<TYPE>( size_t, std::shared_ptr<Variable> );     \
    template Vector::shared_ptr createSimpleVector<TYPE>(                                          \
        size_t, std::shared_ptr<Variable>, AMP_MPI );                                              \
    template Vector::shared_ptr createSimpleVector<TYPE>(                                          \
        std::shared_ptr<Variable>,                                                                 \
        std::shared_ptr<AMP::Discretization::DOFManager>,                                          \
        std::shared_ptr<CommunicationList> );                                                      \
    template Vector::shared_ptr createArrayVector<TYPE>( const ArraySize &, const std::string & ); \
    template Vector::shared_ptr createArrayVector<TYPE>( const ArraySize &,                        \
                                                         std::shared_ptr<Variable> );              \
    template Vector::shared_ptr createArrayVector<TYPE>(                                           \
        const ArraySize &, const ArraySize &, const AMP_MPI &, std::shared_ptr<Variable> );        \
    template Vector::shared_ptr createVectorAdaptor<TYPE>(                                         \
        const std::string &, std::shared_ptr<AMP::Discretization::DOFManager>, TYPE * );           \
    template AMP::LinearAlgebra::Vector::shared_ptr createVector<TYPE>(                            \
        std::shared_ptr<AMP::Discretization::DOFManager> DOFs,                                     \
        std::shared_ptr<AMP::LinearAlgebra::Variable> variable,                                    \
        bool split,                                                                                \
        AMP::Utilities::MemoryType memType );
INSTANTIATE_VECTOR( double )
INSTANTIATE_VECTOR( float )

#ifdef USE_DEVICE
    #define INSTANTIATE_VECTOR_DEVICE( TYPE )                                                    \
        template Vector::shared_ptr                                                              \
            createSimpleVector<TYPE,                                                             \
                               VectorOperationsDevice<TYPE>,                                     \
                               VectorDataDefault<TYPE, AMP::ManagedAllocator<TYPE>>>(            \
                std::shared_ptr<Variable>,                                                       \
                std::shared_ptr<AMP::Discretization::DOFManager>,                                \
                std::shared_ptr<CommunicationList> );                                            \
        template Vector::shared_ptr                                                              \
            createSimpleVector<TYPE,                                                             \
                               VectorOperationsDevice<TYPE>,                                     \
                               VectorDataDefault<TYPE, AMP::DeviceAllocator<TYPE>>>(             \
                std::shared_ptr<Variable>,                                                       \
                std::shared_ptr<AMP::Discretization::DOFManager>,                                \
                std::shared_ptr<CommunicationList> );                                            \
        template Vector::shared_ptr                                                              \
        createVector<TYPE,                                                                       \
                     VectorOperationsDevice<TYPE>,                                               \
                     VectorDataDefault<TYPE, AMP::ManagedAllocator<TYPE>>>(                      \
            std::shared_ptr<AMP::Discretization::DOFManager>, std::shared_ptr<Variable>, bool ); \
        template Vector::shared_ptr                                                              \
        createVector<TYPE,                                                                       \
                     VectorOperationsDevice<TYPE>,                                               \
                     VectorDataDefault<TYPE, AMP::DeviceAllocator<TYPE>>>(                       \
            std::shared_ptr<AMP::Discretization::DOFManager>, std::shared_ptr<Variable>, bool );
INSTANTIATE_VECTOR_DEVICE( float )
INSTANTIATE_VECTOR_DEVICE( double )
#endif

} // namespace AMP::LinearAlgebra
