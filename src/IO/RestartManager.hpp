#ifndef included_AMP_RestartManager_hpp
#define included_AMP_RestartManager_hpp

#include "AMP/IO/RestartManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/UtilityMacros.h"


namespace AMP::Mesh {
class Mesh;
}
namespace AMP::LinearAlgebra {
class Vector;
class Matrix;
} // namespace AMP::LinearAlgebra
namespace AMP::Operator {
class Operator;
}
namespace AMP::Solver {
class SolverStrategy;
}
namespace AMP::TimeIntegrator {
class TimeIntegrator;
}


namespace AMP::IO {


/********************************************************
 *  Register data with the manager                       *
 ********************************************************/
template<class TYPE>
void RestartManager::registerData( const std::string &name, const TYPE &data )
{
    using AMP::LinearAlgebra::Matrix;
    using AMP::LinearAlgebra::Vector;
    using AMP::Mesh::Mesh;
    using AMP::Solver::SolverStrategy;
    using AMP::TimeIntegrator::TimeIntegrator;
    if constexpr ( is_string_v<TYPE> ) {
        registerData( name, std::make_shared<const std::string>( data ) );
    } else if constexpr ( !AMP::is_shared_ptr<TYPE>::value ) {
        registerData( name, std::make_shared<const TYPE>( data ) );
    } else if constexpr ( !std::is_const_v<typename TYPE::element_type> ) {
        registerData( name, std::const_pointer_cast<const typename TYPE::element_type>( data ) );
    } else {
        std::shared_ptr<DataStore> obj;
        typedef typename TYPE::element_type TYPE2;
        if constexpr ( std::is_base_of_v<Mesh, TYPE2> &&
                       !std::is_same_v<TYPE, std::shared_ptr<Mesh>> ) {
            obj = create( name, std::dynamic_pointer_cast<const Mesh>( data ) );
        } else if constexpr ( std::is_base_of_v<Vector, TYPE2> &&
                              !std::is_same_v<TYPE, std::shared_ptr<Vector>> ) {
            obj = create( name, std::dynamic_pointer_cast<const Vector>( data ) );
        } else if constexpr ( std::is_base_of_v<Matrix, TYPE2> &&
                              !std::is_same_v<TYPE, std::shared_ptr<Matrix>> ) {
            obj = create( name, std::dynamic_pointer_cast<const Matrix>( data ) );
        } else if constexpr ( std::is_base_of_v<SolverStrategy, TYPE2> &&
                              !std::is_same_v<TYPE, std::shared_ptr<SolverStrategy>> ) {
            obj = create( name, std::dynamic_pointer_cast<const SolverStrategy>( data ) );
        } else if constexpr ( std::is_base_of_v<TimeIntegrator, TYPE2> &&
                              !std::is_same_v<TYPE, std::shared_ptr<TimeIntegrator>> ) {
            obj = create( name, std::dynamic_pointer_cast<const TimeIntegrator>( data ) );
        } else {
            obj = create( name, data );
        }
        registerData( obj );
        auto hash = obj->getHash();
        if ( !name.empty() ) {
            auto it = d_names.find( name );
            if ( it == d_names.end() )
                d_names[name] = hash;
            else if ( it->second != hash )
                AMP_ERROR( "Two objects registered with the same name" );
        }
    }
}
template<class TYPE>
std::shared_ptr<RestartManager::DataStore>
RestartManager::create( const std::string &name, std::shared_ptr<const TYPE> data )
{
    return std::make_shared<DataStoreType<TYPE>>( name, data, this );
}
template<class TYPE>
std::shared_ptr<TYPE> RestartManager::getData( const std::string &name )
{
    auto it = d_names.find( name );
    AMP_INSIST( it != d_names.end(), "Object not found: " + name );
    auto hash = it->second;
    return getData<TYPE>( hash );
}


} // namespace AMP::IO


#endif
