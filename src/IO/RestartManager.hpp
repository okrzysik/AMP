#ifndef included_AMP_RestartManager_hpp
#define included_AMP_RestartManager_hpp

#include "AMP/IO/RestartManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/UtilityMacros.h"


namespace AMP::Geometry {
class Geometry;
}
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
void RestartManager::registerData( const TYPE &data, const std::string &name )
{
    using AMP::LinearAlgebra::Matrix;
    using AMP::LinearAlgebra::Vector;
    using AMP::Mesh::Mesh;
    using AMP::Solver::SolverStrategy;
    using AMP::TimeIntegrator::TimeIntegrator;
    if constexpr ( is_string_v<TYPE> ) {
        registerData( std::make_shared<const std::string>( data ), name );
    } else if constexpr ( !AMP::is_shared_ptr_v<TYPE> ) {
        registerData( std::make_shared<const TYPE>( data ), name );
    } else if constexpr ( !std::is_const_v<typename TYPE::element_type> ) {
        registerData( std::const_pointer_cast<const typename TYPE::element_type>( data ), name );
    } else {
        std::shared_ptr<const DataStore> obj;
        using TYPE2 = remove_cvref_t<typename TYPE::element_type>;
        if constexpr ( std::is_base_of_v<DataStore, TYPE2> || std::is_same_v<TYPE2, DataStore> ) {
            obj = data;
        } else if constexpr ( std::is_base_of_v<Mesh, TYPE2> && !std::is_same_v<TYPE2, Mesh> ) {
            obj = create( name, std::dynamic_pointer_cast<const Mesh>( data ) );
        } else if constexpr ( std::is_base_of_v<Vector, TYPE2> && !std::is_same_v<TYPE2, Vector> ) {
            obj = create( name, std::dynamic_pointer_cast<const Vector>( data ) );
        } else if constexpr ( std::is_base_of_v<Matrix, TYPE2> && !std::is_same_v<TYPE2, Matrix> ) {
            obj = create( name, std::dynamic_pointer_cast<const Matrix>( data ) );
        } else if constexpr ( std::is_base_of_v<SolverStrategy, TYPE2> &&
                              !std::is_same_v<TYPE2, SolverStrategy> ) {
            obj = create( name, std::dynamic_pointer_cast<const SolverStrategy>( data ) );
        } else if constexpr ( std::is_base_of_v<TimeIntegrator, TYPE2> &&
                              !std::is_same_v<TYPE2, TimeIntegrator> ) {
            obj = create( name, std::dynamic_pointer_cast<const TimeIntegrator>( data ) );
        } else {
            obj = create<TYPE2>( name, data );
        }
        auto hash = obj->getHash();
        AMP_ASSERT( hash != 0 );
        d_data[hash] = obj;
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
RestartManager::DataStorePtr RestartManager::create( const std::string &name,
                                                     std::shared_ptr<const TYPE> data )
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
