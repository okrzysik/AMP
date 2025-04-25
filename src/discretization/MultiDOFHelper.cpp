#include "AMP/discretization/MultiDOFHelper.h"
#include "AMP/IO/HDF5.h"
#include "AMP/vectors/data/VectorData.h"


namespace AMP::Discretization {


/****************************************************************
 * Constructor                                                   *
 ****************************************************************/
multiDOFHelper::multiDOFHelper( const std::vector<std::shared_ptr<DOFManager>> &managers,
                                const AMP::AMP_MPI &comm )
{
    PROFILE( "multiDOFHelper" );
    d_ids.resize( managers.size() );
    d_dofMap.resize( managers.size() );
    d_N_local = 0;
    for ( auto &DOF : managers )
        d_N_local += DOF->numLocalDOF();
    d_local = comm.allGather( d_N_local );
    d_begin = 0;
    for ( int i = 0; i < comm.getRank(); i++ )
        d_begin += d_local[i];
    d_N_global = 0;
    for ( int i = 0; i < comm.getSize(); i++ )
        d_N_global += d_local[i];
    size_t globalBegin = d_begin;
    for ( size_t i = 0; i < managers.size(); i++ ) {
        size_t begin = managers[i]->beginDOF();
        size_t end   = managers[i]->endDOF();
        d_ids[i]     = managers[i]->getComm().rand();
        d_dofMap[i]  = DOFMapStruct( begin, end, globalBegin, d_ids[i] );
        globalBegin += managers[i]->numLocalDOF();
    }
    d_dofMap = comm.allGather( d_dofMap );
}
multiDOFHelper::multiDOFHelper( const std::vector<AMP::LinearAlgebra::VectorData *> &data,
                                const AMP::AMP_MPI &comm )
{
    PROFILE( "multiDOFHelper" );
    d_ids.resize( data.size() );
    d_dofMap.resize( data.size() );
    d_N_local = 0;
    for ( auto &vec : data )
        d_N_local += vec->getLocalSize();
    d_local = comm.allGather( d_N_local );
    for ( int i = 0; i < comm.getRank(); i++ )
        d_begin += d_local[i];
    d_N_global = 0;
    for ( int i = 0; i < comm.getSize(); i++ )
        d_N_global += d_local[i];
    size_t globalBegin = d_begin;
    for ( size_t i = 0; i < data.size(); i++ ) {
        size_t begin = data[i]->getLocalStartID();
        size_t end   = begin + data[i]->getLocalSize();
        d_ids[i]     = data[i]->getID();
        d_dofMap[i]  = DOFMapStruct( begin, end, globalBegin, d_ids[i] );
        globalBegin += data[i]->getLocalSize();
    }
    d_dofMap = comm.allGather( d_dofMap );
}


/****************************************************************
 * Convert between local and global ids                          *
 ****************************************************************/
size_t multiDOFHelper::subToGlobal( int manager, size_t dof ) const
{
    constexpr size_t neg_one = ~( (size_t) 0 );
    for ( const auto &map : d_dofMap ) {
        if ( map.inRangeLocal( dof ) && map.id() == d_ids[manager] )
            return map.toGlobal( dof );
    }
    return neg_one;
}
std::pair<size_t, int> multiDOFHelper::globalToSub( size_t dof ) const
{
    constexpr size_t neg_one = ~( (size_t) 0 );
    for ( const auto &map : d_dofMap ) {
        if ( map.inRangeGlobal( dof ) ) {
            for ( size_t i = 0; i < d_ids.size(); i++ ) {
                if ( d_ids[i] == map.id() )
                    return std::make_pair( map.toLocal( dof ), i );
            }
        }
    }
    return std::make_pair( neg_one, -1 );
}
std::vector<size_t> multiDOFHelper::getGlobalDOF( const int manager,
                                                  const std::vector<size_t> &subDOFs ) const
{
    std::vector<size_t> dofs( subDOFs.size() );
    for ( size_t i = 0; i < dofs.size(); i++ )
        dofs[i] = subToGlobal( manager, subDOFs[i] );
    return dofs;
}
std::vector<size_t> multiDOFHelper::getSubDOF( const int manager,
                                               const std::vector<size_t> &globalDOFs ) const
{
    constexpr size_t neg_one = ~( (size_t) 0 );
    std::vector<size_t> dofs( globalDOFs.size() );
    for ( size_t i = 0; i < dofs.size(); i++ ) {
        auto map = globalToSub( globalDOFs[i] );
        dofs[i]  = map.second == manager ? map.first : neg_one;
    }
    return dofs;
}


/********************************************************
 *  Write/Read from HDF5                                 *
 ********************************************************/
void multiDOFHelper::writeHDF5( size_t gid ) const
{
    AMP::IO::writeHDF5( gid, "N_local", d_N_local );
    AMP::IO::writeHDF5( gid, "N_global", d_N_global );
    AMP::IO::writeHDF5( gid, "begin", d_begin );
    AMP::IO::writeHDF5( gid, "ids", d_ids );
    AMP::IO::writeHDF5( gid, "map", d_dofMap );
}
multiDOFHelper::multiDOFHelper( size_t gid )
{
    AMP::IO::readHDF5( gid, "N_local", d_N_local );
    AMP::IO::readHDF5( gid, "N_global", d_N_global );
    AMP::IO::readHDF5( gid, "begin", d_begin );
    AMP::IO::readHDF5( gid, "ids", d_ids );
    AMP::IO::readHDF5( gid, "map", d_dofMap );
}


} // namespace AMP::Discretization


/********************************************************
 *  Write/Read from HDF5                                 *
 ********************************************************/
#include "AMP/IO/HDF5.hpp"
template<>
void AMP::IO::writeHDF5<AMP::Discretization::multiDOFHelper>(
    hid_t fid, const std::string &name, const AMP::Discretization::multiDOFHelper &data )
{
    auto gid = AMP::IO::createGroup( fid, name );
    data.writeHDF5( gid );
    AMP::IO::closeGroup( fid );
}
template<>
void AMP::IO::readHDF5<AMP::Discretization::multiDOFHelper>(
    hid_t fid, const std::string &name, AMP::Discretization::multiDOFHelper &data )
{
    auto gid = AMP::IO::openGroup( fid, name );
    data     = AMP::Discretization::multiDOFHelper( gid );
    AMP::IO::closeGroup( fid );
}
template<>
void AMP::IO::writeHDF5Array<AMP::Discretization::multiDOFHelper::DOFMapStruct>(
    hid_t fid,
    const std::string &name,
    const AMP::Array<AMP::Discretization::multiDOFHelper::DOFMapStruct> &data )
{
    size_t bytes = data.length() * sizeof( AMP::Discretization::multiDOFHelper::DOFMapStruct );
    AMP::Array<size_t> data2( cat( ArraySize( 4 ), data.size() ) );
    AMP_ASSERT( data2.length() * sizeof( size_t ) == bytes );
    memcpy( data2.data(), data.data(), bytes );
    writeHDF5Array( fid, name, data2 );
}
template<>
void AMP::IO::readHDF5Array<AMP::Discretization::multiDOFHelper::DOFMapStruct>(
    hid_t fid,
    const std::string &name,
    AMP::Array<AMP::Discretization::multiDOFHelper::DOFMapStruct> &data )
{
    AMP::Array<size_t> data2;
    readHDF5Array( fid, name, data2 );
    AMP_ASSERT( data2.size( 0 ) == 4 );
    data.resize( pop( data2.size() ) );
    size_t bytes = data.length() * sizeof( AMP::Discretization::multiDOFHelper::DOFMapStruct );
    AMP_ASSERT( data2.length() * sizeof( size_t ) == bytes );
    memcpy( (void *) data.data(), data2.data(), bytes );
}


/********************************************************
 *  Explicit instantiations                              *
 ********************************************************/
#include "AMP/utils/AMP_MPI.I"
#include "AMP/utils/Array.hpp"
instantiateArrayConstructors( AMP::Discretization::multiDOFHelper::DOFMapStruct );
template std::vector<AMP::Discretization::multiDOFHelper::DOFMapStruct>
AMP::AMP_MPI::allGather<AMP::Discretization::multiDOFHelper::DOFMapStruct>(
    std::vector<AMP::Discretization::multiDOFHelper::DOFMapStruct> const & ) const;
