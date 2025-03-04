#include "AMP/discretization/MultiDOF_Manager.h"
#include "AMP/mesh/MeshElementVectorIterator.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Utilities.h"


namespace AMP::Discretization {


/****************************************************************
 * Constructors                                                  *
 ****************************************************************/
multiDOFManager::multiDOFManager( const AMP_MPI &globalComm,
                                  std::vector<std::shared_ptr<DOFManager>> managers )
    : DOFManager(),
      d_managers( managers ),
      d_ids( managers.size(), 0 ),
      d_localSize( managers.size(), 0 ),
      d_globalSize( managers.size(), 0 )
{
    d_comm = globalComm;
    AMP_ASSERT( !d_comm.isNull() );
    // Compute the total begin, end, and global size
    size_t local_size = 0;
    for ( size_t i = 0; i < managers.size(); i++ ) {
        d_ids[i]        = managers[i]->getComm().rand();
        d_globalSize[i] = managers[i]->numGlobalDOF();
        d_localSize[i]  = managers[i]->numLocalDOF();
        local_size += d_localSize[i];
    }
    d_comm.sumScan( &local_size, &d_end, 1 );
    d_begin  = d_end - local_size;
    d_global = d_comm.bcast( d_end, d_comm.getSize() - 1 );
    // Compute the relationships between the DOFs
    d_dofMap.resize( managers.size() );
    size_t begin = d_begin;
    for ( size_t i = 0; i < managers.size(); i++ ) {
        d_dofMap[i] =
            DOFMapStruct( managers[i]->beginDOF(), managers[i]->endDOF(), begin, d_ids[i] );
        begin += managers[i]->numLocalDOF();
    }
    d_dofMap = d_comm.allGather( d_dofMap );
}


/****************************************************************
 * Deconstructor                                                 *
 ****************************************************************/
multiDOFManager::~multiDOFManager() = default;


/****************************************************************
 * Get the dofs for the element                                  *
 ****************************************************************/
size_t multiDOFManager::appendDOFs( const AMP::Mesh::MeshElementID &id,
                                    size_t *dofs,
                                    size_t index,
                                    size_t capacity ) const
{
    if ( d_managers.empty() )
        return 0;
    if ( d_managers[0]->numGlobalDOF() == this->d_global ) {
        // We are dealing with a multiDOFManager with only 1 sub DOF
        //   this happens with multivectors
        return d_managers[0]->appendDOFs( id, dofs, index, capacity );
    } else {
        size_t N = 0;
        for ( size_t i = 0; i < d_managers.size(); i++ ) {
            size_t N2 = d_managers[i]->appendDOFs( id, dofs, index, capacity );
            for ( size_t j = index; j < std::min( index + N2, capacity ); j++ )
                dofs[j] = subToGlobal( i, dofs[j] );
            N += N2;
            index += N2;
        }
        return N;
    }
}


/****************************************************************
 * Convert between local and global ids                          *
 ****************************************************************/
inline size_t multiDOFManager::subToGlobal( int manager, size_t dof ) const
{
    for ( const auto &map : d_dofMap ) {
        if ( map.inRangeLocal( dof ) && map.id() == d_ids[manager] )
            return map.toGlobal( dof );
    }
    return neg_one;
}
inline std::pair<size_t, int> multiDOFManager::globalToSub( size_t dof ) const
{
    for ( const auto &map : d_dofMap ) {
        if ( map.inRangeGlobal( dof ) ) {
            for ( size_t i = 0; i < d_managers.size(); i++ ) {
                if ( d_ids[i] == map.id() )
                    return std::make_pair( map.toLocal( dof ), i );
            }
        }
    }
    return std::make_pair( neg_one, -1 );
}


/****************************************************************
 * Get the element ID give a dof                                 *
 ****************************************************************/
AMP::Mesh::MeshElementID multiDOFManager::getElementID( size_t dof ) const
{
    auto map = globalToSub( dof );
    AMP_ASSERT( map.second >= 0 );
    return d_managers[map.second]->getElementID( map.first );
}
AMP::Mesh::MeshElement multiDOFManager::getElement( size_t dof ) const
{
    auto map = globalToSub( dof );
    AMP_ASSERT( map.second >= 0 );
    return d_managers[map.second]->getElement( map.first );
}


/****************************************************************
 * Get an entry over the mesh elements associated with the DOFs  *
 * Note: if any sub-DOFManagers are the same, then this will     *
 * iterate over repeated elements.                               *
 ****************************************************************/
AMP::Mesh::MeshIterator multiDOFManager::getIterator() const
{
    if ( d_managers.size() == 1 ) {
        return d_managers[0]->getIterator();
    }
    // Get the iterators for all sub DOFmanagers
    std::vector<std::shared_ptr<AMP::Mesh::MeshIterator>> iterators( d_managers.size() );
    for ( size_t i = 0; i < d_managers.size(); i++ )
        iterators[i] = std::make_shared<AMP::Mesh::MeshIterator>( d_managers[i]->getIterator() );
    // Get the list of unique elements
    size_t N_tot = 0;
    for ( auto &iterator : iterators )
        N_tot += iterator->size();
    auto elements = std::make_shared<std::vector<AMP::Mesh::MeshElement>>( 0 );
    elements->reserve( N_tot );
    for ( auto &iterator : iterators ) {
        for ( const auto &elem : *iterator )
            elements->push_back( elem );
    }
    AMP::Utilities::unique( *elements );
    // Create an iterator over the elements
    return AMP::Mesh::MeshElementVectorIterator( elements );
}


/****************************************************************
 * Return the remote DOFs for a vector                           *
 ****************************************************************/
std::vector<size_t> multiDOFManager::getRemoteDOFs() const
{
    std::vector<size_t> global_dofs;
    for ( size_t i = 0; i < d_managers.size(); i++ ) {
        auto local_dofs = d_managers[i]->getRemoteDOFs();
        global_dofs.reserve( global_dofs.size() + local_dofs.size() );
        for ( size_t j = 0; j < local_dofs.size(); j++ )
            global_dofs.push_back( subToGlobal( i, local_dofs[j] ) );
    }
    AMP::Utilities::quicksort( global_dofs );
    return global_dofs;
}


/****************************************************************
 * Return the global number of D.O.F.s                           *
 ****************************************************************/
size_t multiDOFManager::getRowDOFs( const AMP::Mesh::MeshElementID &id,
                                    size_t *dofs,
                                    size_t N_alloc,
                                    bool sort ) const
{
    size_t N = 0;
    for ( size_t i = 0; i < d_managers.size(); i++ ) {
        size_t N_alloc2 = N >= N_alloc ? 0 : N_alloc - N;
        size_t N_local  = d_managers[i]->getRowDOFs( id, &dofs[N], N_alloc2, false );
        size_t N2       = std::min( N_local, N_alloc2 );
        for ( size_t j = 0; j < N2; j++ )
            dofs[N + j] = subToGlobal( i, dofs[N + j] );
        N += N_local;
    }
    if ( sort )
        AMP::Utilities::quicksort( std::min( N, N_alloc ), dofs );
    return N;
}


/****************************************************************
 * Function to convert DOFs                                      *
 ****************************************************************/
std::vector<size_t> multiDOFManager::getGlobalDOF( const int manager,
                                                   const std::vector<size_t> &subDOFs ) const
{
    std::vector<size_t> dofs( subDOFs.size() );
    for ( size_t i = 0; i < dofs.size(); i++ )
        dofs[i] = subToGlobal( manager, subDOFs[i] );
    return dofs;
}
std::vector<size_t> multiDOFManager::getSubDOF( const int manager,
                                                const std::vector<size_t> &globalDOFs ) const
{
    std::vector<size_t> dofs( globalDOFs.size() );
    for ( size_t i = 0; i < dofs.size(); i++ ) {
        auto map = globalToSub( globalDOFs[i] );
        dofs[i]  = map.second == manager ? map.first : neg_one;
    }
    return dofs;
}


/****************************************************************
 * Function to return the DOFManagers                            *
 ****************************************************************/
std::vector<std::shared_ptr<DOFManager>> multiDOFManager::getDOFManagers() const
{
    return d_managers;
}


/****************************************************************
 * Subset the DOF manager                                        *
 ****************************************************************/
std::shared_ptr<DOFManager> multiDOFManager::subset( const AMP_MPI &comm_in )
{
    // Check if we are dealing with a compatible comm
    if ( comm_in.compare( d_comm ) != 0 )
        return shared_from_this();
    // Get the comm for the new DOFManager
    AMP_MPI comm = AMP_MPI::intersect( comm_in, d_comm );
    // Subset all of the DOFManagers within this DOFManager
    std::vector<std::shared_ptr<DOFManager>> sub_managers;
    for ( auto &elem : d_managers ) {
        auto subset = elem->subset( comm );
        if ( subset != nullptr )
            sub_managers.push_back( subset );
    }
    // Check that we have a valid DOF manager somewhere
    bool valid_DOF = !sub_managers.empty();
    valid_DOF      = comm.anyReduce( valid_DOF );
    if ( !valid_DOF )
        return std::shared_ptr<DOFManager>();
    // Create the new multiDOFManager
    return std::make_shared<multiDOFManager>( comm, sub_managers );
}
std::shared_ptr<DOFManager> multiDOFManager::subset( const std::shared_ptr<AMP::Mesh::Mesh> mesh,
                                                     bool useMeshComm )
{
    // Get the comm for the new DOFManager
    AMP_MPI comm( AMP_COMM_NULL );
    if ( useMeshComm ) {
        if ( mesh )
            comm = mesh->getComm();
    } else {
        comm = d_comm;
    }
    if ( comm.isNull() )
        return std::shared_ptr<DOFManager>();
    // Subset all of the DOFManagers within this DOFManager
    bool changed = false;
    std::vector<std::shared_ptr<DOFManager>> sub_managers;
    for ( auto &elem : d_managers ) {
        auto subset = elem->subset( mesh, useMeshComm );
        if ( subset.get() != elem.get() )
            changed = true;
        if ( subset != nullptr )
            sub_managers.push_back( subset );
    }
    // Check if the DOF manager changed and the comms are compatible
    if ( comm.compare( d_comm ) != 0 ) {
        changed = comm.anyReduce( changed );
        if ( !changed )
            return shared_from_this();
    }
    // Check that we have a valid DOF manager somewhere
    bool valid_DOF = !sub_managers.empty();
    valid_DOF      = comm.anyReduce( valid_DOF );
    if ( !valid_DOF )
        return std::shared_ptr<DOFManager>();
    // Create the new multiDOFManager
    return std::make_shared<multiDOFManager>( comm, sub_managers );
}
std::shared_ptr<DOFManager> multiDOFManager::subset( const AMP::Mesh::MeshIterator &iterator,
                                                     const AMP_MPI &comm_in )
{
    // Get the comm for the new DOFManager
    AMP_MPI comm = AMP_MPI::intersect( comm_in, d_comm );
    // Subset all of the DOFManagers within this DOFManager
    bool changed = false;
    std::vector<std::shared_ptr<DOFManager>> sub_managers;
    for ( auto &elem : d_managers ) {
        auto subset = elem->subset( iterator, comm );
        if ( subset.get() != elem.get() )
            changed = true;
        if ( subset != nullptr )
            sub_managers.push_back( subset );
    }
    // Check if the DOF manager changed and the comms are compatible
    if ( comm.compare( d_comm ) != 0 ) {
        changed = comm.anyReduce( changed );
        if ( !changed )
            return shared_from_this();
    }
    // Check that we have a valid DOF manager somewhere
    bool valid_DOF = !sub_managers.empty();
    valid_DOF      = comm.anyReduce( valid_DOF );
    if ( !valid_DOF )
        return std::shared_ptr<DOFManager>();
    // Create the new multiDOFManager
    return std::make_shared<multiDOFManager>( comm, sub_managers );
}
} // namespace AMP::Discretization


/********************************************************
 *  Explicit instantiations                              *
 ********************************************************/
#include "AMP/utils/AMP_MPI.I"
template std::vector<AMP::Discretization::multiDOFManager::DOFMapStruct>
AMP::AMP_MPI::allGather<AMP::Discretization::multiDOFManager::DOFMapStruct>(
    std::vector<AMP::Discretization::multiDOFManager::DOFMapStruct> const & ) const;
template std::vector<int> AMP::AMP_MPI::bcast<std::vector<int>>( std::vector<int> const &,
                                                                 int ) const;
template std::vector<std::array<double, 3ul>>
AMP::AMP_MPI::bcast<std::vector<std::array<double, 3ul>>>(
    std::vector<std::array<double, 3ul>> const &, int ) const;
