#include "discretization/MultiDOF_Manager.h"

#include "ampmesh/MeshElementVectorIterator.h"
#include "utils/Utilities.h"


namespace AMP {
namespace Discretization {


/****************************************************************
* Constructors                                                  *
****************************************************************/
multiDOFManager::multiDOFManager(
    const AMP_MPI &globalComm, std::vector<DOFManager::shared_ptr> managers )
{
    d_managers = managers;
    d_comm     = globalComm;
    // Compute the total begin, end, and global size
    d_globalSize      = std::vector<size_t>( managers.size(), 0 );
    d_localSize       = std::vector<size_t>( managers.size(), 0 );
    size_t local_size = 0;
    for ( size_t i = 0; i < managers.size(); i++ ) {
        d_globalSize[i] = managers[i]->numGlobalDOF();
        d_localSize[i]  = managers[i]->numLocalDOF();
        local_size += d_localSize[i];
    }
    d_comm.sumScan( &local_size, &d_end, 1 );
    d_begin  = d_end - local_size;
    d_global = d_comm.bcast( d_end, d_comm.getSize() - 1 );
    // Compute the relationships between the DOFs
    d_subToGlobalDOF     = std::vector<std::vector<subDOF_struct>>( managers.size() );
    d_globalToSubDOF     = std::vector<std::vector<subDOF_struct>>( managers.size() );
    size_t manager_begin = d_begin;
    subDOF_struct mypair;
    for ( size_t i = 0; i < managers.size(); i++ ) {
        AMP_MPI comm = managers[i]->getComm();
        comm.barrier(); // Make sure all necessary processors are at the same point
        mypair.DOF1_begin   = managers[i]->beginDOF();
        mypair.DOF2_begin   = manager_begin;
        mypair.DOF1_end     = managers[i]->endDOF();
        mypair.DOF2_end     = manager_begin + d_localSize[i];
        d_subToGlobalDOF[i] = std::vector<subDOF_struct>( comm.getSize() );
        comm.allGather( mypair, &d_subToGlobalDOF[i][0] );
        AMP::Utilities::quicksort( d_subToGlobalDOF[i] );
        d_globalToSubDOF[i] = std::vector<subDOF_struct>( comm.getSize() );
        for ( size_t j = 0; j < d_subToGlobalDOF[i].size(); j++ ) {
            mypair.DOF1_begin      = d_subToGlobalDOF[i][j].DOF2_begin;
            mypair.DOF2_begin      = d_subToGlobalDOF[i][j].DOF1_begin;
            mypair.DOF1_end        = d_subToGlobalDOF[i][j].DOF2_end;
            mypair.DOF2_end        = d_subToGlobalDOF[i][j].DOF1_end;
            d_globalToSubDOF[i][j] = mypair;
        }
        AMP::Utilities::quicksort( d_globalToSubDOF[i] );
        manager_begin += d_localSize[i];
    }
}


/****************************************************************
* Deconstructor                                                 *
****************************************************************/
multiDOFManager::~multiDOFManager() {}


/****************************************************************
* Get the dofs for the element                                  *
****************************************************************/
void multiDOFManager::getDOFs( const AMP::Mesh::MeshElementID &id, std::vector<size_t> &dofs ) const
{
    dofs.resize( 0 );
    if ( d_managers.empty() )
        return;
    if ( d_managers[0]->numGlobalDOF() == this->d_global ) {
        // We are dealing with a multiDOFManager with only 1 sub DOF (this happens with
        // multivectors)
        d_managers[0]->getDOFs( id, dofs );
    } else {
        std::vector<size_t> local_dofs;
        for ( size_t i = 0; i < d_managers.size(); i++ ) {
            d_managers[i]->getDOFs( id, local_dofs );
            if ( !local_dofs.empty() ) {
                std::vector<size_t> tmp_dofs = getGlobalDOF( i, local_dofs );
                dofs.insert( dofs.end(), tmp_dofs.begin(), tmp_dofs.end() );
            }
        }
    }
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
    std::vector<AMP::shared_ptr<AMP::Mesh::MeshIterator>> iterators( d_managers.size() );
    for ( size_t i   = 0; i < d_managers.size(); i++ )
        iterators[i] = AMP::shared_ptr<AMP::Mesh::MeshIterator>(
            new AMP::Mesh::MeshIterator( d_managers[i]->getIterator() ) );
    // Get the list of unique elements
    size_t N_tot = 0;
    for ( size_t i = 0; i < iterators.size(); i++ )
        N_tot += iterators[i]->size();
    AMP::shared_ptr<std::vector<AMP::Mesh::MeshElement>> elements(
        new std::vector<AMP::Mesh::MeshElement>( 0 ) );
    elements->reserve( N_tot );
    for ( size_t i = 0; i < iterators.size(); i++ ) {
        AMP::Mesh::MeshIterator it2 = iterators[i]->begin();
        for ( size_t j = 0; j < it2.size(); j++ ) {
            elements->push_back( *it2 );
            ++it2;
        }
    }
    AMP::Utilities::unique( *elements );
    // Create an iterator over the elements
    return AMP::Mesh::MultiVectorIterator( elements );
}


/****************************************************************
* Return the remote DOFs for a vector                           *
****************************************************************/
std::vector<size_t> multiDOFManager::getRemoteDOFs() const
{
    std::vector<size_t> global_dofs;
    for ( size_t i = 0; i < d_managers.size(); i++ ) {
        std::vector<size_t> local_dofs = d_managers[i]->getRemoteDOFs();
        if ( !local_dofs.empty() ) {
            std::vector<size_t> tmp_dofs = getGlobalDOF( i, local_dofs );
            global_dofs.insert( global_dofs.end(), tmp_dofs.begin(), tmp_dofs.end() );
        }
    }
    AMP::Utilities::quicksort( global_dofs );
    return global_dofs;
}


/****************************************************************
* Return the global number of D.O.F.s                           *
****************************************************************/
std::vector<size_t> multiDOFManager::getRowDOFs( const AMP::Mesh::MeshElement &obj ) const
{
    std::vector<size_t> global_dofs;
    for ( size_t i = 0; i < d_managers.size(); i++ ) {
        std::vector<size_t> local_dofs = d_managers[i]->getRowDOFs( obj );
        if ( !local_dofs.empty() ) {
            std::vector<size_t> tmp_dofs = getGlobalDOF( i, local_dofs );
            global_dofs.insert( global_dofs.end(), tmp_dofs.begin(), tmp_dofs.end() );
        }
    }
    AMP::Utilities::quicksort( global_dofs );
    return global_dofs;
}


/****************************************************************
* Function to convert DOFs                                      *
****************************************************************/
std::vector<size_t> multiDOFManager::getGlobalDOF( const int i,
                                                   const std::vector<size_t> &subDOFs ) const
{
    if ( d_managers.size() == 1 ) {
        if ( d_managers[0]->numGlobalDOF() == this->d_global ) {
            // Special case where we only have 1 DOFManager
            return subDOFs;
        }
    }
    const size_t neg_one           = ~( (size_t) 0 );
    std::vector<size_t> globalDOFs = std::vector<size_t>( subDOFs.size(), neg_one );
    subDOF_struct search( 0, neg_one, neg_one, neg_one );
    for ( size_t j = 0; j < subDOFs.size(); j++ ) {
        search.DOF1_begin = subDOFs[j];
        size_t index      = AMP::Utilities::findfirst( d_subToGlobalDOF[i], search );
        index--;
        AMP_ASSERT( index < d_subToGlobalDOF[i].size() );
        globalDOFs[j] = subDOFs[j] - d_subToGlobalDOF[i][index].DOF1_begin +
                        d_subToGlobalDOF[i][index].DOF2_begin;
    }
    return globalDOFs;
}
std::vector<size_t> multiDOFManager::getSubDOF( const int i, std::vector<size_t> &globalDOFs ) const
{
    if ( d_managers.size() == 1 ) {
        if ( d_managers[0]->numGlobalDOF() == this->d_global ) {
            // Special case where we only have 1 DOFManager
            return globalDOFs;
        }
    }
    size_t neg_one = ~( (size_t) 0 );
    std::vector<size_t> subDOFs( globalDOFs.size(), neg_one );
    for ( size_t j = 0; j < globalDOFs.size(); j++ )
        subDOFs[j] = neg_one;
    subDOF_struct search( 0, neg_one, neg_one, neg_one );
    for ( size_t j = 0; j < globalDOFs.size(); j++ ) {
        search.DOF1_begin = globalDOFs[j];
        size_t index      = AMP::Utilities::findfirst( d_globalToSubDOF[i], search );
        index--;
        if ( index < d_globalToSubDOF[i].size() ) {
            if ( globalDOFs[j] >= d_globalToSubDOF[i][index].DOF1_begin &&
                 globalDOFs[j] < d_globalToSubDOF[i][index].DOF1_end )
                subDOFs[j] = globalDOFs[j] - d_globalToSubDOF[i][index].DOF1_begin +
                             d_globalToSubDOF[i][index].DOF2_begin;
        }
    }
    return subDOFs;
}


/****************************************************************
* Function to return the DOFManagers                            *
****************************************************************/
std::vector<DOFManager::shared_ptr> multiDOFManager::getDOFManagers() const { return d_managers; }


/****************************************************************
* Subset the DOF manager                                        *
****************************************************************/
AMP::shared_ptr<DOFManager> multiDOFManager::subset( const AMP_MPI &comm_in )
{
    // Check if we are dealing with a compatible comm
    if ( comm_in.compare( d_comm ) != 0 )
        return shared_from_this();
    // Get the comm for the new DOFManager
    AMP_MPI comm = AMP_MPI::intersect( comm_in, d_comm );
    // Subset all of the DOFManagers within this DOFManager
    std::vector<DOFManager::shared_ptr> sub_managers;
    for ( size_t i = 0; i < d_managers.size(); i++ ) {
        DOFManager::shared_ptr subset = d_managers[i]->subset( comm );
        if ( subset != nullptr )
            sub_managers.push_back( subset );
    }
    // Check that we have a valid DOF manager somewhere
    bool valid_DOF = !sub_managers.empty();
    valid_DOF      = comm.anyReduce( valid_DOF );
    if ( !valid_DOF )
        return AMP::shared_ptr<DOFManager>();
    // Create the new multiDOFManager
    return AMP::shared_ptr<DOFManager>( new multiDOFManager( comm, sub_managers ) );
}
AMP::shared_ptr<DOFManager> multiDOFManager::subset( const AMP::Mesh::Mesh::shared_ptr mesh,
                                                     bool useMeshComm )
{
    // Get the comm for the new DOFManager
    AMP_MPI comm( AMP_COMM_NULL );
    if ( useMeshComm ) {
        if ( mesh.get() != nullptr )
            comm = mesh->getComm();
    } else {
        comm = d_comm;
    }
    if ( comm.isNull() )
        return AMP::shared_ptr<DOFManager>();
    // Subset all of the DOFManagers within this DOFManager
    bool changed = false;
    std::vector<DOFManager::shared_ptr> sub_managers;
    for ( size_t i = 0; i < d_managers.size(); i++ ) {
        DOFManager::shared_ptr subset = d_managers[i]->subset( mesh, useMeshComm );
        if ( subset.get() != d_managers[i].get() )
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
        return AMP::shared_ptr<DOFManager>();
    // Create the new multiDOFManager
    return AMP::shared_ptr<DOFManager>( new multiDOFManager( comm, sub_managers ) );
}
AMP::shared_ptr<DOFManager> multiDOFManager::subset( const AMP::Mesh::MeshIterator &iterator,
                                                     const AMP_MPI &comm_in )
{
    // Get the comm for the new DOFManager
    AMP_MPI comm = AMP_MPI::intersect( comm_in, d_comm );
    // Subset all of the DOFManagers within this DOFManager
    bool changed = false;
    std::vector<DOFManager::shared_ptr> sub_managers;
    for ( size_t i = 0; i < d_managers.size(); i++ ) {
        DOFManager::shared_ptr subset = d_managers[i]->subset( iterator, comm );
        if ( subset.get() != d_managers[i].get() )
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
        return AMP::shared_ptr<DOFManager>();
    // Create the new multiDOFManager
    return AMP::shared_ptr<DOFManager>( new multiDOFManager( comm, sub_managers ) );
}
}
}
