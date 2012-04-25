#include "ampmesh/MeshElementVectorIterator.h"
#include "discretization/DOF_Manager.h"
#include "discretization/subsetDOFManager.h"
#include "utils/Utilities.h"

#include <set>


namespace AMP {
namespace Discretization {


/****************************************************************
* Constructors                                                  *
****************************************************************/
DOFManager::DOFManager ( size_t N_local, AMP_MPI comm )
{
    d_comm = comm;
    d_comm.sumScan(&N_local,&d_end,1);
    d_begin = d_end - N_local;
    d_global = d_comm.bcast(d_end,d_comm.getSize()-1);
}


/****************************************************************
* Get the DOFs for the element                                  *
****************************************************************/
void DOFManager::getDOFs( const AMP::Mesh::MeshElementID &id, std::vector<size_t> &dofs ) const
{
    AMP_ERROR("getDOFs is not implimented for the base class");
}
void DOFManager::getDOFs( const std::vector<AMP::Mesh::MeshElementID> &ids, std::vector<size_t> &dofs ) const
{
    // This is a simple loop to provide a vector interface.  Ideally this should be overwritten by the user
    dofs.resize(0);
    dofs.reserve(2);
    std::vector<size_t> local_dofs;
    for (size_t i=0; i<ids.size(); i++) {
        getDOFs( ids[i], local_dofs );
        if ( local_dofs.size()+dofs.size() > dofs.capacity() )
            dofs.reserve(2*dofs.capacity());
        for (size_t j=0; j<local_dofs.size(); j++)
            dofs.push_back( local_dofs[j] );
    }
}


/****************************************************************
* Get an entry over the mesh elements associated with the DOFs  *
****************************************************************/
AMP::Mesh::MeshIterator DOFManager::getIterator( ) const
{
    return AMP::Mesh::MeshIterator();
}



/****************************************************************
* Return the first D.O.F. on this core                          *
****************************************************************/
size_t DOFManager::beginDOF( ) const
{
    return d_begin;
}


/****************************************************************
* Return the last D.O.F. on this core                           *
****************************************************************/
size_t DOFManager::endDOF( ) const
{
    return d_end;
}


/****************************************************************
* Return the local number of D.O.F.s                           *
****************************************************************/
size_t DOFManager::numLocalDOF( ) const
{
    return (d_end-d_begin);
}


/****************************************************************
* Return the global number of D.O.F.s                           *
****************************************************************/
size_t DOFManager::numGlobalDOF( ) const
{
    return d_global;
}


/****************************************************************
* Return the communicator                                       *
****************************************************************/
AMP_MPI DOFManager::getComm( ) const
{
    return d_comm;
}


/****************************************************************
* Return the global number of D.O.F.s                           *
****************************************************************/
std::vector<size_t> DOFManager::getRemoteDOFs( ) const
{
    return std::vector<size_t>();
}


/****************************************************************
* Return the global number of D.O.F.s                           *
****************************************************************/
std::vector<size_t> DOFManager::getRowDOFs( const AMP::Mesh::MeshElement &obj ) const
{
    AMP_ERROR("getRowDOFs is not implimented for the base class");
    return std::vector<size_t>();
}


/****************************************************************
* Compare two DOFManagers                                       *
****************************************************************/
bool DOFManager::operator==( const DOFManager &rhs ) const
{
    if ( this == &rhs ) {
        // We are pointing to the same object, it must be ==
        return true;
    }    
    if ( this->d_comm.compare(rhs.d_comm)==0 ) {
        // The comms do not match, the objects cannot be ==
        return false;
    }
    if ( this->numGlobalDOF()!=rhs.numGlobalDOF() || this->beginDOF()!=rhs.beginDOF() || this->endDOF()!=rhs.endDOF() ) {
        // The DOFs do not match, the objects cannot be ==
        return false;
    }    
    if ( this->getIterator()!=rhs.getIterator() ) {
        // The iterators do not match, the objects cannot be ==
        return false;
    }
    return true;
}
bool DOFManager::operator!=( const DOFManager &rhs ) const
{
    return !(this->operator==(rhs));
}



/****************************************************************
* Subset the DOF manager                                        *
****************************************************************/
boost::shared_ptr<DOFManager>  DOFManager::subset( AMP_MPI comm )
{
    if ( comm == d_comm ) 
        return shared_from_this();
    std::vector<size_t> local_dofs(numLocalDOF(),beginDOF());
    for (size_t i=0; i<numLocalDOF(); i++)
        local_dofs[i] += i;
    return subsetDOFManager::create( shared_from_this(), local_dofs, getIterator(), comm );
}
boost::shared_ptr<DOFManager>  DOFManager::subset( const AMP::Mesh::Mesh::shared_ptr mesh, bool useMeshComm )
{
    // Get a list of the elements in the mesh
    AMP::Mesh::MeshIterator iterator = getIterator();
    std::vector<AMP::Mesh::MeshID> meshIDs;
    if ( mesh.get() != NULL )
        meshIDs = mesh->getLocalBaseMeshIDs();
    std::set<AMP::Mesh::MeshElement> element_list;
    for (size_t i=0; i<iterator.size(); i++) {
        AMP::Mesh::MeshElement elem = *iterator;
        AMP::Mesh::MeshID meshID = elem.globalID().meshID();
        for (size_t j=0; j<meshIDs.size(); j++) {
            if ( meshID == meshIDs[j] ) {
                element_list.insert(elem);
                break;
            }
        }
        ++iterator;
    }
    // Create the element iterator
    boost::shared_ptr<std::vector<AMP::Mesh::MeshElement> > elements( 
        new std::vector<AMP::Mesh::MeshElement>(element_list.begin(),element_list.end()) );
    AMP::Mesh::MeshIterator  subsetIterator = AMP::Mesh::MultiVectorIterator( elements, 0 );
    // Get the DOFs
    std::vector<AMP::Mesh::MeshElementID> id_list(elements->size());
    for (size_t i=0; i<elements->size(); i++)
        id_list[i] = elements->operator[](i).globalID();
    std::vector<size_t> dofs;
    getDOFs( id_list, dofs );
    // Sort and check the DOFs for errors
    AMP::Utilities::quicksort(dofs);
    for (size_t i=0; i<dofs.size(); i++) {
        if ( dofs[i]<d_begin || dofs[i]>=d_end )
            AMP_ERROR("Internal error subsetting DOF manager (out of range)");
    }
    for (size_t i=1; i<dofs.size(); i++) {
        if ( dofs[i]==dofs[i-1] )
            AMP_ERROR("Internal error subsetting DOF manager (duplicate)");
    }
    // Create the subset DOF Manager
    AMP_MPI comm(AMP_COMM_NULL);
    if ( useMeshComm ) {
        if ( mesh.get()!=NULL ) 
            comm = mesh->getComm();
    } else {
        comm = d_comm;
    }
    if ( comm.isNull() )
        return boost::shared_ptr<DOFManager>();
    return subsetDOFManager::create( shared_from_this(), dofs, subsetIterator, comm );
}
boost::shared_ptr<DOFManager>  DOFManager::subset( const AMP::Mesh::MeshIterator &iterator, AMP_MPI comm )
{
    // Get the intesection of the current iterator with the given iterator
    AMP::Mesh::MeshIterator intersection = AMP::Mesh::Mesh::getIterator( AMP::Mesh::Intersection, iterator, getIterator() );
    if ( intersection.size()==getIterator().size() )
        intersection = getIterator();
    // Get the list of element we want
    std::vector<AMP::Mesh::MeshElementID> element_list(iterator.size());
    for (size_t i=0; i<intersection.size(); i++) {
        element_list[i] = intersection->globalID();
        ++intersection;
    }
    // Get the DOFs
    std::vector<size_t> dofs;
    getDOFs( element_list, dofs );
    // Sort and check the DOFs for errors
    AMP::Utilities::quicksort(dofs);
    for (size_t i=0; i<dofs.size(); i++) {
        if ( dofs[i]<d_begin || dofs[i]>=d_end )
            AMP_ERROR("Internal error subsetting DOF manager (out of range)");
    }
    for (size_t i=1; i<dofs.size(); i++) {
        if ( dofs[i]==dofs[i-1] )
            AMP_ERROR("Internal error subsetting DOF manager (duplicate)");
    }
    // Create the subset DOF Manager    
    if ( comm.isNull() )
        return boost::shared_ptr<DOFManager>();
    return subsetDOFManager::create( shared_from_this(), dofs, intersection, comm );
}


}
}

