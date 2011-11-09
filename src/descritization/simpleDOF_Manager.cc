#include "descritization/simpleDOF_Manager.h"
#include "utils/Utilities.h"

#include "vectors/petsc/ManagedPetscVector.h"
#include "vectors/trilinos/EpetraVectorEngine.h"

namespace AMP {
namespace Discretization {


/****************************************************************
* Constructors                                                  *
****************************************************************/
simpleDOFManager::simpleDOFManager( boost::shared_ptr<AMP::Mesh::Mesh> mesh, AMP::Mesh::GeomType type, int gcw )
{
    d_mesh = mesh;
    d_type = type;
    d_gcw = gcw;
    // Create a sorted list of the local and remote types
    d_local_id.resize(mesh->numLocalElements(type));
    d_remote_id.resize(mesh->numGhostElements(type,gcw));
    AMP::Mesh::MeshIterator pos = mesh->getIterator(type,gcw);
    AMP::Mesh::MeshIterator end = pos.end();
    int i=0;
    int j=0;
    while ( pos != end ) {
        AMP::Mesh::MeshElementID id = pos->globalID();
        if ( id.is_local ) {
            d_local_id[i] = pos->globalID();
            i++;
        } else {
            d_remote_id[j] = pos->globalID();
            j++;
        }
        pos++;
    }
    AMP::Utilities::quicksort(d_local_id);
    AMP::Utilities::quicksort(d_remote_id);
}


/****************************************************************
* Get the entry indices of nodal values given a mesh element    *
* Note:  this function is likely temporary, we are assuming     *
* all data will be stored on nodes.                             *
****************************************************************/
void simpleDOFManager::getDOFs( const AMP::Mesh::MeshElement &obj, std::vector <unsigned int> &ids, unsigned int which ) const
{
    std::vector<AMP::Mesh::MeshElement> elements = obj.getElements(AMP::Mesh::Vertex);
    if ( which==static_cast<unsigned int>(-1) ) {
        // Return all dofs
        ids.resize(elements.size());
        for (size_t i=0; i<elements.size(); i++) {
            AMP::Mesh::MeshElementID local_id = elements[i].globalID();
            ids[i] = AMP::Utilities::findfirst(d_local_id,local_id);
            AMP_INSIST(local_id==d_local_id[ids[i]],"Internal Error: id not found");
        }
    } else {
        // Return only the desired dof
        ids.resize(1);
        AMP::Mesh::MeshElementID local_id = elements[which].globalID();
        ids[0] = AMP::Utilities::findfirst(d_local_id,local_id);
        AMP_INSIST(local_id==d_local_id[ids[0]],"Internal Error: id not found");
    }
}


/****************************************************************
* Retrun the first D.O.F. on this core                          *
****************************************************************/
size_t simpleDOFManager::beginDOF( )
{
    AMP_ERROR("beginDOF is not implimented yet");
    return 0;
}


/****************************************************************
* Retrun the last D.O.F. on this core                           *
****************************************************************/
size_t simpleDOFManager::endDOF( )
{
    AMP_ERROR("endDOF is not implimented yet");
    return 0;
}


/****************************************************************
* Retrun the first D.O.F. on this core                          *
****************************************************************/
AMP::LinearAlgebra::Vector::shared_ptr simpleDOFManager::createVector( AMP::LinearAlgebra::Variable::shared_ptr variable )
{
    AMP::Mesh::GeomType type = (AMP::Mesh::GeomType) variable->variableID();
    if ( type != d_type )
        AMP_ERROR("The variableID must match the element type specified at construction");
    // Create the communication list
    AMP::LinearAlgebra::CommunicationList::shared_ptr comm_list;
    // Create the vector parameters
    boost::shared_ptr<AMP::LinearAlgebra::ManagedPetscVectorParameters> mvparams(
        new AMP::LinearAlgebra::ManagedPetscVectorParameters() );
    boost::shared_ptr<AMP::LinearAlgebra::EpetraVectorEngineParameters> eveparams(
        new AMP::LinearAlgebra::EpetraVectorEngineParameters( d_mesh->numLocalElements(type), d_mesh->numGlobalElements(type), d_mesh->getComm() ) );
    int i = 0;
    for (size_t local_start=this->beginDOF(); local_start!=this->endDOF(); local_start++, i++ ) {
        eveparams->addMapping ( i , local_start );
    }
    AMP::LinearAlgebra::VectorEngine::BufferPtr t_buffer ( new AMP::LinearAlgebra::VectorEngine::Buffer ( d_mesh->numLocalElements(type) ) );
    AMP::LinearAlgebra::VectorEngine::shared_ptr epetra_engine( new AMP::LinearAlgebra::EpetraVectorEngine( eveparams, t_buffer ) );
    mvparams->d_Engine = epetra_engine;
    mvparams->d_CommList = comm_list;
    // Create the vector
    AMP::LinearAlgebra::Vector::shared_ptr vector = AMP::LinearAlgebra::Vector::shared_ptr( new AMP::LinearAlgebra::ManagedPetscVector(mvparams) );
    return vector;
}


}
}

