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
    d_local_id.resize(mesh->numLocalElements(AMP::Mesh::Vertex));
    AMP::Mesh::MeshIterator pos = mesh->getIterator(type,0);
    for (size_t i=0; i<d_local_id.size(); i++) {
        d_local_id[i] = pos->globalID();
        pos++;
    }
    AMP::Utilities::quicksort(d_local_id);
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
    if ( type != AMP::Mesh::Vertex )
        AMP_ERROR("Not yet programmed for elements other than nodes");
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
    mvparams->d_Engine = AMP::LinearAlgebra::VectorEngine::shared_ptr ( new AMP::LinearAlgebra::EpetraVectorEngine( eveparams, t_buffer ) );
    //mvparams->d_CommList = d_vDOFMapCache[var->variableID()]->getCommunicationList();
    AMP_ERROR("Not finished yet");
    AMP::LinearAlgebra::Vector::shared_ptr vector = AMP::LinearAlgebra::Vector::shared_ptr( new AMP::LinearAlgebra::ManagedPetscVector(mvparams) );
    
    return vector;
}


}
}

